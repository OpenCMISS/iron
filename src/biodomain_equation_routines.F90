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

!>This module handles all bioelectric domain equation routines.
MODULE BiodomainEquationsRoutines

  USE BaseRoutines
  USE BasisRoutines
  USE BasisAccessRoutines
  USE BoundaryConditionsRoutines
  USE Constants
  USE ControlLoopRoutines
  USE ControlLoopAccessRoutines
  USE CoordinateSystemRoutines
  USE DecompositionAccessRoutines
  USE DistributedMatrixVector
  USE DistributedMatrixVectorAccessRoutines
  USE DomainMappings
  USE ElectrophysiologyCellRoutines
  USE EquationsRoutines
  USE EquationsAccessRoutines
  USE EquationsMappingRoutines
  USE EquationsMappingAccessRoutines
  USE EquationsMatricesRoutines
  USE EquationsMatricesAccessRoutines
  USE EquationsSetAccessRoutines
  USE FieldRoutines
  USE FieldAccessRoutines
  USE FIELD_IO_ROUTINES
  USE InputOutput
  USE ISO_VARYING_STRING
  USE Kinds
  USE MatrixVector
  USE ProblemAccessRoutines
  USE RegionAccessRoutines
  USE Strings
  USE SolverRoutines
  USE SolverAccessRoutines
  USE SolverMappingAccessRoutines
  USE Timer
  USE Types

#include "macros.h"

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC Biodomain_PostLoop
  
  PUBLIC Biodomain_EquationsSetSetup

  PUBLIC Biodomain_EquationsSetSolutionMethodSet
  
  PUBLIC Biodomain_EquationsSetSpecificationSet

  PUBLIC Biodomain_FiniteElementCalculate

  PUBLIC Biodomain_PostSolve
  
  PUBLIC Biodomain_PreSolve
  
  PUBLIC Biodomain_ProblemSetup

  PUBLIC Biodomain_ProblemSpecificationSet
  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Runs after each control loop iteration
  SUBROUTINE Biodomain_PostLoop(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: currentIteration,equationsSetIdx,inputIteration,loopType,numberOfEquationsSets,outputIteration, &
      & parentIteration,regionUserNumber
    REAL(DP) :: currentTime,startTime,stopTime,timeIncrement
    TYPE(ControlLoopType), POINTER :: parentLoop
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldsType), POINTER :: fields
    TYPE(RegionType), POINTER :: region   
    TYPE(SolverType), POINTER :: solver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: filename,localError,method

    ENTERS("Biodomain_PostLoop",err,error,*999)

    CALL ControlLoop_TypeGet(controlLoop,loopType,err,error,*999)
    SELECT CASE(loopType)
    CASE(CONTROL_SIMPLE_TYPE)
      !do nothing
    CASE(CONTROL_FIXED_LOOP_TYPE)
      !do nothing
    CASE(CONTROL_TIME_LOOP_TYPE)
      !Export the dependent field for this time step
      CALL ControlLoop_CurrentTimeInformationGet(controlLoop,currentTime,timeIncrement,startTime,stopTime,currentIteration, &
        & outputIteration,inputIteration,err,error,*999)
      IF(outputIteration>0) THEN
        IF(MOD(currentIteration,outputIteration)==0) THEN
          NULLIFY(parentLoop)
          CALL ControlLoop_ParentLoopExists(controlLoop,parentLoop,err,error,*999)
          IF(ASSOCIATED(parentLoop)) THEN
            CALL ControlLoop_IterationNumberGet(parentLoop,parentIteration,err,error,*999)
          ELSE
            parentIteration=0
          ENDIF
          NULLIFY(solvers)
          !Get the solver. For Biodomain problems of any split the 2nd solver will contain the dependent field equation set
          CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)            
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
          NULLIFY(solverEquations)
          CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
          NULLIFY(solverMapping)
          CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
          CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
          DO equationsSetIdx=1,numberOfEquationsSets
            NULLIFY(equationsSet)
            CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
            NULLIFY(region)
            CALL EquationsSet_RegionGet(equationsSet,region,err,error,*999)
            CALL Region_UserNumberGet(region,regionUserNumber,err,error,*999)        
            filename="Time_"//TRIM(NumberToVString(regionUserNumber,"*",err,error))//"_"// &
              & TRIM(NumberToVString(parentIteration,"*",err,error))//"_"// &
              & TRIM(NumberToVString(currentIteration,"*",err,error))
            method="FORTRAN"
            NULLIFY(fields)
            CALL Region_FieldsGet(region,fields,err,error,*999)
            CALL FIELD_IO_NODES_EXPORT(fields,filename,method,err,error,*999)
          ENDDO !equationsSetIdx
        ENDIF
      ENDIF
    CASE(CONTROL_WHILE_LOOP_TYPE)
      !do nothing
    CASE(CONTROL_LOAD_INCREMENT_LOOP_TYPE)
      !do nothing
    CASE DEFAULT
      localError="The control loop type of "//TRIM(NumberToVString(loopType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
   
    EXITS("Biodomain_PostLoop")
    RETURN
999 ERRORSEXITS("Biodomain_PostLoop",err,error)
    RETURN 1
    
  END SUBROUTINE Biodomain_PostLoop

  !
  !================================================================================================================================
  !

  !>Sets up the bioelectric domain equation type of a bioelectric equations set class.
  SUBROUTINE Biodomain_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to setup a bioelectric domain equation on.
    TYPE(EquationsSetSetupType), INTENT(INOUT) :: equationsSetSetup !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,dimensionIdx,dimensionMultiplier,esSpecification(3),esSpecificationType,esSpecificationSubtype, &
      & geometricComponentNumber,geometricMeshComponent,geometricScalingType,lumpingType,numberOfDimensions, &
      & numberOfIndependentComponents,numberOfMaterialsComponents,solutionMethod,sparsityType
    TYPE(DecompositionType), POINTER :: geometricDecomposition
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsSetIndependentType), POINTER :: equationsIndependent
    TYPE(EquationsSetMaterialsType), POINTER :: equationsMaterials
    TYPE(FieldType), POINTER :: geometricField
    TYPE(RegionType), POINTER :: region
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Biodomain_EquationsSetSetup",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)

    esSpecificationType=esSpecification(2)
    esSpecificationSubtype=esSpecification(3)
    
    SELECT CASE(esSpecificationType)
    CASE(EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE)
      SELECT CASE(esSpecificationSubtype)
      CASE(EQUATIONS_SET_MONODOMAIN_CELLML_SUBTYPE)
        !OK
      CASE(EQUATIONS_SET_MONODOMAIN_BUENOOROVIO_SUBTYPE, &
        & EQUATIONS_SET_MONODOMAIN_TENTUSSCHER06_SUBTYPE)
        !OK
      CASE(EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE, &
        & EQUATIONS_SET_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE, &
        & EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
        & EQUATIONS_SET_1D3D_MONODOMAIN_ACTIVE_STRAIN_SUBTYPE)
        !OK
      CASE DEFAULT
        localError="The third equations set specification of "//TRIM(NumberToVString(esSpecificationSubtype,"*",err,error))// &
          & " is not valid for a monodomain equation set."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE)
      SELECT CASE(esSpecificationSubtype)
      CASE(EQUATIONS_SET_FIRST_BIDOMAIN_SUBTYPE, &
        & EQUATIONS_SET_SECOND_BIDOMAIN_SUBTYPE)
        !OK
      CASE DEFAULT
        localError="The equations set subtype of "//TRIM(NumberToVString(esSpecificationSubtype,"*",err,error))// &
          & " is invalid for a bidomain equations set type."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The equation set type of "//TRIM(NumberToVString(esSpecificationType,"*",err,error))// &
        & " is invalid for a biodomain equations set class."
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
        CALL Biodomain_EquationsSetSolutionMethodSet(equationsSet,EQUATIONS_SET_FEM_SOLUTION_METHOD,err,error,*999)
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
!!Todo: CHECK VALID SETUP
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a bioelectric domain equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
      !-----------------------------------------------------------------
      ! G e o m e t r i c   f i e l d
      !-----------------------------------------------------------------
      !\todo Check geometric dimension
    CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
      !-----------------------------------------------------------------
      ! D e p e n d e n t   f i e l d
      !-----------------------------------------------------------------
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        SELECT CASE(esSpecificationType)
        CASE(EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE)
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
            SELECT CASE(esSpecificationSubtype)
            CASE(EQUATIONS_SET_MONODOMAIN_CELLML_SUBTYPE, &
              & EQUATIONS_SET_MONODOMAIN_BUENOOROVIO_SUBTYPE, &
              & eQUATIONS_SET_MONODOMAIN_TENTUSSCHER06_SUBTYPE)
              CALL Field_NumberOfVariablesSetAndLock(equationsSet%dependent%dependentField,2,err,error,*999)
              CALL Field_VariableTypesSetAndLock(equationsSet%dependent%dependentField,[FIELD_U_VARIABLE_TYPE, &
                & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
            CASE(EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE,EQUATIONS_SET_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE, &
              & EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE,EQUATIONS_SET_1D3D_MONODOMAIN_ACTIVE_STRAIN_SUBTYPE)
              CALL Field_NumberOfVariablesSetAndLock(equationsSet%dependent%dependentField,3,err,error,*999)
              CALL Field_VariableTypesSetAndLock(equationsSet%dependent%dependentField,[FIELD_U_VARIABLE_TYPE, &
                & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE],err,error,*999)
              CALL Field_VariableLabelSet(equationsSet%dependent%dependentField,FIELD_V_VARIABLE_TYPE,"GeometryM3D",ERR, &
                & ERROR,*999)
              CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_V_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_V_VARIABLE_TYPE,3, &
                & err,error,*999)
            CASE DEFAULT
              localError="The third equations set specification of "// &
                & TRIM(NumberToVString(esSpecificationSubtype,"*",err,error))// &
                & " is not valid for a monodomain equation set."
              CALL FlagError(localError,err,error,*999)
            END SELECT
            CALL Field_VariableLabelSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,"Vm",err,error,*999)
            CALL Field_VariableLabelSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,"dVm/dn", &
              & err,error,*999)
            CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
              & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_DP_TYPE,err,error,*999)
            CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
              & FIELD_DP_TYPE,err,error,*999)
            CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,1, &
              & err,error,*999)
            CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,1, &
              & err,error,*999)
            !Default to the geometric interpolation setup
            CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
            CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,1, &
              & geometricMeshComponent,err,error,*999)
            CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,1, &
              & geometricMeshComponent,err,error,*999)
            CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
            SELECT CASE(solutionMethod)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,1, &
                & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
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
            SELECT CASE(esSpecificationSubtype)
            CASE(EQUATIONS_SET_NO_SUBTYPE)
              CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,2,err,error,*999)
              CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE], &
                & err,error,*999)
            CASE(EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE,EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
              & EQUATIONS_SET_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE,EQUATIONS_SET_1D3D_MONODOMAIN_ACTIVE_STRAIN_SUBTYPE)
              CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,3,err,error,*999)
              CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_V_VARIABLE_TYPE],err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,3,err,error,*999)
            CASE DEFAULT
              localError="The third equations set specification of "// &
                & TRIM(NumberToVString(esSpecificationSubtype,"*",err,error))// &
                & " is not valid for a monodomain equation set."
              CALL FlagError(localError,err,error,*999)
            END SELECT
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
              & err,error,*999)
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
              & err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,1,err,error,*999)
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,1,err,error,*999)
            CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
            SELECT CASE(solutionMethod)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,1, &
                & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,1, &
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
              localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ENDIF
        CASE(EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE)
          SELECT CASE(esSpecificationSubtype)
          CASE(EQUATIONS_SET_FIRST_BIDOMAIN_SUBTYPE)
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
              CALL Field_NumberOfVariablesSetAndLock(equationsSet%dependent%dependentField,4,err,error,*999)
!! \todo allow for no rhs variable and so eliminate one of the flux variables
              CALL Field_VariableTypesSetAndLock(equationsSet%dependent%dependentField,[FIELD_U_VARIABLE_TYPE, &
                & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE,FIELD_DELVDELN_VARIABLE_TYPE],err,error,*999)
              CALL Field_VariableLabelSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,"Vm",err,error,*999)
              CALL Field_VariableLabelSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,"dVm/dn", &
                & err,error,*999)
              CALL Field_VariableLabelSet(equationsSet%dependent%dependentField,FIELD_V_VARIABLE_TYPE,"Phie",err,error,*999)
              CALL Field_VariableLabelSet(equationsSet%dependent%dependentField,FIELD_DELVDELN_VARIABLE_TYPE,"dPhie/dn", &
                & err,error,*999)
              CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_V_VARIABLE_TYPE, &
                & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_DELVDELN_VARIABLE_TYPE, &
                & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_V_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_DELVDELN_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,1, &
                & err,error,*999)
              CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                & err,error,*999)
              CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_V_VARIABLE_TYPE,1, &
                & err,error,*999)
              CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_DELVDELN_VARIABLE_TYPE,1, &
                & err,error,*999)
              !Default to the geometric interpolation setup
              CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                & geometricMeshComponent,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,1, &
                & geometricMeshComponent,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_V_VARIABLE_TYPE,1, &
                & geometricMeshComponent,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_DELVDELN_VARIABLE_TYPE,1, &
                & geometricMeshComponent,err,error,*999)
              CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
             SELECT CASE(solutionMethod)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField,FIELD_V_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField,FIELD_DELVDELN_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
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
!! \todo allow for no rhs variable and so eliminate one of the flux variables
              CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,4,err,error,*999)
              CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_V_VARIABLE_TYPE,FIELD_DELVDELN_VARIABLE_TYPE],err,error,*999)
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                & err,error,*999)
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                & err,error,*999)
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                & err,error,*999)
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
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_DELVDELN_VARIABLE_TYPE,1, &
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
                localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDIF
          CASE(EQUATIONS_SET_SECOND_BIDOMAIN_SUBTYPE)
            IF(equationsSet%dependent%dependentFieldAutoCreated) THEN
              CALL FlagError("The dependent field for the second bidomain equation cannot be auto-created. "// &
                & "You must pass in the field from the first bidomain equation.",err,error,*999)
            ELSE
              CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
              CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_DEPENDENT_TYPE,err,error,*999)
!! \todo allow for no rhs variable and so eliminate one of the flux variables
              CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,4,err,error,*999)
              CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_V_VARIABLE_TYPE,FIELD_DELVDELN_VARIABLE_TYPE],err,error,*999)
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                & err,error,*999)
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                & err,error,*999)
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                & err,error,*999)
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
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_DELVDELN_VARIABLE_TYPE,1, &
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
                localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDIF
          CASE DEFAULT
            localError="The equations set subtype of "//TRIM(NumberToVString(esSpecificationSubtype,"*",err,error))// &
              & " is invalid for a bidomain equations set type."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The equation set type of "//TRIM(NumberToVString(esSpecificationType,"*",err,error))// &
            & " is invalid for a biodomain equations set class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        IF(equationsSet%dependent%dependentFieldAutoCreated) THEN
          CALL Field_CreateFinish(equationsSet%dependent%dependentField,err,error,*999)
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a bioelectric domain equation"
        CALL FlagError(localError,err,error,*999)
      END SELECT      
    CASE(EQUATIONS_SET_SETUP_INDEPENDENT_TYPE)
      !-----------------------------------------------------------------
      ! I n d e p e n d e n t   f i e l d
      !-----------------------------------------------------------------
      NULLIFY(equationsIndependent)
      CALL EquationsSet_IndependentGet(equationsSet,equationsIndependent,err,error,*999)
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        SELECT CASE(esSpecificationType)
        CASE(EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE)
          SELECT CASE(esSpecificationSubtype)
          CASE(EQUATIONS_SET_MONODOMAIN_CELLML_SUBTYPE)
            !Do nothing
          CASE(EQUATIONS_SET_MONODOMAIN_BUENOOROVIO_SUBTYPE, &
            & EQUATIONS_SET_MONODOMAIN_TENTUSSCHER06_SUBTYPE)
            IF(esSpecificationSubtype==EQUATIONS_SET_MONODOMAIN_BUENOOROVIO_SUBTYPE) THEN
              numberOfIndependentComponents=4
            ELSE
              numberOfIndependentComponents=19
            ENDIF
            IF(equationsIndependent%independentFieldAutoCreated) THEN
              CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsIndependent%independentField,err,error,*999)
              CALL Field_TypeSetAndLock(equationsIndependent%independentField,FIELD_GENERAL_TYPE,err,error,*999)
              CALL Field_DependentTypeSetAndLock(equationsIndependent%independentField,FIELD_INDEPENDENT_TYPE,err,error,*999)
              NULLIFY(geometricDecomposition)
              CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
              CALL Field_DecompositionSetAndLock(equationsIndependent%independentField,geometricDecomposition,err,error,*999)
              CALL Field_GeometricFieldSetAndLock(equationsIndependent%independentField,geometricField,err,error,*999)
              CALL Field_NumberOfVariablesSetAndLock(equationsIndependent%independentField,1,err,error,*999)
              CALL Field_VariableTypesSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE,err,error,*999)
              
              CALL Field_DimensionSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE, &
                & err,error,*999)         
              CALL Field_NumberOfComponentsSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE,&
                & numberOfIndependentComponents,err,error,*999)
              CALL Field_DOFOrderTypeSet(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE,&
              &  FIELD_CONTIGUOUS_COMPONENT_DOF_ORDER,err,error,*999) ! dofs continuous, so first + (x-1) is x'th component index
              !Default to the geometric interpolation setup
              CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE,1, &
                & geometricMeshComponent,err,error,*999)
              CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
              SELECT CASE(solutionMethod)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField, &
                  & FIELD_U_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                !Default the scaling to the geometric field scaling
                CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
                CALL Field_ScalingTypeSet(equationsIndependent%independentField,geometricScalingType,err,error,*999)
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
              ! user specified field
              CALL FlagError("No user specified field supported!",err,error,*999)
            ENDIF
         CASE(EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE, &
            & EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
            & EQUATIONS_SET_1D3D_MONODOMAIN_ACTIVE_STRAIN_SUBTYPE)
            IF(equationsIndependent%independentFieldAutoCreated) THEN
              !Create the auto created independent field
              CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsIndependent%independentField,err,error,*999)
              CALL Field_LabelSet(equationsIndependent%independentField,"Independent Field",err,error,*999)
              CALL Field_TypeSetAndLock(equationsIndependent%independentField,FIELD_GENERAL_TYPE,err,error,*999)
              CALL Field_DependentTypeSetAndLock(equationsIndependent%independentField,FIELD_INDEPENDENT_TYPE,err,error,*999)
              NULLIFY(geometricDecomposition)
              CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
              CALL Field_DecompositionSetAndLock(equationsIndependent%independentField,geometricDecomposition,err,error,*999)
              CALL Field_GeometricFieldSetAndLock(equationsIndependent%independentField,geometricField,err,error,*999)
              CALL Field_NumberOfVariablesSetAndLock(equationsIndependent%independentField,4,err,error,*999)
              CALL Field_VariableTypesSetAndLock(equationsIndependent%independentField,[FIELD_U_VARIABLE_TYPE, &
                & FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE,FIELD_U2_VARIABLE_TYPE],err,error,*999)
              IF(esSpecificationSubtype==EQUATIONS_SET_1D3D_MONODOMAIN_ACTIVE_STRAIN_SUBTYPE) THEN
                CALL Field_VariableLabelSet(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE,"XB_state_variables", &
                  & err,error,*999)
              ELSE                
                CALL Field_VariableLabelSet(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE,"Active_stress", &
                  & err,error,*999)
              ENDIF
              CALL Field_VariableLabelSet(equationsIndependent%independentField,FIELD_U1_VARIABLE_TYPE,"sarcomere half length", &
                & err,error,*999)
              CALL Field_VariableLabelSet(equationsIndependent%independentField,FIELD_U2_VARIABLE_TYPE,"contraction velocity" &
                & ,err,error,*999)
              IF(esSpecificationSubtype==EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE) THEN
                CALL Field_DimensionSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
              ENDIF
              CALL Field_DataTypeSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE, &
                & err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsIndependent%independentField,FIELD_U1_VARIABLE_TYPE,FIELD_DP_TYPE, &
                & err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsIndependent%independentField,FIELD_U2_VARIABLE_TYPE,FIELD_DP_TYPE, &
                & err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsIndependent%independentField,FIELD_V_VARIABLE_TYPE,FIELD_INTG_TYPE, &
                & err,error,*999)
              IF(esSpecificationSubtype==EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE) THEN
                CALL Field_NumberOfComponentsSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                  & 1,err,error,*999)
              ELSEIF(esSpecificationSubtype==EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE) THEN
                CALL Field_NumberOfComponentsSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                  & 6,err,error,*999)
              ELSEIF(esSpecificationSubtype==EQUATIONS_SET_1D3D_MONODOMAIN_ACTIVE_STRAIN_SUBTYPE) THEN
                CALL Field_NumberOfComponentsSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                  & 4,err,error,*999)
              ENDIF
              CALL Field_NumberOfComponentsSetAndLock(equationsIndependent%independentField,FIELD_V_VARIABLE_TYPE,5, &
                & err,error,*999)
              IF(esSpecificationSubtype==EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE) THEN
                CALL Field_NumberOfComponentsSetAndLock(equationsIndependent%independentField,FIELD_U1_VARIABLE_TYPE, &
                  & 4,err,error,*999)
              ELSE
                CALL Field_NumberOfComponentsSetAndLock(equationsIndependent%independentField,FIELD_U1_VARIABLE_TYPE, &
                  & 3,err,error,*999)
              ENDIF
              CALL Field_NumberOfComponentsSetAndLock(equationsIndependent%independentField,FIELD_U2_VARIABLE_TYPE, &
                & 6,err,error,*999)
              !Default to the geometric interpolation setup
              CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE,1, &
                & geometricMeshComponent,err,error,*999)
              IF(esSpecificationSubtype==EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE) THEN
                CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE,2, &
                  & geometricMeshComponent,err,error,*999)
                CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE,3, &
                  & geometricMeshComponent,err,error,*999)
                CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE,4, &
                  & geometricMeshComponent,err,error,*999)
                CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE,5, &
                  & geometricMeshComponent,err,error,*999)
                CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE,6, &
                  & geometricMeshComponent,err,error,*999)
              ENDIF
              CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_V_VARIABLE_TYPE,1, &
                & geometricMeshComponent,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_V_VARIABLE_TYPE,2, &
                & geometricMeshComponent,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_V_VARIABLE_TYPE,3, &
                & geometricMeshComponent,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_V_VARIABLE_TYPE,4, &
                & geometricMeshComponent,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_V_VARIABLE_TYPE,5, &
                & geometricMeshComponent,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_U1_VARIABLE_TYPE,1, &
                & geometricMeshComponent,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_U1_VARIABLE_TYPE,2, &
                & geometricMeshComponent,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_U1_VARIABLE_TYPE,3, &
                & geometricMeshComponent,err,error,*999)
              IF(esSpecificationSubtype==EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE) THEN
                CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_U1_VARIABLE_TYPE, &
                  & 4,geometricMeshComponent,err,error,*999)
              ENDIF
              CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_U2_VARIABLE_TYPE,1, &
                & geometricMeshComponent,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_U2_VARIABLE_TYPE,2, &
                & geometricMeshComponent,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_U2_VARIABLE_TYPE,3, &
                & geometricMeshComponent,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_U2_VARIABLE_TYPE,4, &
                & geometricMeshComponent,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_U2_VARIABLE_TYPE,5, &
                & geometricMeshComponent,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_U2_VARIABLE_TYPE,6, &
                & geometricMeshComponent,err,error,*999)
              CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
              SELECT CASE(solutionMethod)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField, &
                  & FIELD_U_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                IF(esSpecificationSubtype==EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE) THEN
                  CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField, &
                    & FIELD_U_VARIABLE_TYPE,2,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField, &
                    & FIELD_U_VARIABLE_TYPE,3,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField, &
                    & FIELD_U_VARIABLE_TYPE,4,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField, &
                    & FIELD_U_VARIABLE_TYPE,5,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField, &
                    & FIELD_U_VARIABLE_TYPE,6,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                ENDIF
                CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField, &
                  & FIELD_V_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField, &
                  & FIELD_V_VARIABLE_TYPE,2,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField, &
                  & FIELD_V_VARIABLE_TYPE,3,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField, &
                  & FIELD_V_VARIABLE_TYPE,4,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField, &
                  & FIELD_V_VARIABLE_TYPE,5,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField, &
                  & FIELD_U1_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField, &
                  & FIELD_U1_VARIABLE_TYPE,2,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField, &
                  & FIELD_U1_VARIABLE_TYPE,3,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                IF(esSpecificationSubtype==EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE) THEN
                  CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField, &
                    & FIELD_U1_VARIABLE_TYPE,4,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                ENDIF
                CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField, &
                  & FIELD_U2_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField, &
                  & FIELD_U2_VARIABLE_TYPE,2,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField, &
                  & FIELD_U2_VARIABLE_TYPE,3,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField, &
                  & FIELD_U2_VARIABLE_TYPE,4,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField, &
                  & FIELD_U2_VARIABLE_TYPE,5,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField, &
                  & FIELD_U2_VARIABLE_TYPE,6,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                !Default the scaling to the geometric field scaling
                CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
                CALL Field_ScalingTypeSet(equationsIndependent%independentField,geometricScalingType,err,error,*999)
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
              CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
              CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,4,err,error,*999)
              CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE, &
                & FIELD_U1_VARIABLE_TYPE,FIELD_U2_VARIABLE_TYPE],err,error,*999)
              IF(esSpecificationSubtype==EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE) THEN
                CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                  & err,error,*999)
              ENDIF
              CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,FIELD_INTG_TYPE,err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U1_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U2_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              IF(esSpecificationSubtype==EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE) THEN
                CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,1,err,error,*999)
              ELSE IF(esSpecificationSubtype==EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE) THEN
                CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,6,err,error,*999)
              ELSEIF(esSpecificationSubtype==EQUATIONS_SET_1D3D_MONODOMAIN_ACTIVE_STRAIN_SUBTYPE) THEN
                CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,4,err,error,*999)
              ENDIF
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,5,err,error,*999)
              IF(esSpecificationSubtype==EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE) THEN
                CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U1_VARIABLE_TYPE,4,err,error,*999)
              ELSE
                CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U1_VARIABLE_TYPE,3,err,error,*999)
              ENDIF
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U2_VARIABLE_TYPE,6,err,error,*999)
              CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
              SELECT CASE(solutionMethod)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                IF(esSpecificationSubtype==EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE) THEN
                  CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,2, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,3, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,4, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,5, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,6, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                ENDIF
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,2, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,3, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,4, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,5, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U1_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U1_VARIABLE_TYPE,2, &
                  & FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U1_VARIABLE_TYPE,3, &
                  & FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                IF(esSpecificationSubtype==EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE) THEN
                  CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U1_VARIABLE_TYPE,4, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                ENDIF
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U2_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U2_VARIABLE_TYPE,2, &
                  & FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U2_VARIABLE_TYPE,3, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U2_VARIABLE_TYPE,4, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U2_VARIABLE_TYPE,5, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U2_VARIABLE_TYPE,6, &
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
                localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDIF
          CASE(EQUATIONS_SET_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE)
            IF(equationsIndependent%independentFieldAutoCreated) THEN
              !Create the auto created independent field
              CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsSet%INDEPENDENT%independentField, &
                & err,error,*999)
              CALL Field_LabelSet(equationsIndependent%independentField,"Independent Field",err,error,*999)
              CALL Field_TypeSetAndLock(equationsIndependent%independentField,FIELD_GENERAL_TYPE,err,error,*999)
              CALL Field_DependentTypeSetAndLock(equationsIndependent%independentField,FIELD_INDEPENDENT_TYPE,err,error,*999)
              NULLIFY(geometricDecomposition)
              CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
              CALL Field_DecompositionSetAndLock(equationsIndependent%independentField,geometricDecomposition, err,error,*999)
              CALL Field_GeometricFieldSetAndLock(equationsIndependent%independentField,geometricField,err,error,*999)
              CALL Field_NumberOfVariablesSetAndLock(equationsIndependent%independentField,2,err,error,*999)
              CALL Field_VariableTypesSetAndLock(equationsIndependent%independentField,[FIELD_U_VARIABLE_TYPE, &
                & FIELD_V_VARIABLE_TYPE],err,error,*999)
              CALL Field_VariableLabelSet(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE,"Active_stress", &
                & err,error,*999)
              CALL Field_VariableLabelSet(equationsIndependent%independentField,FIELD_V_VARIABLE_TYPE,"fibre_info", &
                & err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE, &
                & err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsIndependent%independentField,FIELD_V_VARIABLE_TYPE,FIELD_INTG_TYPE, &
                & err,error,*999)
              CALL Field_NumberOfComponentsSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                & 1,err,error,*999)
              CALL Field_NumberOfComponentsSetAndLock(equationsIndependent%independentField,FIELD_V_VARIABLE_TYPE,5, &
                & err,error,*999)
              !Default to the geometric interpolation setup
              CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE,1, &
                & geometricMeshComponent,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_V_VARIABLE_TYPE,1, &
                & geometricMeshComponent,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_V_VARIABLE_TYPE,2, &
                & geometricMeshComponent,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_V_VARIABLE_TYPE,3, &
                & geometricMeshComponent,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_V_VARIABLE_TYPE,4, &
                & geometricMeshComponent,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_V_VARIABLE_TYPE,5, &
                & geometricMeshComponent,err,error,*999)
              CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
              SELECT CASE(solutionMethod)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField, &
                  & FIELD_U_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField, &
                  & FIELD_V_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField, &
                  & FIELD_V_VARIABLE_TYPE,2,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField, &
                  & FIELD_V_VARIABLE_TYPE,3,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField, &
                  & FIELD_V_VARIABLE_TYPE,4,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField, &
                  & FIELD_V_VARIABLE_TYPE,5,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                !Default the scaling to the geometric field scaling
                CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
                CALL Field_ScalingTypeSet(equationsIndependent%independentField,geometricScalingType,err,error,*999)
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
              CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
              CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,2,err,error,*999)
              CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE],err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,FIELD_INTG_TYPE,err,error,*999)
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,1,err,error,*999)
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,5,err,error,*999)
              CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
              SELECT CASE(solutionMethod)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,2, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,3, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,4, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,5, &
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
                localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDIF
          CASE DEFAULT
            localError="The equations set subtype of "//TRIM(NumberToVString(esSpecificationSubtype,"*",err,error))// &
              " is not implemented for an equations set setup independent type."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE)
          localError="Equations set setup independent type is not implemented for an equations set bidomain equation type"
          CALL FlagError(localError,err,error,*999)
        CASE DEFAULT
          localError="The equation set type of "//TRIM(NumberToVString(esSpecificationType,"*",err,error))// &
            & " is invalid for a biodomain equations set class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        IF(equationsIndependent%independentFieldAutoCreated) THEN
          CALL Field_CreateFinish(equationsIndependent%independentField,err,error,*999)
          SELECT CASE(esSpecificationType)
          CASE(EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE)
            SELECT CASE(esSpecificationSubtype)
            CASE(EQUATIONS_SET_MONODOMAIN_CELLML_SUBTYPE)
              !Do nothing
            CASE(EQUATIONS_SET_MONODOMAIN_BUENOOROVIO_SUBTYPE)
              CALL Electrophysiology_BuenoOrovioInitialise(equationsIndependent%independentField,err,error,*999)
            CASE(EQUATIONS_SET_MONODOMAIN_TENTUSSCHER06_SUBTYPE)
              CALL Electrophysiology_TenTusscher06Initialise(equationsIndependent%independentField,err,error,*999)
            CASE(EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE, &
              & EQUATIONS_SET_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE, &
              & EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
              & EQUATIONS_SET_1D3D_MONODOMAIN_ACTIVE_STRAIN_SUBTYPE)
              !Do nothing
            CASE DEFAULT
              localError="The third equations set specification of "// &
                & TRIM(NumberToVString(esSpecificationSubtype,"*",err,error))// &
                & " is not valid for a monodomain equation set."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE(EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE)
            SELECT CASE(esSpecificationSubtype)
            CASE(EQUATIONS_SET_FIRST_BIDOMAIN_SUBTYPE, &
              & EQUATIONS_SET_SECOND_BIDOMAIN_SUBTYPE)
              !Do nothing
            CASE DEFAULT
              localError="The equations set subtype of "//TRIM(NumberToVString(esSpecificationSubtype,"*",err,error))// &
                & " is invalid for a bidomain equations set type."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The equation set type of "//TRIM(NumberToVString(esSpecificationType,"*",err,error))// &
              & " is invalid for a biodomain equations set class."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ENDIF        
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a bioelectric domain equation"
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
      SELECT CASE(esSpecificationType)
      CASE(EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE)
        SELECT CASE(esSpecificationSubtype)
        CASE(EQUATIONS_SET_MONODOMAIN_CELLML_SUBTYPE, &
          & EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE, &
          & EQUATIONS_SET_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE, &
          & EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
          & EQUATIONS_SET_1D3D_MONODOMAIN_ACTIVE_STRAIN_SUBTYPE)
          !Monodomain. Materials field components are 2 plus the conductivity tensor i.e., Am, Cm and \sigma
          numberOfMaterialsComponents=2+NUMBER_OF_VOIGT(numberOfDimensions)
          dimensionMultiplier=1
        CASE(EQUATIONS_SET_MONODOMAIN_BUENOOROVIO_SUBTYPE, &
          & EQUATIONS_SET_MONODOMAIN_TENTUSSCHER06_SUBTYPE)
          !Materials field components are
          ! 1. activation  factor (usually 0.0 or 1.0)
          ! 2,3 for fiber/transverse conductivity   . defaults to constant interpolation 
          ! 4,5[,6] : fiber unit vector in dimension
          ! 7: out - activation times
          numberOfMaterialsComponents=4+numberOfDimensions
          dimensionMultiplier=1
        CASE DEFAULT
          localError="The third equations set specification of "//TRIM(NumberToVString(esSpecificationSubtype,"*",err,error))// &
            & " is not valid for a monodomain equation set."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE)
        !Bidomain. Materials field components are 2 plus two conductivity tensors i.e., Am, C, \sigma_i and \sigma_e
!!TODO or make the extracellular field in the V variable?          
        numberOfMaterialsComponents=2+2*NUMBER_OF_VOIGT(numberOfDimensions)
        dimensionMultiplier=2
      CASE DEFAULT
        localError="The equations set type of "//TRIM(NumberToVString(esSpecificationType,"*",err,error))// &
          & " is invalid for a bioelectrics class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        CALL EquationsSet_AssertDependentIsFinished(equationsSet,err,error,*999)
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
          CALL Field_NumberOfVariablesSetAndLock(equationsMaterials%materialsField,1,err,error,*999)
          CALL Field_VariableTypesSetAndLock(equationsMaterials%materialsField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
          CALL Field_VariableLabelSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,"Materials",err,error,*999)
          CALL Field_DimensionSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
            & err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
          !Set the number of materials components
          CALL Field_NumberOfComponentsSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
            & numberOfMaterialsComponents,err,error,*999)
          CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricComponentNumber,err,error,*999)
          IF((esSpecification(2)==EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE).AND. &
            & (esSpecification(3)==EQUATIONS_SET_MONODOMAIN_BUENOOROVIO_SUBTYPE.OR. &
            & esSpecification(3)==EQUATIONS_SET_MONODOMAIN_TENTUSSCHER06_SUBTYPE)) THEN
            ! 1st = activation = node based
            CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,1, &
              & geometricComponentNumber,err,error,*999)
            CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,1, &
              & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            CALL Field_ComponentLabelSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,1,"Activation", &
              & err,error,*999)
            ! 2 3 diffusion constants
            CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,2, &
              & geometricComponentNumber,err,error,*999)
            CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,2, &
              & FIELD_CONSTANT_INTERPOLATION,err,error,*999)
            CALL Field_ComponentLabelSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,2,"Df",err,error,*999)
            CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,3, &
              & geometricComponentNumber,err,error,*999)
            CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,3, &
              & FIELD_CONSTANT_INTERPOLATION,err,error,*999)
            CALL Field_ComponentLabelSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,2,"Dt",err,error,*999)
            ! 4 (5 (6)) fiber unit vector
            CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,4, &
              & geometricComponentNumber,err,error,*999)
            CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,4, &
              & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            CALL Field_ComponentLabelSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,4,"fx",err,error,*999)
            IF(numberOfDimensions>1) THEN
              CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,5, &
                & geometricComponentNumber,err,error,*999)
              CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,5, &
                & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentLabelSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,5,"fy",err,error,*999)
              IF(numberOfDimensions>2) THEN
                CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,6, &
                  & geometricComponentNumber,err,error,*999)
                CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,6, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentLabelSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,6,"fZ",err,error,*999)
              ENDIF
            ENDIF
            CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
              & numberOfMaterialsComponents,geometricComponentNumber,err,error,*999)
            CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
              & numberOfMaterialsComponents,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            CALL Field_ComponentLabelSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
              & numberOfMaterialsComponents,"Activation times",err,error,*999)
          ELSE
            !Default the Am and Cm materials components to the first component geometric interpolation with const interpolation
            DO componentIdx=1,2
              CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,componentIdx, &
                & geometricComponentNumber,err,error,*999)
              CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,componentIdx, &
                & FIELD_CONSTANT_INTERPOLATION,err,error,*999)
            ENDDO !components_idx
            CALL Field_ComponentLabelSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,1,"Am",err,error,*999)
            CALL Field_ComponentLabelSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,2,"Cm",err,error,*999)
            !Default the \sigma materials components to the first geometric interpolation setup with constant interpolation
            CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricComponentNumber,err,error,*999)
            DO componentIdx=1,NUMBER_OF_VOIGT(numberOfDimensions)
              DO dimensionIdx=1,dimensionMultiplier
                CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                  & 2+componentIdx+(dimensionIdx-1)*numberOfDimensions,geometricComponentNumber,err,error,*999)
                CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                  & 2+componentIdx+(dimensionIdx-1)*numberOfDimensions,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
              ENDDO !dimensionIdx
            ENDDO !componentIdx
          ENDIF
          !Default the field scaling to that of the geometric field
          CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
          CALL Field_ScalingTypeSet(equationsMaterials%materialsField,geometricScalingType,err,error,*999)
        ELSE
          !Check the user specified field
          CALL Field_TypeCheck(equationsSetSetup%field,FIELD_MATERIAL_TYPE,err,error,*999)
          CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
          CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,1,err,error,*999)
          CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE],err,error,*999)
          CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfMaterialsComponents, &
            & err,error,*999)
        ENDIF
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        IF(equationsMaterials%materialsFieldAutoCreated) THEN
          !Finish creating the materials field
          CALL Field_CreateFinish(equationsMaterials%materialsField,err,error,*999)
          !Set the default values for the materials field
          IF((esSpecification(2)==EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE).AND. &
            & (esSpecification(3)==EQUATIONS_SET_MONODOMAIN_BUENOOROVIO_SUBTYPE.OR. &
            & esSpecification(3)==EQUATIONS_SET_MONODOMAIN_TENTUSSCHER06_SUBTYPE)) THEN
            DO componentIdx=1,numberOfMaterialsComponents
              CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,componentIdx,1.0_DP,err,error,*999)
            ENDDO !componentIdx
         ELSE
            !First set Am
            CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VALUES_SET_TYPE,1,200.0_DP,err,error,*999)
            !Now set Cm
            CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VALUES_SET_TYPE,2,0.0025_DP,err,error,*999)
            !Now set the sigmas to be 1.0 diagonal
            DO componentIdx=1,numberOfDimensions
              DO dimensionIdx=1,dimensionMultiplier
                CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,2+componentIdx+(dimensionIdx-1)*NUMBER_OF_VOIGT(numberOfDimensions),1.0_DP, &
                  & err,error,*999)
              ENDDO !dimensionIdx
            ENDDO !componentIdx
          ENDIF
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a bioelectric domain equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
      !-----------------------------------------------------------------
      ! S o u r c e   f i e l d
      !-----------------------------------------------------------------
!!TODO: allow for source current fields
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        !Do nothing
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        !Do nothing
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a bioelectric domain equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
      !-----------------------------------------------------------------
      ! A n a l y t i c  T y p e
      !-----------------------------------------------------------------
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a bioelectric domain equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
      !-----------------------------------------------------------------
      ! E q u a t i o n s    t y p e
      !-----------------------------------------------------------------     
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        CALL EquationsSet_AssertDependentIsFinished(equationsSet,err,error,*999)
        CALL EquationsSet_AssertMaterialsIsFInished(equationsSet,err,error,*999)
        !Create the equations
        NULLIFY(equations)
        CALL Equations_CreateStart(equationsSet,equations,err,error,*999)
        CALL Equations_LinearityTypeSet(equations,EQUATIONS_LINEAR,err,error,*999)
        SELECT CASE(esSpecificationType)
        CASE(EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE)
          CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
        CASE(EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE)
          SELECT CASE(esSpecificationSubtype)
          CASE(EQUATIONS_SET_FIRST_BIDOMAIN_SUBTYPE)
            CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
          CASE(EQUATIONS_SET_SECOND_BIDOMAIN_SUBTYPE)
            CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_STATIC,err,error,*999)
          CASE DEFAULT
            localError="The third equations set specification of "// &
              & TRIM(NumberToVString(esSpecificationSubtype,"*",err,error))// &
              & " is invalid for a bidomain equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The second equations set specification of "// &
            & TRIM(NumberToVString(esSpecificationType,"*",err,error))// &
            & " is invalid for a bioelectrics class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
        SELECT CASE(solutionMethod)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          !Finish the creation of the equations
          NULLIFY(equations)
          CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
          CALL Equations_CreateFinish(equations,err,error,*999)
          NULLIFY(vectorEquations)
          CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
          !Create the equations mapping.
          NULLIFY(vectorMapping)
          CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_U_VARIABLE_TYPE,vectorMapping,err,error,*999)
          SELECT CASE(esSpecificationType)
          CASE(EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE)
            CALL EquationsMappingVector_DynamicMatricesSet(vectorMapping,.TRUE.,.TRUE.,err,error,*999)
            CALL EquationsMappingVector_DynamicVariableTypeSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
            CALL EquationsMappingVector_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
          CASE(EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE)
            SELECT CASE(esSpecificationSubtype)
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
              localError="The equations set subtype of "//TRIM(NumberToVString(esSpecificationSubtype,"*",err,error))// &
                & " is invalid for a bidomain equation type."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The equations set type of "//TRIM(NumberToVString(esSpecificationType,"*",err,error))// &
              & " is invalid for a bioelectrics class."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
          !Create the equations matrices
          NULLIFY(vectorMatrices)
          CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
          CALL Equations_LumpingTypeGet(equations,lumpingType,err,error,*999)
          CALL Equations_SparsityTypeGet(equations,sparsityType,err,error,*999)
          SELECT CASE(esSpecificationType)
          CASE(EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE)
            !Set up matrix storage and structure
            IF(lumpingType==EQUATIONS_LUMPED_MATRICES) THEN
              !Set up lumping
              CALL EquationsMatricesVector_DynamicLumpingTypeSet(vectorMatrices, &
                & [EQUATIONS_MATRIX_UNLUMPED,EQUATIONS_MATRIX_LUMPED],err,error,*999)
              CALL EquationsMatricesVector_DynamicStorageTypeSet(vectorMatrices, &
                & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE], &
                & err,error,*999)
              CALL EquationsMatricesVector_DynamicStructureTypeSet(vectorMatrices, &
                & [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_DIAGONAL_STRUCTURE],err,error,*999)
            ELSE              
              SELECT CASE(sparsityType)
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
                localError="The equations matrices sparsity type of "// &
                  & TRIM(NumberToVString(sparsityType,"*",err,error))//" is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDIF
          CASE(EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE)
            SELECT CASE(esSpecificationSubtype)
            CASE(EQUATIONS_SET_FIRST_BIDOMAIN_SUBTYPE)
              !Set up matrix storage and structure
              IF(lumpingType==EQUATIONS_LUMPED_MATRICES) THEN
                !Set up lumping
                CALL EquationsMatricesVector_DynamicLumpingTypeSet(vectorMatrices, &
                  & [EQUATIONS_MATRIX_UNLUMPED,EQUATIONS_MATRIX_LUMPED],err,error,*999)
                CALL EquationsMatricesVector_DynamicStorageTypeSet(vectorMatrices, &
                  & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE], &
                  & err,error,*999)
                CALL EquationsMatricesVector_DynamicStructureTypeSet(vectorMatrices, &
                  & [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_DIAGONAL_STRUCTURE],err,error,*999)
              ELSE
                SELECT CASE(sparsityType)
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
                  localError="The equations matrices sparsity type of "// &
                    & TRIM(NumberToVString(sparsityType,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              ENDIF
            CASE(EQUATIONS_SET_SECOND_BIDOMAIN_SUBTYPE)
              SELECT CASE(sparsityType)
              CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,[MATRIX_BLOCK_STORAGE_TYPE], &
                  & err,error,*999)
              CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,[MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
                  & err,error,*999)
                CALL EquationsMatricesVector_LinearStructureTypeSet(vectorMatrices,[EQUATIONS_MATRIX_FEM_STRUCTURE], &
                  & err,error,*999)
              CASE DEFAULT
                localError="The equations matrices sparsity type of "// &
                  & TRIM(NumberToVString(sparsityType,"*",err,error))//" is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE DEFAULT
              localError="The equations set subtype of "//TRIM(NumberToVString(esSpecificationSubtype,"*",err,error))// &
                & " is invalid for a bidomain equation type."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The equations set type of "//TRIM(NumberToVString(esSpecificationType,"*",err,error))// &
              & " is invalid for a bioelectrics class."
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
          localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a bioelectric domain equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
        & " is invalid for a bioelectric domain equation."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("Biodomain_EquationsSetSetup")
    RETURN
999 ERRORSEXITS("Biodomain_EquationsSetSetup",err,error)
    RETURN 1
    
  END SUBROUTINE Biodomain_EquationsSetSetup

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a bioelectric domain equation type of an bioelectrics equations set class.
  SUBROUTINE Biodomain_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: solutionMethod !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: esSpecification(3),esSpecificationType,esSpecificationSubtype
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Biodomain_EquationsSetSolutionMethodSet",err,error,*999)
    
    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
    
    esSpecificationType=esSpecification(2)
    esSpecificationSubtype=esSpecification(3)
    SELECT CASE(esSpecificationType)
    CASE(EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE)        
      SELECT CASE(esSpecificationSubtype)
      CASE(EQUATIONS_SET_MONODOMAIN_CELLML_SUBTYPE, &
        & EQUATIONS_SET_MONODOMAIN_BUENOOROVIO_SUBTYPE, &
        & EQUATIONS_SET_MONODOMAIN_TENTUSSCHER06_SUBTYPE, &
        & EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE, &
        & EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
        & EQUATIONS_SET_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE, &
        & EQUATIONS_SET_1D3D_MONODOMAIN_ACTIVE_STRAIN_SUBTYPE)
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
        localError="Equations set subtype of "//TRIM(NumberToVString(esSpecificationSubtype,"*",err,error))// &
          & " is not valid for a bioelectric monodomain equation type of an bioelectrics equations set class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE)        
      SELECT CASE(esSpecificationSubtype)
      CASE(EQUATIONS_SET_FIRST_BIDOMAIN_SUBTYPE, &
        & EQUATIONS_SET_SECOND_BIDOMAIN_SUBTYPE)
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
        localError="Equations set subtype of "//TRIM(NumberToVString(esSpecificationSubtype,"*",err,error))// &
          & " is not valid for a bioelectric bidomain equation type of an bioelectrics equations set class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="Equations set type of "//TRIM(NumberToVString(esSpecificationType,"*",err,error))// &
        & " is not valid for a bioelectrics equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
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
    INTEGER(INTG) :: equationsSetType,equationsSetSubtype
    TYPE(VARYING_STRING) :: localError

    ENTERS("Biodomain_EquationsSetSpecificationSet",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(ALLOCATED(equationsSet%specification)) CALL FlagError("Equations set specification is already allocated.",err,error,*999)    
    IF(SIZE(specification,1)<3) THEN
      localError="The size of the specified specification array of "// &
        & TRIM(NumberToVString(SIZE(specification,1),"*",err,error))//" is invalid. The size should be >= 3."
      CALL FlagError(localError,err,error,*999)
    END IF
    
    equationsSetType=specification(2)
    equationsSetSubtype=specification(3)
    SELECT CASE(equationsSetType)
    CASE(EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE)
      SELECT CASE(equationsSetSubtype)
      CASE(EQUATIONS_SET_MONODOMAIN_CELLML_SUBTYPE)
        !ok
      CASE(EQUATIONS_SET_MONODOMAIN_BUENOOROVIO_SUBTYPE, &
        & EQUATIONS_SET_MONODOMAIN_TENTUSSCHER06_SUBTYPE)
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
      CASE(EQUATIONS_SET_FIRST_BIDOMAIN_SUBTYPE, &
        & EQUATIONS_SET_SECOND_BIDOMAIN_SUBTYPE)
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
    ALLOCATE(equationsSet%specification(3),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
    equationsSet%specification(1:3)=[EQUATIONS_SET_BIOELECTRICS_CLASS,equationsSetType,equationsSetSubtype]

    EXITS("Biodomain_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("Biodomain_EquationsSetSpecificationSet",err,error)
    EXITS("Biodomain_EquationsSetSpecificationSet")
    RETURN 1
    
  END SUBROUTINE Biodomain_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !

  !>Performs post-solve actions for mono- and bi-domain problems.
  SUBROUTINE Biodomain_PostSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to perform the post-solve actions for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: esSpecification(3),pSpecification(3)
    REAL(DP) :: currentTime,timeIncrement    
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatrixType), POINTER :: dampingMatrix,stiffnessMatrix
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField,geometricField,independentField,materialsField
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError

    ENTERS("Biodomain_PostSolve",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)

    CALL ControlLoop_CurrentTimesGet(controlLoop,currentTime,timeIncrement,err,error,*999)
    
    SELECT CASE(pSpecification(2))
    CASE(PROBLEM_MONODOMAIN_EQUATION_TYPE)
      SELECT CASE(pSpecification(3))
      CASE(PROBLEM_MONODOMAIN_GUDUNOV_SPLIT_SUBTYPE)
        !Do nothing
      CASE(PROBLEM_MONODOMAIN_STRANG_SPLIT_SUBTYPE)
        !Do nothing
      CASE(PROBLEM_MONODOMAIN_DIRECT_MODEL_SUBTYPE)
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        NULLIFY(solverMapping)
        CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
        NULLIFY(equations)
        CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
        NULLIFY(vectorEquations)
        CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
        NULLIFY(vectorMatrices)
        CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
        NULLIFY(dynamicMatrices)
        CALL EquationsMatricesVector_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*999)
        NULLIFY(stiffnessMatrix)
        CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,1,stiffnessMatrix,err,error,*999)
        NULLIFY(dampingMatrix)
        CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,2,dampingMatrix,err,error,*999)
        NULLIFY(rhsVector)
        CALL EquationsMatricesVector_RHSVectorGet(vectorMatrices,rhsVector,err,error,*999)
        CALL EquationsMatrix_UpdateMatrixSet(stiffnessMatrix,.FALSE.,err,error,*999)
        CALL EquationsMatrix_UpdateMatrixSet(dampingMatrix,.FALSE.,err,error,*999)
        CALL EquationsMatricesRHS_UpdateVectorSet(rhsVector,.FALSE.,err,error,*999)
       
        NULLIFY(geometricField)
        CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
        NULLIFY(dependentField)
        CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
        NULLIFY(materialsField)
        CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
        NULLIFY(independentField)
        CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)

        CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
        SELECT CASE(esSpecification(2))
        CASE(EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE)
          SELECT CASE(esSpecification(3))
          CASE(EQUATIONS_SET_MONODOMAIN_CELLML_SUBTYPE, &
            & EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE, &
            & EQUATIONS_SET_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE, &
            & EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
            & EQUATIONS_SET_1D3D_MONODOMAIN_ACTIVE_STRAIN_SUBTYPE)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(EQUATIONS_SET_MONODOMAIN_BUENOOROVIO_SUBTYPE)
            CALL Electrophysiology_BuenoOrovioIntegrate(independentField,materialsField,currentTime-timeIncrement,currentTime, &
              & err,error,*999)  ! from t-dt to t
          CASE(EQUATIONS_SET_MONODOMAIN_TENTUSSCHER06_SUBTYPE)
            CALL Electrophysiology_TenTusscher06Integrate(independentField,materialsField,currentTime-timeIncrement,currentTime, &
              & err,error,*999) ! from t-dt to t
          CASE DEFAULT
            localError="The third equations set specification of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
              & " is not valid for a monodomain equation set."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE)
          CALL FlagError("Not implemented",err,error,*999)
        CASE DEFAULT
          localError="The equation set type of "//TRIM(NumberToVString(esSpecification(2),"*",err,error))// &
            & " is invalid for a bioelectrics equations set class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The problem subtype of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
          & " is invalid for a monodomain problem type."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_BIDOMAIN_EQUATION_TYPE)
      SELECT CASE(pSpecification(3))
      CASE(PROBLEM_BIDOMAIN_GUDUNOV_SPLIT_SUBTYPE)
        !Do nothing
      CASE(PROBLEM_BIDOMAIN_STRANG_SPLIT_SUBTYPE)
        !Do nothing
      CASE DEFAULT
        localError="The problem subtype of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
          & " is invalid for a bidomain problem type."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE)
      SELECT CASE(pSpecification(3))
      CASE(PROBLEM_GUDUNOV_MONODOMAIN_SIMPLE_ELASTICITY_SUBTYPE,PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE, &
        & PROBLEM_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE,EQUATIONS_SET_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE, &
        & PROBLEM_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE)
        !Do nothing
      CASE DEFAULT
        localError="The problem subtype of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
          & " is invalid for a monodomain problem type."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The problem type of "//TRIM(NumberToVString(pSpecification(2),"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
      
    EXITS("Biodomain_PostSolve")
    RETURN
999 ERRORSEXITS("Biodomain_PostSolve",err,error)
    RETURN 1
    
  END SUBROUTINE Biodomain_PostSolve

  !
  !================================================================================================================================
  !

  !>Performs pre-solve actions for mono- and bi-domain problems.
  SUBROUTINE Biodomain_PreSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to perform the pre-solve actions for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(3),solverNumber
    REAL(DP) :: currentTime,timeIncrement
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldType), POINTER :: dependentField,independentField
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError

    ENTERS("Biodomain_PreSolve",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)

    CALL Solver_GlobalNumberGet(solver,solverNumber,err,error,*999)
    CALL ControlLoop_CurrentTimesGet(controlLoop,currentTime,timeIncrement,err,error,*999)
    
    SELECT CASE(pSpecification(2))
    CASE(PROBLEM_MONODOMAIN_EQUATION_TYPE)
      SELECT CASE(pSpecification(3))
      CASE(PROBLEM_MONODOMAIN_GUDUNOV_SPLIT_SUBTYPE)
        SELECT CASE(solverNumber)
        CASE(1)
          CALL Solver_DAETimesSet(solver,currentTime,currentTime+timeIncrement,err,error,*999)
        CASE(2)
          !Do nothing
        CASE DEFAULT
          localError="The solver global number of "//TRIM(NumberToVString(solverNumber,"*",err,error))// &
            & " is invalid for a Gudunov split monodomain problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PROBLEM_MONODOMAIN_STRANG_SPLIT_SUBTYPE)
        SELECT CASE(solverNumber)
        CASE(1)
          CALL Solver_DAETimesSet(solver,currentTime,currentTime+timeIncrement/2.0_DP,err,error,*999)
        CASE(2)
          !Do nothing
        CASE(3)
          CALL Solver_DAETimesSet(solver,currentTime+timeIncrement/2.0_DP,currentTime+timeIncrement,err,error,*999)
        CASE DEFAULT
          localError="The solver global number of "//TRIM(NumberToVString(solverNumber,"*",err,error))// &
            & " is invalid for a Strang split monodomain problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PROBLEM_MONODOMAIN_DIRECT_MODEL_SUBTYPE)
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        NULLIFY(solverMapping)
        CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
        NULLIFY(dependentField)
        CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
        NULLIFY(independentField)
        CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)
        CALL Field_ParametersToFieldParametersCopy(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1, &
          & dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,err,error,*999)
        CALL Field_ParametersToFieldParametersCopy(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1, &
          & dependentField,FIELD_U_VARIABLE_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE,1,err,error,*999) ! also to prev.
      CASE DEFAULT
        localError="The problem subtype of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
          & " is invalid for a monodomain problem type."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_BIDOMAIN_EQUATION_TYPE)
      SELECT CASE(pSpecification(3))
      CASE(PROBLEM_BIDOMAIN_GUDUNOV_SPLIT_SUBTYPE)
        SELECT CASE(solverNumber)
        CASE(1)
          CALL Solver_DAETimesSet(solver,currentTime,currentTime+timeIncrement,err,error,*999)
        CASE(2)
          !Do nothing
        CASE(3)
          !Do nothing
        CASE DEFAULT
          localError="The solver global number of "//TRIM(NumberToVString(solverNumber,"*",err,error))// &
            & " is invalid for a Gudunov split bidomain problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PROBLEM_BIDOMAIN_STRANG_SPLIT_SUBTYPE)
        SELECT CASE(solverNumber)
        CASE(1)
          CALL Solver_DAETimesSet(solver,currentTime,currentTime+timeIncrement/2.0_DP,err,error,*999)
        CASE(2)
          !Do nothing
        CASE(3)
          !Do nothing
        CASE(4)
          CALL Solver_DAETimesSet(solver,currentTime+timeIncrement/2.0_DP,currentTime+timeIncrement,err,error,*999)
        CASE DEFAULT
          localError="The solver global number of "//TRIM(NumberToVString(solverNumber,"*",err,error))// &
            & " is invalid for a Gudunov split bidomain problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The problem subtype of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
          & " is invalid for a bidomain problem type."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE)
      SELECT CASE(pSpecification(3))
      CASE(PROBLEM_GUDUNOV_MONODOMAIN_SIMPLE_ELASTICITY_SUBTYPE,PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE, &
        & PROBLEM_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE,EQUATIONS_SET_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE, &
        & PROBLEM_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE)
        SELECT CASE(solverNumber)
        CASE(1)
          CALL Solver_DAETimesSet(solver,currentTime,currentTime+timeIncrement,err,error,*999)
        CASE(2)
          !Do nothing
        CASE DEFAULT
          localError="The solver global number of "//TRIM(NumberToVString(solverNumber,"*",err,error))// &
            & " is invalid for a bioelectrics finite elasticity problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The problem subtype of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
          & " is invalid for a monodomain problem type."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The problem type of "//TRIM(NumberToVString(pSpecification(2),"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
      
    EXITS("Biodomain_PreSolve")
    RETURN
999 ERRORSEXITS("Biodomain_PreSolve",err,error)
    RETURN 1
    
  END SUBROUTINE Biodomain_PreSolve

  !
  !================================================================================================================================
  !
 
  !>Sets up the bioelectric domain problem.
  SUBROUTINE Biodomain_ProblemSetup(problem,problemSetup,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem set to setup a bioelectric domain equation on.
    TYPE(ProblemSetupType), INTENT(INOUT) :: problemSetup !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(3)
    TYPE(CellMLEquationsType), POINTER :: cellMLEquations
    TYPE(ControlLoopType), POINTER :: controlLoop,controlLoopRoot
    TYPE(SolverType), POINTER :: solver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError
   
    ENTERS("Biodomain_ProblemSetup",err,error,*999)

    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    
    SELECT CASE(pSpecification(2))
    CASE(PROBLEM_MONODOMAIN_EQUATION_TYPE)
      SELECT CASE(pSpecification(3))
      CASE(PROBLEM_MONODOMAIN_GUDUNOV_SPLIT_SUBTYPE, &
        & PROBLEM_MONODOMAIN_STRANG_SPLIT_SUBTYPE, &
        & PROBLEM_MONODOMAIN_DIRECT_MODEL_SUBTYPE)
        !OK
      CASE DEFAULT
        localError="The problem subtype of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
          & " is invalid for a monodomain problem type of a bioelectric problem class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_BIDOMAIN_EQUATION_TYPE)
      SELECT CASE(pSpecification(3))


        
      CASE(PROBLEM_BIDOMAIN_GUDUNOV_SPLIT_SUBTYPE, &
        & PROBLEM_BIDOMAIN_STRANG_SPLIT_SUBTYPE)
        !OK
      CASE DEFAULT
        localError="The problem subtype of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))//  &
          & " is invalid for a bidomain problem type."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The problem type of "//TRIM(NumberToVString(pSpecification(2),"*",err,error))// &
        & " is invalid for a bioelectric problem class."
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
          & " is invalid for a bioelectric domain equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_CONTROL_TYPE)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Set up a time control loop
        NULLIFY(controlLoop)
        CALL ControlLoop_CreateStart(problem,controlLoop,err,error,*999)
        CALL ControlLoop_TypeSet(controlLoop,CONTROL_TIME_LOOP_TYPE,err,error,*999)
        CALL ControlLoop_LabelSet(controlLoop,"Time Loop",err,error,*999)
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
          & " is invalid for a bioelectric domain equation."
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
        SELECT CASE(pSpecification(2))
        CASE(PROBLEM_MONODOMAIN_EQUATION_TYPE)
          SELECT CASE(pSpecification(3))
          CASE(PROBLEM_MONODOMAIN_GUDUNOV_SPLIT_SUBTYPE)
            CALL Solvers_NumberOfSolversSet(solvers,2,err,error,*999)
            !Set the first solver to be a differential-algebraic equations solver
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_DAE_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"ODE Solver",err,error,*999)
            !Set solver defaults
            CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
            !Set the second solver to be a dynamic solver 
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_DYNAMIC_TYPE,err,error,*999)
            CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
            CALL Solver_LabelSet(solver,"Parabolic solver",err,error,*999)
            !Set solver defaults
            CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
            CALL Solver_DynamicSchemeSet(solver,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
            CALL Solver_DynamicRestartSet(solver,.TRUE.,err,error,*999)
            CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
          CASE(PROBLEM_MONODOMAIN_STRANG_SPLIT_SUBTYPE)
            CALL Solvers_NumberOfSolversSet(solvers,3,err,error,*999)
            !Set the first solver to be a differential-algebraic equations solver
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_DAE_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"First ODE solver",err,error,*999)
            !Set solver defaults
            CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
            !Set the second solver to be a dynamic solver 
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_DYNAMIC_TYPE,err,error,*999)
            CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
            CALL Solver_LabelSet(solver,"Parabolic solver",err,error,*999)
            !Set solver defaults
            CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
            CALL Solver_DynamicSchemeSet(solver,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
            CALL Solver_DynamicRestartSet(solver,.TRUE.,err,error,*999)
            CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
            !Set the third solver to be a differential-algebraic equations solver
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,3,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_DAE_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"Second ODE solver",err,error,*999)
            !Set solver defaults
            CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
          CASE(PROBLEM_MONODOMAIN_DIRECT_MODEL_SUBTYPE)
            CALL Solvers_NumberOfSolversSet(solvers,1,err,error,*999)
            !Set the solver to be a dynamic solver 
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_DYNAMIC_TYPE,err,error,*999)
            CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
            CALL Solver_LabelSet(solver,"Parabolic solver",err,error,*999)
            !Set solver defaults
            CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
            CALL Solver_DynamicSchemeSet(solver,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
            CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)           
          CASE DEFAULT
            localError="The problem subtype of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
              & " is invalid for a monodomain problem type of a bioelectric problem class."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_BIDOMAIN_EQUATION_TYPE)
          SELECT CASE(pSpecification(3))
          CASE(PROBLEM_BIDOMAIN_GUDUNOV_SPLIT_SUBTYPE)
            CALL Solvers_NumberOfSolversSet(solvers,3,err,error,*999)
            !Set the first solver to be a differential-algebraic equations solver
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_DAE_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"ODE solver",err,error,*999)
            !Set solver defaults
            CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
            !Set the second solver to be a dynamic solver 
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_DYNAMIC_TYPE,err,error,*999)
            CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
            CALL Solver_LabelSet(solver,"Parabolic solver",err,error,*999)
            !Set solver defaults
            CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
            CALL Solver_DynamicSchemeSet(solver,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
            CALL Solver_DynamicRestartSet(solver,.TRUE.,err,error,*999)
            CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
            !Set the third solver to be a linear solver
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,3,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_LINEAR_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"Elliptic solver",err,error,*999)
            !Set solver defaults
            CALL Solver_LibraryTypeSet(solver,SOLVER_PETSC_LIBRARY,err,error,*999)
          CASE(PROBLEM_BIDOMAIN_STRANG_SPLIT_SUBTYPE)
            CALL Solvers_NumberOfSolversSet(solvers,4,err,error,*999)
            !Set the first solver to be a differential-algebraic equations solver
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_DAE_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"First ODE solver",err,error,*999)
            !Set solver defaults
            CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
            !Set the second solver to be a dynamic solver 
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_DYNAMIC_TYPE,err,error,*999)
            CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
            CALL Solver_LabelSet(solver,"Parabolic solver",err,error,*999)
            !Set solver defaults
            CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
            CALL Solver_DynamicSchemeSet(solver,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
            CALL Solver_DynamicRestartSet(solver,.TRUE.,err,error,*999)
            CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
            !Set the third solver to be a differential-algebraic equations solver
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,3,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_DAE_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"Second ODE solver",err,error,*999)
            !Set solver defaults
            CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
            !Set the fourth solver to be a linear solver
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,4,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_LINEAR_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"Elliptic solver",err,error,*999)
            !Set solver defaults
            CALL Solver_LibraryTypeSet(solver,SOLVER_PETSC_LIBRARY,err,error,*999)
          CASE DEFAULT
            localError="The problem subtype of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
              & " is invalid for a monodomain problem type of a bioelectric problem class."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The problem type of "//TRIM(NumberToVString(pSpecification(2),"*",err,error))// &
            & " is invalid for a bioelectric problem class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Get the solvers
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
        !Finish the solvers creation
        CALL Solvers_CreateFinish(solvers,err,error,*999)
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a bioelectric equation."
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
        SELECT CASE(pSpecification(2))
        CASE(PROBLEM_MONODOMAIN_EQUATION_TYPE)
          NULLIFY(solver)
          SELECT CASE(pSpecification(3))
          CASE(PROBLEM_MONODOMAIN_GUDUNOV_SPLIT_SUBTYPE, &
            & PROBLEM_MONODOMAIN_STRANG_SPLIT_SUBTYPE)
            !Create the solver equations for the second (parabolic) solver
            CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
          CASE(PROBLEM_MONODOMAIN_DIRECT_MODEL_SUBTYPE)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          CASE DEFAULT
            localError="The problem subtype of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
              & " is invalid for a monodomain problem type of a bioelectric problem class."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          NULLIFY(solverEquations)
          CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
          CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_LINEAR,err,error,*999)
          CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
          CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
        CASE(PROBLEM_BIDOMAIN_EQUATION_TYPE)
          !Create the solver equations for the second (parabolic) solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
          NULLIFY(solverEquations)
          CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
          CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_LINEAR,err,error,*999)
          CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
          CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
          !Create the solver equations for the elliptic solver
          NULLIFY(solver)
          SELECT CASE(pSpecification(3))
          CASE(PROBLEM_BIDOMAIN_GUDUNOV_SPLIT_SUBTYPE)
            CALL Solvers_SolverGet(solvers,3,solver,err,error,*999)
          CASE(PROBLEM_BIDOMAIN_STRANG_SPLIT_SUBTYPE)
            CALL Solvers_SolverGet(solvers,4,solver,err,error,*999)
          CASE DEFAULT
            localError="The problem subtype of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))//  &
              & " is invalid for a bidomain problem type."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          NULLIFY(solverEquations)
          CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
          CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_LINEAR,err,error,*999)
          CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_STATIC,err,error,*999)
          CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
        CASE DEFAULT
          localError="The problem type of "//TRIM(NumberToVString(pSpecification(2),"*",err,error))//  &
            & " is invalid for a bioelectric problem class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Get the control loop
        NULLIFY(controlLoopRoot)
        CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
        NULLIFY(controlLoop)
        CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
        SELECT CASE(pSpecification(2))
        CASE(PROBLEM_MONODOMAIN_EQUATION_TYPE)
          NULLIFY(solver)
          SELECT CASE(pSpecification(3))
          CASE(PROBLEM_MONODOMAIN_GUDUNOV_SPLIT_SUBTYPE, &
            & PROBLEM_MONODOMAIN_STRANG_SPLIT_SUBTYPE)
            !Create the solver equations for the second (parabolic) solver
            CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
          CASE(PROBLEM_MONODOMAIN_DIRECT_MODEL_SUBTYPE)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          CASE DEFAULT
            localError="The problem subtype of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
              & " is invalid for a monodomain problem type of a bioelectric problem class."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          NULLIFY(solverEquations)
          CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
          !Finish the solver equations creation
          CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)             
        CASE(PROBLEM_BIDOMAIN_EQUATION_TYPE)
          !Get the solver equations for the second (parabolic) solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
          NULLIFY(solverEquations)
          CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
          !Finish the solver equations creation
          CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)             
          !Get the solver equations for the elliptic solver
          NULLIFY(solver)
          SELECT CASE(pSpecification(3))
          CASE(PROBLEM_BIDOMAIN_GUDUNOV_SPLIT_SUBTYPE)
            CALL Solvers_SolverGet(solvers,3,solver,err,error,*999)
          CASE(PROBLEM_BIDOMAIN_STRANG_SPLIT_SUBTYPE)
            CALL Solvers_SolverGet(solvers,4,solver,err,error,*999)
          CASE DEFAULT
            localError="The problem subtype of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))//  &
              & " is invalid for a bidomain problem type."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          NULLIFY(solverEquations)
          CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
          !Finish the solver equations creation
          CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)             
        CASE DEFAULT
          localError="The problem type of "//TRIM(NumberToVString(pSpecification(2),"*",err,error))//  &
            & " is invalid for a bioelectric problem class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a bioelectric equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_CELLML_EQUATIONS_TYPE)
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
        SELECT CASE(pSpecification(2))
        CASE(PROBLEM_MONODOMAIN_EQUATION_TYPE)
          !Create the CellML equations for the first DAE solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          NULLIFY(cellMLEquations)
          CALL CellMLEquations_CreateStart(solver,cellMLEquations,err,error,*999)
          !Set the time dependence
          CALL CellMLEquations_TimeDependenceTypeSet(cellMLEquations,CELLML_EQUATIONS_DYNAMIC,err,error,*999)
          !Set the linearity
          CALL CellMLEquations_LinearityTypeSet(cellMLEquations,CELLML_EQUATIONS_NONLINEAR,err,error,*999)
          IF(pSpecification(3)==PROBLEM_MONODOMAIN_STRANG_SPLIT_SUBTYPE) THEN
            !Create the CellML equations for the second DAE solver
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,3,solver,err,error,*999)
            NULLIFY(cellMLEquations)
            CALL CellMLEquations_CreateStart(solver,cellMLEquations,err,error,*999)
            !Set the time dependence
            CALL CellMLEquations_TimeDependenceTypeSet(cellMLEquations,CELLML_EQUATIONS_DYNAMIC,err,error,*999)
            !Set the linearity
            CALL CellMLEquations_LinearityTypeSet(cellMLEquations,CELLML_EQUATIONS_NONLINEAR,err,error,*999)
          ENDIF
        CASE(PROBLEM_BIDOMAIN_EQUATION_TYPE)
          !Create the CellML equations for the first DAE solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          NULLIFY(cellMLEquations)
          CALL CellMLEquations_CreateStart(solver,cellMLEquations,err,error,*999)
          !Set the time dependence
          CALL CellMLEquations_TimeDependenceTypeSet(cellMLEquations,CELLML_EQUATIONS_DYNAMIC,err,error,*999)
          !Set the linearity
          CALL CellMLEquations_LinearityTypeSet(cellMLEquations,CELLML_EQUATIONS_NONLINEAR,err,error,*999)
          IF(pSpecification(3)==PROBLEM_BIDOMAIN_STRANG_SPLIT_SUBTYPE) THEN
            !Create the CellML equations for the second DAE solver
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,4,solver,err,error,*999)
            NULLIFY(cellMLEquations)
            CALL CellMLEquations_CreateStart(SOLVER,cellMLEquations,err,error,*999)
            !Set the time dependence
            CALL CellMLEquations_TimeDependenceTypeSet(cellMLEquations,CELLML_EQUATIONS_DYNAMIC,err,error,*999)
            !Set the linearity
            CALL CellMLEquations_LinearityTypeSet(cellMLEquations,CELLML_EQUATIONS_NONLINEAR,err,error,*999)
          ENDIF
        CASE DEFAULT
          localError="The problem type of "//TRIM(NumberToVString(pSpecification(2),"*",err,error))//  &
            & " is invalid for a bioelectric problem class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Get the control loop
        NULLIFY(controlLoopRoot)
        CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
        NULLIFY(controlLoop)
        CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
        SELECT CASE(pSpecification(2))
        CASE(PROBLEM_MONODOMAIN_EQUATION_TYPE)
          !Get the CellML equations for the first DAE solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          NULLIFY(cellMLEquations)
          CALL Solver_CellMLEquationsGet(solver,cellMLEquations,err,error,*999)
          !Finish the CellML equations creation
          CALL CellMLEquations_CreateFinish(cellMLEquations,err,error,*999)
          IF(pSpecification(3)==PROBLEM_MONODOMAIN_STRANG_SPLIT_SUBTYPE) THEN
            !Get the CellML equations for the second DAE solver
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,3,solver,err,error,*999)
            NULLIFY(cellMLEquations)
            CALL Solver_CellMLEquationsGet(solver,cellMLEquations,err,error,*999)
            !Finish the CellML equations creation
            CALL CellMLEquations_CreateFinish(cellMLEquations,err,error,*999)
          ENDIF
        CASE(PROBLEM_BIDOMAIN_EQUATION_TYPE)
          !Get the CellML equations for the first DAE solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          NULLIFY(cellMLEquations)
          CALL Solver_CellMLEquationsGet(solver,cellMLEquations,err,error,*999)
          !Finish the CellML equations creation
          CALL CellMLEquations_CreateFinish(cellMLEquations,err,error,*999)
          IF(pSpecification(3)==PROBLEM_BIDOMAIN_STRANG_SPLIT_SUBTYPE) THEN
            !Get the CellML equations for the second DAE solver
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,4,solver,err,error,*999)
            NULLIFY(cellMLEquations)
            CALL Solver_CellMLEquationsGet(solver,cellMLEquations,err,error,*999)
            !Finish the CellML equations creation
            CALL CellMLEquations_CreateFinish(cellMLEquations,err,error,*999)
          ENDIF
        CASE DEFAULT
          localError="The problem type of "//TRIM(NumberToVString(pSpecification(2),"*",err,error))//  &
            & " is invalid for a bioelectric problem class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a bioelectric equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
        & " is invalid for a bioelectric domain equation."
      CALL FlagError(localError,err,error,*999)
    END SELECT
      
    EXITS("Biodomain_ProblemSetup")
    RETURN
999 ERRORSEXITS("Biodomain_ProblemSetup",err,error)
    RETURN 1
    
  END SUBROUTINE Biodomain_ProblemSetup
  
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

    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
    IF(ALLOCATED(problem%specification)) CALL FlagError("Problem specification is already allocated.",err,error,*999)
    IF(SIZE(problemSpecification,1)<3) THEN
      localError="The size of the specified problem specification array of "// &
        & TRIM(NumberToVString(SIZE(problemSpecification,1),"*",err,error))// &
        & " is invalid. The size should be >= 3."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
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
    
    ALLOCATE(problem%specification(3),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
    problem%specification(1:3)=[PROBLEM_BIOELECTRICS_CLASS,problemType,problemSubtype]

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
  SUBROUTINE Biodomain_FiniteElementCalculate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG), PARAMETER :: MAX_NUMBER_OF_COMPONENTS=3
    INTEGER(INTG) colsVariableType,columnComponentIdx,columnElementDOFIdx,columnElementParameterIdx,columnXiIdx,componentIdx, &
      & componentIdx2,esSpecification(3),gaussPointIdx,numberOfColsComponents, &
      & numberOfColumnElementParameters(MAX_NUMBER_OF_COMPONENTS),numberOfDimensions,numberOfGauss, &
      & numberOfRowElementParameters(MAX_NUMBER_OF_COMPONENTS), numberOfRowsComponents,numberOfXi,rowComponentIdx, &
      & rowElementDOFIdx,rowElementParameterIdx,rowXiIdx,rowsVariableType,scalingType,variableType,xiIdx
    REAL(DP) :: am,cm,columnPhi,columndPhidXi(MAX_NUMBER_OF_COMPONENTS), &
      & conductivity(MAX_NUMBER_OF_COMPONENTS,MAX_NUMBER_OF_COMPONENTS),Df,dPhidX(MAX_NUMBER_OF_COMPONENTS,64),Dt, &
      & extraConductivity(MAX_NUMBER_OF_COMPONENTS,MAX_NUMBER_OF_COMPONENTS),f(MAX_NUMBER_OF_COMPONENTS),fnorm, &
      & gaussWeight,jacobian,jacobianGaussWeight,intraConductivity(MAX_NUMBER_OF_COMPONENTS,MAX_NUMBER_OF_COMPONENTS), &
      & rowPhi,rowdPhidXi(MAX_NUMBER_OF_COMPONENTS),sourceParam,sum
    LOGICAL :: extracellular,uFibres,update,updateDamping,updateMatrices,updateMatrix,updateRHS,updateSource,updateStiffness, &
      & useFibre,uSource,vFibres    
    TYPE(BasisType), POINTER :: dependentBasis,geometricBasis,fibreBasis
    TYPE(BasisPtrType) :: columnBasis(MAX_NUMBER_OF_COMPONENTS),rowBasis(MAX_NUMBER_OF_COMPONENTS)
    TYPE(DecompositionType), POINTER :: dependentDecomposition,geometricDecomposition
    TYPE(DomainType), POINTER :: columnDomain,dependentDomain,geometricDomain,rowDomain
    TYPE(DomainElementsType), POINTER :: columnDomainElements,dependentDomainElements,geometricDomainElements,rowDomainElements
    TYPE(DomainTopologytype), POINTER  :: columnDomainTopology,dependentDomainTopology,geometricDomainTopology,rowDomainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(EquationsMappingSourceType), POINTER :: sourceMapping
    TYPE(EquationsMappingSourcesType), POINTER :: sourcesMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors
    TYPE(EquationsMatrixType), POINTER :: dampingMatrix,linearMatrix,stiffnessMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField,geometricField,fibreField,materialsField,sourceField
    TYPE(FieldInterpolationParametersType), POINTER :: colsInterpParameters,geometricInterpParameters,materialsInterpParameters, &
      & rowsInterpParameters,sourceInterpParameters,uFibreInterpParameters,vFibreInterpParameters
    TYPE(FieldInterpolatedPointType), POINTER :: geometricInterpPoint,materialsInterpPoint,sourceInterpPoint,uFibreInterpPoint, &
      & vFibreInterpPoint
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpPointMetrics
    TYPE(FieldVariableType), POINTER :: colsVariable,geometricVariable,rowsVariable
    TYPE(QuadratureSchemeType), POINTER :: dependentQuadratureScheme,geometricQuadratureScheme
    TYPE(QuadratureSchemePtrType) :: columnQuadratureScheme(MAX_NUMBER_OF_COMPONENTS),rowQuadratureScheme(MAX_NUMBER_OF_COMPONENTS)
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Biodomain_FiniteElementCalculate",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
    
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE)
      !OK
    CASE(EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE)
      SELECT CASE(esSpecification(3))
      CASE(EQUATIONS_SET_FIRST_BIDOMAIN_SUBTYPE, &
        & EQUATIONS_SET_SECOND_BIDOMAIN_SUBTYPE)
        !OK
      CASE DEFAULT
        localError="The equations set subtype of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
          & " is not valid for a bioelectric domain type of a bioelectrics equations set class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The equations set type of "//TRIM(NumberToVString(esSpecification(2),"*",err,error))// &
        & " is not valid for a bioelectric domain type of a bioelectrics equations set class."
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
    NULLIFY(rowsVariable)
    CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,rowsVariable,err,error,*999)
    NULLIFY(sourcesMapping)
    CALL EquationsMappingVector_SourcesMappingExists(vectorMapping,sourcesMapping,err,error,*999)
    IF(ASSOCIATED(sourcesMapping)) THEN
      CALL EquationsMappingSources_SourceMappingGet(sourcesMapping,1,sourceMapping,err,error,*999)
    ENDIF
    NULLIFY(rhsMapping)
    CALL EquationsMappingVector_RHSMappingExists(vectorMapping,rhsMapping,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(sourceVectors)
    NULLIFY(sourceVector)
    updateSource=.FALSE.
    IF(ASSOCIATED(sourcesMapping)) THEN
      CALL EquationsSet_SourceFieldGet(equationsSet,sourceField,err,error,*999)
      CALL EquationsMatricesVector_SourceVectorsGet(vectorMatrices,sourceVectors,err,error,*999)
      CALL EquationsMatricesSources_SourceVectorGet(sourceVectors,1,sourceVector,err,error,*999)
      CALL EquationsMatricesSource_UpdateVectorGet(sourceVector,updateSource,err,error,*999)
    ENDIF
    NULLIFY(rhsVector)
    updateRHS=.FALSE.
    IF(ASSOCIATED(rhsMapping)) THEN
      CALL EquationsMatricesVector_RHSVectorGet(vectorMatrices,rhsVector,err,error,*999)
      CALL EquationsMatricesRHS_UpdateVectorGet(rhsVector,updateRHS,err,error,*999)
      updateRHS=rhsVector%updateVector
    ENDIF
    NULLIFY(dynamicMapping)
    NULLIFY(linearMapping)    
    NULLIFY(dynamicMatrices)
    NULLIFY(linearMatrices)
    NULLIFY(stiffnessMatrix)
    NULLIFY(dampingMatrix)
    NULLIFY(linearMatrix)
    updateMatrices=.FALSE.
    updateStiffness=.FALSE.
    updateDamping=.FALSE.
    updateMatrix=.FALSE.
    extracellular=.FALSE.
    SELECT CASE(esSpecification(2))
    CASE(EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE)
      CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
      CALL EquationsMatricesVector_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*999)
      CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,1,stiffnessMatrix,err,error,*999)
      CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,2,dampingMatrix,err,error,*999)
      CALL EquationsMatrix_UpdateMatrixGet(stiffnessMatrix,updateStiffness,err,error,*999)
      CALL EquationsMatrix_UpdateMatrixGet(dampingMatrix,updateDamping,err,error,*999)
      updateMatrices=(updateStiffness.OR.updateDamping)
    CASE(EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE)
      SELECT CASE(esSpecification(3))
      CASE(EQUATIONS_SET_FIRST_BIDOMAIN_SUBTYPE)
        CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
        CALL EquationsMatricesVector_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*999)
        CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,1,stiffnessMatrix,err,error,*999)
        CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,2,dampingMatrix,err,error,*999)
        CALL EquationsMatrix_UpdateMatrixGet(stiffnessMatrix,updateStiffness,err,error,*999)
        CALL EquationsMatrix_UpdateMatrixGet(dampingMatrix,updateDamping,err,error,*999)
        updateMatrices=(updateStiffness.OR.updateDamping)
!!TODO also need a linear extracellular matrix?
      CASE(EQUATIONS_SET_SECOND_BIDOMAIN_SUBTYPE)
        extracellular=.TRUE.
        CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
        CALL EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*999)
        CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,1,linearMatrix,err,error,*999)
        CALL EquationsMatrix_UpdateMatrixGet(linearMatrix,updateMatrix,err,error,*999)
        updateMatrices=updateMatrix
      CASE DEFAULT
        localError="The equations set subtype of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
          & " is not valid for a bioelectric domain type of a bioelectrics equations set class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The equations set type of "//TRIM(NumberToVString(esSpecification(2),"*",err,error))// &
        & " is not valid for a bioelectric domain type of a bioelectrics equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    update=(updateMatrices.OR.updateSource.OR.updateRHS)

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
     
      NULLIFY(fibreField)
      CALL EquationsSet_FibreFieldExists(equationsSet,fibreField,err,error,*999)

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
      
      NULLIFY(rowsVariable)
      CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,rowsVariable,err,error,*999)
      CALL FieldVariable_VariableTypeGet(rowsVariable,rowsVariableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(rowsVariable,numberOfRowsComponents,err,error,*999)
      
      NULLIFY(colsVariable)
      IF(esSpecification(2)==EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE.AND. &
        & esSpecification(3)==EQUATIONS_SET_SECOND_BIDOMAIN_SUBTYPE) THEN
        CALL EquationsMappingLinear_LinearMatrixVariableGet(linearMapping,1,colsVariable,err,error,*999)
      ELSE
        CALL EquationsMappingDynamic_DynamicVariableGet(dynamicMapping,colsVariable,err,error,*999)
      ENDIF
      CALL FieldVariable_VariableTypeGet(colsVariable,colsVariableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(colsVariable,numberOfColsComponents,err,error,*999)

      NULLIFY(materialsField)
      CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
      
      NULLIFY(sourceField)
      IF(ASSOCIATED(sourceMapping)) CALL EquationsSet_SourceFieldExists(equationsSet,sourceField,err,error,*999)
      
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

      NULLIFY(materialsInterpParameters)
      CALL EquationsInterpolation_MaterialsParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & materialsInterpParameters,err,error,*999)
      NULLIFY(materialsInterpPoint)
      CALL EquationsInterpolation_MaterialsPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,materialsInterpPoint, &
        & err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,materialsInterpParameters,err,error,*999)
      
      NULLIFY(uFibreInterpParameters)
      NULLIFY(uFibreInterpPoint)
      NULLIFY(vFibreInterpParameters)
      NULLIFY(vFibreInterpPoint)
      uFibres=.FALSE.
      vFibres=.FALSE.
      IF(ASSOCIATED(fibreField)) THEN
        uFibres=.TRUE.
        CALL EquationsInterpolation_FibreParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,uFibreInterpParameters, &
          & err,error,*999)
        CALL EquationsInterpolation_FibrePointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,uFibreInterpPoint,err,error,*999)
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,uFibreInterpParameters,err,error,*999)
        IF(esSpecification(2)==EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE) THEN          
          vFibres=.TRUE.
          CALL EquationsInterpolation_FibreParametersGet(equationsInterpolation,FIELD_V_VARIABLE_TYPE,vFibreInterpParameters, &
            & err,error,*999)
          CALL EquationsInterpolation_FibrePointGet(equationsInterpolation,FIELD_V_VARIABLE_TYPE,vFibreInterpPoint,err,error,*999)
          CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,vFibreInterpParameters,err,error,*999)
        ENDIF
      ENDIF
      
      NULLIFY(sourceInterpParameters)
      NULLIFY(sourceInterpPoint)
      uSource=.FALSE.
      IF(ASSOCIATED(sourceField)) THEN
        uSource=.TRUE.
        CALL EquationsInterpolation_SourceParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,sourceInterpParameters, &
          & err,error,*999)
        CALL EquationsInterpolation_SourcePointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,sourceInterpPoint,err,error,*999)
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,sourceInterpParameters,err,error,*999)
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
      IF(updateMatrix) THEN
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
        CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,geometricInterpPointMetrics,err,error,*999)
        CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,materialsInterpPoint, &
          & err,error,*999)
        IF(uFibres) THEN
          CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,uFibreInterpPoint, &
            & err,error,*999)
          IF(vFibres) CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx, &
            & vFibreInterpPoint,err,error,*999)
        ENDIF
        IF(uSource) THEN
          CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,sourceInterpPoint, &
            & err,error,*999)
          sourceParam=sourceInterpPoint%values(1,NO_PART_DERIV)
        ENDIF

        IF((esSpecification(2)==EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE).AND. &
          & (esSpecification(3)==EQUATIONS_SET_MONODOMAIN_BUENOOROVIO_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_MONODOMAIN_TENTUSSCHER06_SUBTYPE)) THEN
          
          ! Diffusion tensor  D  =  Dt I + (Df - Dt) f f^T  where Dt and Df are diffusivity/conductivity in fiber/transverse
          ! directions
          Df = materialsInterpPoint%values(2,NO_PART_DERIV) ! 2 = Df
          Dt = materialsInterpPoint%values(3,NO_PART_DERIV) ! 3 = Dt
          fnorm = 0.0
          DO componentIdx=1,numberOfDimensions
            f(componentIdx)=materialsInterpPoint%values(3+componentIdx,NO_PART_DERIV) ! 4,5[,6] = f
            fnorm=fnorm+f(componentIdx)*f(componentIdx)
          ENDDO !componentIdx
          ! normalize f, and fill in default for 0,0,0 -> 1,0,0
          fnorm=SQRT(fnorm)
          IF(fnorm<ZERO_TOLERANCE) THEN
            f=[ 1.0, 0.0, 0.0 ] ! default
          ELSE
            f=f/fnorm
          ENDIF
          DO componentIdx=1,numberOfDimensions
            intraConductivity(componentIdx,:)=0.0
            intraConductivity(componentIdx,componentIdx)=Dt
            DO componentIdx2=1,numberOfDimensions
              intraConductivity(componentIdx,componentIdx2)=intraConductivity(componentIdx,componentIdx2)+ &
                & (Df-Dt)*f(componentIdx)*f(componentIdx2)
            ENDDO !componentIdx2
          ENDDO !componentIdx
        ELSE
          am=materialsInterpPoint%values(1,NO_PART_DERIV)
          cm=materialsInterpPoint%values(2,NO_PART_DERIV)
          IF((am*cm)<ZERO_TOLERANCE) THEN
            localError="The Am value of "//TRIM(NumberToVString(am,"*",err,error))//" and the value of Cm of "// &
              & TRIM(NumberToVString(cm,"*",err,error))//" for Gauss point number "// &
              & TRIM(NumberToVString(gaussPointIdx,"*",err,error))//" is invalid as their product <= zero."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          !Calculate (intracellular) conductivity tensor
          CALL CoordinateSystem_MaterialTransformSymTensor2(geometricInterpPointMetrics,uFibreInterpPoint, &
            & materialsInterpPoint%values(3:3+NUMBER_OF_VOIGT(numberOfDimensions),NO_PART_DERIV),intraConductivity, &
            & err,error,*999)
          intraConductivity(1:numberOfXi,1:numberOfXi)=intraConductivity(1:numberOfXi,1:numberOfXi)/(am*cm)
          IF(extracellular) THEN
            !Calculate extracellular conductivity tensor
            CALL CoordinateSystem_MaterialTransformSymTensor2(geometricInterpPointMetrics,vFibreInterpPoint, &
              & materialsInterpPoint%values(3+NUMBER_OF_VOIGT(numberOfDimensions):2+2*NUMBER_OF_VOIGT(numberOfDimensions), &
              & NO_PART_DERIV),extraConductivity,err,error,*999)
          ENDIF
        ENDIF
                
        !Calculate jacobianGaussWeight.
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
            IF(updateMatrices) THEN
              DO xiIdx=1,numberOfXi
                CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme(rowComponentIdx)%ptr,rowElementParameterIdx, &
                  & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx),gaussPointIdx,rowdPhidXi(xiIdx),err,error,*999)
              ENDDO !xiIdx
!!TODO: The cols variable will be the one mapped to the columns of the appropriate matrix. When bidomain is done then we will need to select either V_m or \phi_e dependening on what dynamic/linear matrices we have. For now just do monodomain.
              SELECT CASE(esSpecification(2))
              CASE(EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE)
                !Loop over element columns
                columnElementDOFIdx=0
                DO columnComponentIdx=1,numberOfColsComponents
                  DO columnElementParameterIdx=1,numberOfColumnElementParameters(columnComponentIdx)
                    columnElementDOFIdx=columnElementDOFIdx+1
                    CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme(columnComponentIdx)%ptr, &
                      & columnElementParameterIdx,NO_PART_DERIV,gaussPointIdx,columnPhi,err,error,*999)
                    IF(updateStiffness) THEN
                      DO xiIdx=1,numberOfXi
                        CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme(columnComponentIdx)%ptr, &
                          & columnElementParameterIdx,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx),gaussPointIdx, &
                          & columndPhidXi(xiIdx),err,error,*999)
                      ENDDO !xiIdx
                      sum=0.0_DP
                      DO rowXiIdx=1,numberOfXi
                        DO columnXiIdx=1,numberOfXi
                          DO xiIdx=1,numberOfXi
                            sum=sum+intraConductivity(rowXiIdx,columnXiIdx)*rowdPhidXi(xiIdx)*columndPhidXi(columnXiIdx)* &
                              & geometricInterpPointMetrics%gu(rowXiIdx,xiIdx)
                          ENDDO !xiIdx
                        ENDDO !columnXiIdx
                      ENDDO !rowXiIdx
                      stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                        & stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)+ &
                        & sum*jacobianGaussWeight                      
                    ENDIF
                    IF(updateDamping) THEN
                      sum=rowPhi*columnPhi
                      dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                        & dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)+ &
                        & sum*jacobianGaussWeight
                    ENDIF
                  ENDDO !columnElementParameterIdx
                ENDDO !columnComponentIdx                
              CASE(EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE)
                SELECT CASE(esSpecification(3))
                CASE(EQUATIONS_SET_FIRST_BIDOMAIN_SUBTYPE)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(EQUATIONS_SET_SECOND_BIDOMAIN_SUBTYPE)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  localError="The equations set subtype of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
                    & " is not valid for a bioelectric domain type of a bioelectrics equations set class."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE DEFAULT
                localError="The equations set type of "//TRIM(NumberToVString(esSpecification(2),"*",err,error))// &
                  & " is not valid for a bioelectric domain type of a bioelectrics equations set class."
                CALL FlagError(localError,err,error,*999)
              END SELECT              
            ENDIF !updateMatrices       
            IF(updateSource) THEN
              sourceVector%elementVector%vector(rowElementDOFIdx)=sourceVector%elementVector%vector(rowElementDOFIdx)+ &
                & sourceParam*rowPhi*jacobianGaussWeight               
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
            ENDIF !update matrix
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
        
    EXITS("Biodomain_FiniteElementCalculate")
    RETURN
999 ERRORSEXITS("Biodomain_FiniteElementCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE Biodomain_FiniteElementCalculate

  !
  !================================================================================================================================
  !
   
END MODULE BiodomainEquationsRoutines
