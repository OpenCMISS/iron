!> \file
!> \author Ting Yu
!> \brief This module set the boundary conditions for the given equation set
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
!> Contributor(s): Ting Yu, Chris Bradley
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

!> \defgroup OpenCMISS_BoundaryConditions OpenCMISS::Iron::BoundaryConditions
!> This module handles all boundary conditions routines.
MODULE BoundaryConditionsRoutines

  USE BaseRoutines
  USE BasisRoutines
  USE BasisAccessRoutines
  USE BoundaryConditionAccessRoutines
  USE CmissMPI
  USE ComputationRoutines
  USE ComputationAccessRoutines
  USE Constants
  USE CoordinateSystemRoutines
  USE CoordinateSystemAccessRoutines
  USE DecompositionAccessRoutines
  USE DistributedMatrixVector
  USE DistributedMatrixVectorAccessRoutines
  USE DomainMappings
  USE EquationsAccessRoutines
  USE EquationsMappingAccessRoutines
  USE EquationsMatricesAccessRoutines
  USE EquationsSetAccessRoutines
  USE InterfaceAccessRoutines
  USE InterfaceConditionAccessRoutines
  USE InterfaceEquationsAccessRoutines
  USE InterfaceMappingAccessRoutines
  USE FieldRoutines
  USE FieldAccessRoutines
  USE InputOutput
  USE ISO_VARYING_STRING
  USE Kinds
  USE MeshAccessRoutines
#ifndef NOMPIMOD
  USE MPI
#endif
  USE NodeRoutines
  USE RegionAccessRoutines
  USE SolverAccessRoutines
  USE SolverMappingAccessRoutines
  USE Strings
  USE Timer
  USE Types
  USE Lists
  USE LINKEDLIST_ROUTINES

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

  INTERFACE BoundaryConditions_AddConstant
    MODULE PROCEDURE BoundaryConditions_AddConstantField
    MODULE PROCEDURE BoundaryConditions_AddConstantVariable
  END INTERFACE BoundaryConditions_AddConstant

  !>Adds to the value of the specified local DOF and sets this as a boundary condition on the specified local DOF.
  INTERFACE BoundaryConditions_AddLocalDOF
    MODULE PROCEDURE BoundaryConditions_AddLocalDOF0
    MODULE PROCEDURE BoundaryConditions_AddLocalDOF1
  END INTERFACE BoundaryConditions_AddLocalDOF

  INTERFACE BoundaryConditions_AddElement
    MODULE PROCEDURE BoundaryConditions_AddElementField
    MODULE PROCEDURE BoundaryConditions_AddElementVariable
  END INTERFACE BoundaryConditions_AddElement

  INTERFACE BoundaryConditions_AddNode
    MODULE PROCEDURE BoundaryConditions_AddNodeField
    MODULE PROCEDURE BoundaryConditions_AddNodeVariable
  END INTERFACE BoundaryConditions_AddNode

  INTERFACE BoundaryConditions_SetConstant
    MODULE PROCEDURE BoundaryConditions_SetConstantField
    MODULE PROCEDURE BoundaryConditions_SetConstantVariable
  END INTERFACE BoundaryConditions_SetConstant

  !>Sets a boundary condition on the specified local DOF.
  INTERFACE BoundaryConditions_SetLocalDOF
    MODULE PROCEDURE BoundaryConditions_SetLocalDOF0
    MODULE PROCEDURE BoundaryConditions_SetLocalDOF1
  END INTERFACE BoundaryConditions_SetLocalDOF

  INTERFACE BoundaryConditions_SetElement
    MODULE PROCEDURE BoundaryConditions_SetElementField
    MODULE PROCEDURE BoundaryConditions_SetElementVariable
  END INTERFACE BoundaryConditions_SetElement

  INTERFACE BoundaryConditions_SetNode
    MODULE PROCEDURE BoundaryConditions_SetNodeField
    MODULE PROCEDURE BoundaryConditions_SetNodeVariable
  END INTERFACE BoundaryConditions_SetNode

  PUBLIC BoundaryConditions_AddConstant

  PUBLIC BoundaryConditions_AddLocalDOF

  PUBLIC BoundaryConditions_AddElement

  PUBLIC BoundaryConditions_AddNode

  PUBLIC BoundaryConditions_ConstrainNodeDOFsEqual

  PUBLIC BoundaryConditions_CreateFinish,BoundaryConditions_CreateStart

  PUBLIC BoundaryConditions_Destroy

  PUBLIC BoundaryConditions_LHSVariableExists

  PUBLIC BoundaryConditions_LHSVariableGet

  PUBLIC BoundaryConditions_NeumannSparsityTypeSet

  PUBLIC BoundaryConditions_RHSVariableExists

  PUBLIC BoundaryConditions_RHSVariableGet

  PUBLIC BoundaryConditions_RowVariableExists

  PUBLIC BoundaryConditions_RowVariableGet

  PUBLIC BoundaryConditions_SetAnalyticBoundaryNode

  PUBLIC BoundaryConditions_SetConstant

  PUBLIC BoundaryConditions_SetLocalDOF
  
  PUBLIC BoundaryConditions_SetElement

  PUBLIC BoundaryConditions_SetNode

  PUBLIC BoundaryConditions_UpdateAnalyticNode

  PUBLIC BoundaryConditions_VariableExists

  PUBLIC BoundaryConditions_VariableGet

  PUBLIC BoundaryConditionsRowVariable_RowConditionTypeSet

  PUBLIC BoundaryConditionsVariable_NeumannIntegrate

 CONTAINS

  !
  !================================================================================================================================
  !

  !>Adds to the value of the specified constant and sets this as a boundary condition on the specified constant. \see OpenCMISS::Iron::cmfe_BoundaryConditions_AddConstant
  SUBROUTINE BoundaryConditions_AddConstantField(boundaryConditions,field,variableType,componentNumber,condition,bcValue, &
    & err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FieldType), POINTER :: field !<The dependent field to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: variableType !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The component number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: condition !<The boundary condition type to set \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
    REAL(DP), INTENT(IN) :: bcValue !<The value of the boundary condition to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(FieldVariableType), POINTER :: fieldVariable

    ENTERS("BoundaryConditions_AddConstantField",err,error,*999)

    NULLIFY(fieldVariable)
    CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
    CALL BoundaryConditions_AddConstantVariable(boundaryConditions,fieldVariable,componentNumber,condition,bcValue,err,error,*999)

    EXITS("BoundaryConditions_AddConstantField")
    RETURN
999 ERRORSEXITS("BoundaryConditions_AddConstantField",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_AddConstantField

  !
  !================================================================================================================================
  !

  !>Adds to the value of the specified constant and sets this as a boundary condition on the specified constant.
  SUBROUTINE BoundaryConditions_AddConstantVariable(boundaryConditions,fieldVariable,componentNumber,condition,bcValue,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FieldVariableType), POINTER :: fieldVariable !<The dependent field variable to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The component number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: condition !<The boundary condition type to set \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
    REAL(DP), INTENT(IN) :: bcValue !<The value of the boundary condition to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: localDOF,globalDOF
    TYPE(VARYING_STRING) :: localError

    ENTERS("BoundaryConditions_AddConstantVariable",err,error,*999)

    !Note: This routine is for constant interpolation
    CALL BoundaryConditions_AssertNotFinished(boundaryConditions,err,error,*999)
    CALL FieldVariable_ConstantDOFGet(fieldVariable,componentNumber,localDOF,err,error,*999)
    CALL BoundaryConditions_CheckInterpolationType(condition,fieldVariable,componentNumber,err,error,*999)
    CALL BoundaryConditions_AddLocalDOF(boundaryConditions,fieldVariable,localDOF,condition,bcValue,err,error,*999)

    EXITS("BoundaryConditions_AddConstantVariable")
    RETURN
999 ERRORSEXITS("BoundaryConditions_AddConstantVariable",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_AddConstantVariable

  !
  !================================================================================================================================
  !

  !>Adds to the value of the specified constant and sets this as a boundary condition on the specified user element. \see OpenCMISS:Iron::BoundaryConditions::cmfe_BoundaryConditions_AddElement
  SUBROUTINE BoundaryConditions_AddElementField(boundaryConditions,field,variableType,userElementNumber,componentNumber, &
    & condition,bcValue,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FieldType), POINTER :: field !<The dependent field to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: variableType !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The component number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: condition !<The boundary condition type to set \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
    REAL(DP), INTENT(IN) :: bcValue !<The value of the boundary condition to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(FieldVariableType), POINTER :: fieldVariable

    ENTERS("BoundaryConditions_AddElementField",err,error,*999)

    !Note: this routine is for element based interpolation
    NULLIFY(fieldVariable)
    CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
    CALL BoundaryConditions_AddElementVariable(boundaryConditions,fieldVariable,userElementNumber,componentNumber, &
      & condition,bcValue,err,error,*999)

    EXITS("BoundaryConditions_AddElementField")
    RETURN
999 ERRORSEXITS("BoundaryConditions_AddElementField",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_AddElementField

  !
  !================================================================================================================================
  !

  !>Adds to the value of the specified constant and sets this as a boundary condition on the specified user element. 
  SUBROUTINE BoundaryConditions_AddElementVariable(boundaryConditions,fieldVariable,userElementNumber,componentNumber, &
    & condition,bcValue,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FieldVariableType), POINTER :: fieldVariable !<The dependent field variable to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The component number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: condition !<The boundary condition type to set \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
    REAL(DP), INTENT(IN) :: bcValue !<The value of the boundary condition to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: localDOF
    LOGICAL :: ghostDOF
    TYPE(VARYING_STRING) :: localError

    ENTERS("BoundaryConditions_AddElementVariable",err,error,*999)

    !Note: this routine is for element based interpolation
    CALL BoundaryConditions_AssertNotFinished(boundaryConditions,err,error,*999)
    CALL FieldVariable_UserElementDOFGet(fieldVariable,userElementNumber,componentNumber,localDOF,ghostDOF,err,error,*999)
    CALL BoundaryConditions_CheckInterpolationType(condition,fieldVariable,componentNumber,err,error,*999)
    CALL BoundaryConditions_AddLocalDOF(boundaryConditions,fieldVariable,localDOF,condition,bcValue,err,error,*999)

    EXITS("BoundaryConditions_AddElementVariable")
    RETURN
999 ERRORSEXITS("BoundaryConditions_AddElementVariable",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_AddElementVariable

  !
  !================================================================================================================================
  !

  !>Adds to the value of the specified DOF and sets this as a boundary condition on the specified DOF.
  SUBROUTINE BoundaryConditions_AddLocalDOF0(boundaryConditions,fieldVariable,dofIndex,condition,bcValue,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FieldVariableType), POINTER :: fieldVariable !<The dependent field variable to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: dofIndex !<The local dof index to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: condition !<The boundary condition type to set \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
    REAL(DP), INTENT(IN) :: bcValue !<The value of the boundary condition to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("BoundaryConditions_AddLocalDOF0",err,error,*999)

    CALL BoundaryConditions_AddLocalDOF1(boundaryConditions,fieldVariable,[dofIndex],[condition],[bcValue],err,error,*999)

    EXITS("BoundaryConditions_AddLocalDOF0")
    RETURN
999 ERRORSEXITS("BoundaryConditions_AddLocalDOF0",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_AddLocalDOF0

  !
  !================================================================================================================================
  !

  !>Adds to the value of the specified DOF and sets this as a boundary condition on the specified DOFs.
  SUBROUTINE BoundaryConditions_AddLocalDOF1(boundaryConditions,fieldVariable,dofIndices,conditions,bcValues,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FieldVariableType), POINTER :: fieldVariable !<The dependent field variable to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: dofIndices(:) !<dofIndices(:). The local dof index for the i'th dof to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: conditions(:) !<conditions(:). The boundary condition type to set for the i'th dof \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
    REAL(DP), INTENT(IN) :: bcValues(:) !<bcValues(:). The value of the boundary condition for the i'th dof to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dofIdx,globalDOF,localDOF,numberOfLocalDOFs
    REAL(DP) :: initialValue
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable
    TYPE(DomainMappingType), POINTER :: domainMapping
    TYPE(VARYING_STRING) :: localError

    ENTERS("BoundaryConditions_AddLocalDOF1",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(SIZE(dofIndices,1)/=SIZE(conditions,1)) THEN
      localError="The size of the DOF indices array of "//TRIM(NumberToVString(SIZE(dofIndices,1),"*",err,error))// &
        & " does not match the size of the fixed conditions array of "//TRIM(NumberToVString(SIZE(conditions,1),"*",err,error))// &
        & "."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(dofIndices,1)/=SIZE(bcValues,1)) THEN
      localError="The size of the DOF indices array of "//TRIM(NumberToVString(SIZE(dofIndices,1),"*",err,error))// &
        & " does not match the size of the values array of "//TRIM(NumberToVString(SIZE(bcValues,1),"*",err,error))//")."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    CALL BoundaryConditions_AssertNotFinished(boundaryConditions,err,error,*999)
    NULLIFY(boundaryConditionsVariable)
    CALL BoundaryConditions_VariableGet(boundaryConditions,fieldVariable,boundaryConditionsVariable,err,error,*999)
    NULLIFY(domainMapping)
    CALL FieldVariable_DomainMappingGet(fieldVariable,domainMapping,err,error,*999)
    CALL DomainMapping_NumberOfLocalGet(domainMapping,numberOfLocalDOFs,err,error,*999)
    DO dofIdx=1,SIZE(dofIndices,1)
      localDOF=dofIndices(dofIdx)
      IF(localDOF<1.OR.localDOF>numberOfLocalDOFs) THEN
        localError="The local DOF of  "//TRIM(NumberToVString(localDOF,"*",err,error))//" at DOF index "// &
          & TRIM(NumberToVString(dofIdx,"*",err,error))//" is invalid. The DOF should be >= 1 and <= "// &
          & TRIM(NumberToVString(numberOfLocalDOFs,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      CALL DomainMapping_LocalToGlobalGet(domainMapping,localDOF,globalDOF,err,error,*999)
      ! Set boundary condition and dof type, and make sure parameter sets are created
      CALL BoundaryConditionsVariable_ConditionTypeSet(boundaryConditionsVariable,globalDOF,conditions(dofIdx),err,error,*999)
      ! Update field sets by adding boundary condition values
      SELECT CASE(conditions(dofIdx))
      CASE(BOUNDARY_CONDITION_NONE)
        ! No field update
      CASE(BOUNDARY_CONDITION_FIXED)
        CALL FieldVariable_ParameterSetAddLocalDOF(fieldVariable,FIELD_VALUES_SET_TYPE,localDOF,bcValues(dofIdx),err,error,*999)
      CASE(BOUNDARY_CONDITION_FIXED_INLET)
        CALL FieldVariable_ParameterSetAddLocalDOF(fieldVariable,FIELD_VALUES_SET_TYPE,localDOF,bcValues(dofIdx),err,error,*999)
      CASE(BOUNDARY_CONDITION_FIXED_OUTLET)
        CALL FieldVariable_ParameterSetAddLocalDOF(fieldVariable,FIELD_VALUES_SET_TYPE,localDOF,bcValues(dofIdx),err,error,*999)
      CASE(BOUNDARY_CONDITION_FIXED_WALL)
        CALL FieldVariable_ParameterSetAddLocalDOF(fieldVariable,FIELD_VALUES_SET_TYPE,localDOF,bcValues(dofIdx),err,error,*999)
      CASE(BOUNDARY_CONDITION_MOVED_WALL)
        CALL FieldVariable_ParameterSetAddLocalDOF(fieldVariable,FIELD_VALUES_SET_TYPE,localDOF,bcValues(dofIdx),err,error,*999)
      CASE(BOUNDARY_CONDITION_FREE_WALL)
        CALL FieldVariable_ParameterSetAddLocalDOF(fieldVariable,FIELD_VALUES_SET_TYPE,localDOF,bcValues(dofIdx),err,error,*999)
      CASE(BOUNDARY_CONDITION_NEUMANN_POINT)
        ! Point value is stored in boundary conditions field set, and is then integrated to get the RHS variable value
        CALL FieldVariable_ParameterSetAddLocalDOF(fieldVariable,FIELD_BOUNDARY_CONDITIONS_SET_TYPE,localDOF,bcValues(dofIdx), &
          & err,error,*999)
      CASE(BOUNDARY_CONDITION_NEUMANN_INTEGRATED)
        ! For integrated Neumann condition, integration is already done, so set the RHS dof value directly
        CALL FieldVariable_ParameterSetAddLocalDOF(fieldVariable,FIELD_VALUES_SET_TYPE, &
          & localDOF,bcValues(dofIdx),err,error,*999)
      CASE(BOUNDARY_CONDITION_DIRICHLET)
         CALL FieldVariable_ParameterSetAddLocalDOF(fieldVariable,FIELD_VALUES_SET_TYPE,localDOF,bcValues(dofIdx),err,error,*999)
     CASE(BOUNDARY_CONDITION_CAUCHY)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(BOUNDARY_CONDITION_ROBIN)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(BOUNDARY_CONDITION_FIXED_INCREMENTED)
        ! For increment loops, we need to set the full BC parameter set value by
        ! getting the current value from the values parameter set
        CALL FieldVariable_ParameterSetGetLocalDOF(fieldVariable,FIELD_VALUES_SET_TYPE,localDOF,initialValue,err,error,*999)
        CALL FieldVariable_ParameterSetUpdateLocalDOF(fieldVariable,FIELD_BOUNDARY_CONDITIONS_SET_TYPE,localDOF, &
          & initialValue+bcValues(dofIdx),err,error,*999)
      CASE(BOUNDARY_CONDITION_PRESSURE)
        CALL FieldVariable_ParameterSetAddLocalDOF(fieldVariable,FIELD_PRESSURE_VALUES_SET_TYPE,localDOF,bcValues(dofIdx), &
          & err,error,*999)
      CASE(BOUNDARY_CONDITION_PRESSURE_INCREMENTED)
        ! For pressure incremented, adding to the values_set parameter value doesn't make sense,
        ! so just increment the value in the pressure values parameter set
        CALL FieldVariable_ParameterSetAddLocalDOF(fieldVariable,FIELD_PRESSURE_VALUES_SET_TYPE,localDOF,bcValues(dofIdx), &
          & err,error,*999)
      CASE(BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED)
        CALL FieldVariable_ParameterSetAddLocalDOF(fieldVariable,FIELD_VALUES_SET_TYPE,localDOF,bcValues(dofIdx),err,error,*999)
        CALL FieldVariable_ParameterSetAddLocalDOF(fieldVariable,FIELD_BOUNDARY_CONDITIONS_SET_TYPE,localDOF,bcValues(dofIdx), &
          & err,error,*999)
      CASE(BOUNDARY_CONDITION_CORRECTION_MASS_INCREASE)
        ! No field update
      CASE(BOUNDARY_CONDITION_IMPERMEABLE_WALL)
        CALL FieldVariable_ParameterSetAddLocalDOF(fieldVariable,FIELD_IMPERMEABLE_FLAG_VALUES_SET_TYPE,localDOF, &
          & bcValues(dofIdx),err,error,*999)
      CASE(BOUNDARY_CONDITION_NEUMANN_INTEGRATED_ONLY)
        ! For integrated Neumann condition, integration is already done, so set the RHS dof value directly
        CALL FieldVariable_ParameterSetAddLocalDOF(fieldVariable,FIELD_VALUES_SET_TYPE, &
          & localDOF,bcValues(dofIdx),err,error,*999)
      CASE(BOUNDARY_CONDITION_LINEAR_CONSTRAINT)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)
        CALL FieldVariable_ParameterSetAddLocalDOF(fieldVariable,FIELD_BOUNDARY_CONDITIONS_SET_TYPE,localDOF,bcValues(dofIdx), &
          & err,error,*999)
      CASE(BOUNDARY_CONDITION_FIXED_FITTED,BOUNDARY_CONDITION_FIXED_NONREFLECTING, &
        &  BOUNDARY_CONDITION_FIXED_CELLML,BOUNDARY_CONDITION_FIXED_STREE)
        CALL FieldVariable_ParameterSetAddLocalDOF(fieldVariable,FIELD_VALUES_SET_TYPE,localDOF,bcValues(dofIdx),err,error,*999)
      CASE(BOUNDARY_CONDITION_COUPLING_FLOW)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(BOUNDARY_CONDITION_COUPLING_STRESS)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(BOUNDARY_CONDITION_FIXED_PRESSURE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The specified boundary condition type for DOF index "//TRIM(NumberToVString(dofIdx,"*",err,error))//" of "// &
          & TRIM(NumberToVString(conditions(dofIdx),"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDDO !dofIdx

    EXITS("BoundaryConditions_AddLocalDOF1")
    RETURN
999 ERRORSEXITS("BoundaryConditions_AddLocalDOF1",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_AddLocalDOF1

  !
  !================================================================================================================================
  !

  !>Adds to the value of the specified constant and sets this as a boundary condition on the specified user node. \see OpenCMISS::Iron::cmfe_BoundaryConditions_AddNode
  SUBROUTINE BoundaryConditions_AddNodeField(boundaryConditions,field,variableType,versionNumber,derivativeNumber, &
    & userNodeNumber,componentNumber,condition,bcValue,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FieldType), POINTER :: field !<A pointer to the field to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: variableType !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: versionNumber !<The derivative version to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: derivativeNumber !<The derivative to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: userNodeNumber !<The user node number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The component number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: condition !<The boundary condition type to set \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
    REAL(DP), INTENT(IN) :: bcValue !<The value of the boundary condition to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(FieldVariableType), POINTER :: fieldVariable

    ENTERS("BoundaryConditions_AddNodeField",err,error,*999)

    NULLIFY(fieldVariable)
    CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
    CALL BoundaryConditions_AddNodeVariable(boundaryConditions,fieldVariable,versionNumber,derivativeNumber, &
      & userNodeNumber,componentNumber,condition,bcValue,err,error,*999)
    
    EXITS("BoundaryConditions_AddNodeField")
    RETURN
999 ERRORSEXITS("BoundaryConditions_AddNodeField",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_AddNodeField

  !
  !================================================================================================================================
  !

  !>Adds to the value of the specified constant and sets this as a boundary condition on the specified user node.
  SUBROUTINE BoundaryConditions_AddNodeVariable(boundaryConditions,fieldVariable,versionNumber,derivativeNumber, &
    & userNodeNumber,componentNumber,condition,bcValue,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: versionNumber !<The derivative version to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: derivativeNumber !<The derivative to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: userNodeNumber !<The user node number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The component number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: condition !<The boundary condition type to set \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
    REAL(DP), INTENT(IN) :: bcValue !<The value of the boundary condition to add
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: localDOF
    LOGICAL :: ghostDOF

    ENTERS("BoundaryConditions_AddNode",err,error,*999)

    CALL BoundaryConditions_AssertNotFinished(boundaryConditions,err,error,*999)
    CALL FieldVariable_UserNodeDOFGet(fieldVariable,versionNumber,derivativeNumber,userNodeNumber,componentNumber, &
      & localDOF,ghostDOF,err,error,*999)
    CALL BoundaryConditions_CheckInterpolationType(condition,fieldVariable,componentNumber,err,error,*999)
    CALL BoundaryConditions_AddLocalDOF(boundaryConditions,fieldVariable,localDOF,condition,bcValue,err,error,*999)
    
    EXITS("BoundaryConditions_AddNodeVariable")
    RETURN
999 ERRORSEXITS("BoundaryConditions_AddNodeVariable",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_AddNodeVariable

  !
  !================================================================================================================================
  !

  !>Checks that the specified boundary condition is appropriate for the field variable interpolation type
  SUBROUTINE BoundaryConditions_CheckInterpolationType(condition,fieldVariable,componentNumber,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: condition !<The boundary condition type being set
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable the boundary condition is set on
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The component number the boundary condition is set on
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: interpolationType
    LOGICAL :: validCondition

    ENTERS("BoundaryConditions_CheckInterpolationType",err,error,*999)

    CALL FieldVariable_ComponentInterpolationGet(fieldVariable,componentNumber,interpolationType,err,error,*999)

    validCondition=.TRUE.
    SELECT CASE(condition)
    CASE(BOUNDARY_CONDITION_NONE, &
      & BOUNDARY_CONDITION_FIXED, &
      & BOUNDARY_CONDITION_FIXED_INCREMENTED)
      !Valid for all interpolation types
    CASE(BOUNDARY_CONDITION_FIXED_INLET, &
      & BOUNDARY_CONDITION_FIXED_OUTLET)
      IF(interpolationType/=FIELD_NODE_BASED_INTERPOLATION) validCondition=.FALSE.
    CASE(BOUNDARY_CONDITION_FIXED_WALL, &
      & BOUNDARY_CONDITION_MOVED_WALL, &
      & BOUNDARY_CONDITION_FREE_WALL, &
      & BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED)
      IF(interpolationType/=FIELD_NODE_BASED_INTERPOLATION) validCondition=.FALSE.
    CASE(BOUNDARY_CONDITION_PRESSURE, &
      & BOUNDARY_CONDITION_PRESSURE_INCREMENTED)
      IF(interpolationType/=FIELD_NODE_BASED_INTERPOLATION) validCondition=.FALSE.
    CASE(BOUNDARY_CONDITION_CORRECTION_MASS_INCREASE)
      IF(interpolationType/=FIELD_NODE_BASED_INTERPOLATION) validCondition=.FALSE.
    CASE(BOUNDARY_CONDITION_IMPERMEABLE_WALL)
      IF(interpolationType/=FIELD_NODE_BASED_INTERPOLATION) validCondition=.FALSE.
    CASE(BOUNDARY_CONDITION_NEUMANN_POINT, &
      & BOUNDARY_CONDITION_NEUMANN_INTEGRATED, &
      & BOUNDARY_CONDITION_NEUMANN_INTEGRATED_ONLY, &
      & BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)
      IF(interpolationType/=FIELD_NODE_BASED_INTERPOLATION) validCondition=.FALSE.
    CASE(BOUNDARY_CONDITION_FIXED_FITTED,BOUNDARY_CONDITION_FIXED_NONREFLECTING,BOUNDARY_CONDITION_FIXED_CELLML, &
      & BOUNDARY_CONDITION_FIXED_STREE)
      IF(interpolationType/=FIELD_NODE_BASED_INTERPOLATION) validCondition=.FALSE.
    CASE DEFAULT
      CALL FlagError("The specified boundary condition type of "//TRIM(NumberToVString(condition,"*",err,error))//" is invalid.", &
        & err,error,*999)
    END SELECT
    IF(.NOT.validCondition) THEN
      CALL FlagError("The specified boundary condition type of "// &
        & TRIM(NumberToVString(condition,"*",err,error))//" is not valid for the field component "// &
        & "interpolation type of "//TRIM(NumberToVString(interpolationType,"*",err,error))//".", &
        & err,error,*999)
    ENDIF

    EXITS("BoundaryConditions_CheckInterpolationType")
    RETURN
999 ERRORSEXITS("BoundaryConditions_CheckInterpolationType",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_CheckInterpolationType

  !
  !================================================================================================================================
  !

  !>Constrain multiple equations dependent field DOFs to be a single solver DOF in the solver equations
  SUBROUTINE BoundaryConditions_ConstrainDOFsEqual(boundaryConditions,fieldVariable,globalDOFs,coefficient,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER, INTENT(IN) :: boundaryConditions !<The boundary conditions for the solver equations in which to constrain the DOF.
    TYPE(FieldVariableType), POINTER, INTENT(IN) :: fieldVariable !<A pointer to the field variable containing the DOFs.
    INTEGER(INTG), INTENT(IN) :: globalDOFs(:) !<The global DOFs to be constrained to be equal.
    REAL(DP), INTENT(IN) :: coefficient !<The coefficient of constraint.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message.
    !Local variables
    INTEGER(INTG) :: numberOfDOFs,dofIdx,dofIdx2
    TYPE(VARYING_STRING) :: localError

    ENTERS("BoundaryConditions_ConstrainDOFsEqual",err,error,*999)

    numberOfDOFs=SIZE(globalDOFs,1)
    IF(numberOfDOFs<2) CALL FlagError("Cannot constrain zero or 1 DOF to be equal.",err,error,*999)

    !Check for duplicate DOFs
    DO dofIdx=1,numberOfDOFs
      DO dofIdx2=dofIdx+1,numberOfDOFs
        IF(globalDOFs(dofIdx)==globalDOFs(dofIdx2)) THEN
          localError="DOF number "//TRIM(NumberToVstring(globalDOFs(dofIdx),"*",err,error))// &
            & " is duplicated in the DOFs constrained to be equal."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDDO !dofIdx2
    ENDDO !dofIdx

    !Add new DOF constraints
    !We set all DOFs except the first to be equal to coefficient * the first DOF
    !The first DOF is left unconstrained
    DO dofIdx=2,numberOfDOFs
      CALL BoundaryConditions_DOFConstraintSet(boundaryConditions,fieldVariable,globalDOFs(dofIdx),[globalDOFs(1)], &
        & [coefficient],err,error,*999)
    ENDDO !dofIdx

    EXITS("BoundaryConditions_ConstrainDOFsEqual")
    RETURN
999 ERRORSEXITS("BoundaryConditions_ConstrainDOFsEqual",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_ConstrainDOFsEqual

  !
  !================================================================================================================================
  !

  !>Constrain multiple nodal equations dependent field DOFs to be a single solver DOF in the solver equations
  SUBROUTINE BoundaryConditions_ConstrainNodeDOFsEqual(boundaryConditions,field,variableType,versionNumber,derivativeNumber, &
    & component,nodes,coefficient,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER, INTENT(IN) :: boundaryConditions !<The solver equations boundary conditions to constrain the DOFs for.
    TYPE(FieldType), POINTER, INTENT(IN) :: field !<The equations dependent field containing the field DOFs to be constrained.
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type of the DOFs to be constrained. \see OPENCMISS_FieldVariableTypes
    INTEGER(INTG), INTENT(IN) :: versionNumber !<The derivative version number.
    INTEGER(INTG), INTENT(IN) :: derivativeNumber !<The derivative number.
    INTEGER(INTG), INTENT(IN) :: component !<The field component number of the DOFs to be constrained.
    INTEGER(INTG), INTENT(IN) :: nodes(:) !<The user numbers of the nodes to be constrained to be equal.
    REAL(DP), INTENT(IN) :: coefficient !<The coefficient of constraint, applied to all but the first node.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message.
    !Local variables
    TYPE(FieldVariableType), POINTER :: fieldVariable
    INTEGER(INTG) :: numberOfNodes, nodeIdx, dof
    INTEGER(INTG), ALLOCATABLE :: globalDofs(:)

    ENTERS("BoundaryConditions_ConstrainNodeDOFsEqual",err,error,*998)

    IF(.NOT.ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions are not associated.",err,error,*998)
    !Get the field variable and boundary conditions variable for the field
    NULLIFY(fieldVariable)
    CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)

    numberOfNodes=SIZE(nodes,1)
    ALLOCATE(globalDOFs(numberOfNodes),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate equal global DOFs array.",err,error,*998)
    !Get field DOFs for nodes
    DO nodeIdx=1,numberOfNodes
      CALL Field_ComponentDOFGetUserNode(field,variableType,versionNumber,derivativeNumber,nodes(nodeIdx), &
        & component,dof,globalDofs(nodeIdx),err,error,*999)
    ENDDO !nodeIdx

    !Now set DOF constraint
    CALL BoundaryConditions_ConstrainDOFsEqual(boundaryConditions,fieldVariable,globalDofs,coefficient,err,error,*999)

    IF(ALLOCATED(globalDOFs)) DEALLOCATE(globalDOFs)

    EXITS("BoundaryConditions_ConstrainNodeDOFsEqual")
    RETURN
999 IF(ALLOCATED(globalDOFs)) DEALLOCATE(globalDOFs)
998 ERRORSEXITS("BoundaryConditions_ConstrainNodeDOFsEqual",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_ConstrainNodeDOFsEqual

  !
  !================================================================================================================================
  !

  !>Finish the creation of boundary conditions.
  SUBROUTINE BoundaryConditions_CreateFinish(boundaryConditions,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: boundaryConditionsVarType,colIdx,count,currentRowCondition,dirichletDOF,dirichletIdx,dofIdx,domainNumber, &
      & dummy,equationsMatrixIdx,equationsSetIdx,groupCommunicator,last,lhsVariableType,localDOFIdx,mpiIError, &
      & myGroupComputationNodeNumber,neumannIdx,numberOfBoundaryConditionsVariables,numberOfBoundaryConditionsRowVariables, &
      & numberOfDynamicMatrices,numberOfEquationsSets,numberOfGlobal,numberOfGroupComputationNodes,numberOfLinearMatrices, &
      & numberOfNonZeros,numberOfRows,parameterSetIdx,pressureIdx,rowIdx,sendCount,storageType,totalNumberOfLocal,variableIdx
    INTEGER(INTG), ALLOCATABLE:: columnArray(:)
    TYPE(BoundaryConditionsDirichletType), POINTER :: boundaryConditionsDirichlet
    TYPE(BoundaryConditionsNeumannType), POINTER :: boundaryConditionsNeumann
    TYPE(BoundaryConditionsPressureIncrementedType), POINTER :: boundaryConditionsPressureIncremented
    TYPE(BoundaryConditionsSparsityIndicesType), POINTER :: sparsityIndices
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable
    TYPE(BoundaryConditionsRowVariableType), POINTER :: boundaryConditionsRHSRowVariable,boundaryConditionsRowVariable
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix
    TYPE(DomainMappingType), POINTER :: domainMapping
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: field
    TYPE(FieldVariableType), POINTER :: fieldVariable,lhsVariable,matrixVariable
    TYPE(LinkedList), POINTER :: list(:)
    TYPE(ListType), POINTER :: sparseIndices
    TYPE(SolverType), POINTER :: solver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(VARYING_STRING) :: localError
    TYPE(WorkGroupType), POINTER :: workGroup

    ENTERS("BoundaryConditions_CreateFinish",err,error,*999)

    CALL BoundaryConditions_AssertNotFinished(boundaryConditions,err,error,*999)
    IF(.NOT.ALLOCATED(boundaryConditions%boundaryConditionsVariables)) &
      & CALL FlagError("Boundary conditions variables array is not allocated.",err,error,*999)
    
    NULLIFY(solverEquations)
    CALL BoundaryConditions_SolverEquationsGet(boundaryConditions,solverEquations,err,error,*999)
    NULLIFY(solver)
    CALL SolverEquations_SolverGet(solverEquations,solver,err,error,*999)
    NULLIFY(workGroup)
    CALL Solver_WorkGroupGet(solver,workGroup,err,error,*999)
    CALL WorkGroup_GroupCommunicatorGet(workGroup,groupCommunicator,err,error,*999)
    CALL WorkGroup_NumberOfGroupNodesGet(workGroup,numberOfGroupComputationNodes,err,error,*999)
    CALL WorkGroup_GroupNodeNumberGet(workGroup,myGroupComputationNodeNumber,err,error,*999)
    
    CALL BoundaryConditions_NumberOfVariablesGet(boundaryConditions,numberOfBoundaryConditionsVariables,err,error,*999)
    IF(numberOfGroupComputationNodes>0) THEN
      !Transfer all the boundary conditions to all the computation nodes.
      !\todo Look at this.
      DO variableIdx=1,numberOfBoundaryConditionsVariables
        NULLIFY(boundaryConditionsVariable)
        CALL BoundaryConditions_VariableIndexGet(boundaryConditions,variableIdx,boundaryConditionsVariable,err,error,*999)
        NULLIFY(fieldVariable)
        CALL BoundaryConditionsVariable_FieldVariableGet(boundaryConditionsVariable,fieldVariable,err,error,*999)
        NULLIFY(domainMapping)
        CALL FieldVariable_DomainMappingGet(fieldVariable,domainMapping,err,error,*999)
        CALL DomainMapping_NumberOfGlobalGet(domainMapping,sendCount,err,error,*999)
        NULLIFY(field)
        CALL FieldVariable_FieldGet(fieldVariable,field,err,error,*999)
        IF(numberOfGroupComputationNodes>1) THEN
          !\todo This operation is a little expensive as we are doing an unnecessary sum across all the ranks in order to combin
          !\todo the data from each rank into all ranks. We will see how this goes for now.
          CALL MPI_ALLREDUCE(MPI_IN_PLACE,boundaryConditionsVariable%DOFTypes,sendCount,MPI_INTEGER,MPI_SUM, &
            & groupCommunicator,mpiIError)
          CALL MPI_ErrorCheck("MPI_ALLREDUCE",mpiIError,err,error,*999)
          CALL MPI_ALLREDUCE(MPI_IN_PLACE,boundaryConditionsVariable%conditionTypes,sendCount,MPI_INTEGER,MPI_SUM, &
            & groupCommunicator,mpiIError)
          CALL MPI_ErrorCheck("MPI_ALLREDUCE",mpiIError,err,error,*999)
        ENDIF !mpi_in_place bug workaround - only do this when num comp nodes > 1
        IF(numberOfGroupComputationNodes>1) THEN
          ! Update the total number of boundary condition types by summing across all nodes
          CALL MPI_ALLREDUCE(MPI_IN_PLACE,boundaryConditionsVariable%dofCounts,MAX_BOUNDARY_CONDITION_NUMBER,MPI_INTEGER, &
            & MPI_SUM,groupCommunicator,mpiIError)
          CALL MPI_ErrorCheck("MPI_ALLREDUCE",mpiIError,err,error,*999)
          CALL MPI_ALLREDUCE(MPI_IN_PLACE,boundaryConditionsVariable%numberOfDirichletConditions,1,MPI_INTEGER,MPI_SUM, &
            & groupCommunicator,mpiIError)
          CALL MPI_ErrorCheck("MPI_ALLREDUCE",mpiIError,err,error,*999)
        ENDIF !mpi_in_place bug workaround - only do this when num comp nodes > 1
        
        !Check that the boundary conditions set are appropriate for equations sets
        CALL BoundaryConditionsVariable_CheckEquations(boundaryConditionsVariable,err,error,*999)

        IF(numberOfGroupComputationNodes>1) THEN
          !Make sure the required parameter sets are created on all computational nodes and begin updating them
          CALL MPI_ALLREDUCE(MPI_IN_PLACE,boundaryConditionsVariable%parameterSetRequired,FIELD_NUMBER_OF_SET_TYPES,MPI_LOGICAL, &
            & MPI_LOR,groupCommunicator,mpiIError)
          CALL MPI_ErrorCheck("MPI_ALLREDUCE",mpiIError,err,error,*999)
          DO parameterSetIdx=1,FIELD_NUMBER_OF_SET_TYPES
            IF(boundaryConditionsVariable%parameterSetRequired(parameterSetIdx)) THEN
              CALL FieldVariable_ParameterSetEnsureCreated(fieldVariable,parameterSetIdx,err,error,*999)
              CALL FieldVariable_ParameterSetUpdateStart(fieldVariable,parameterSetIdx,err,error,*999)
            ENDIF
          ENDDO
        ENDIF !mpi_in_place bug workaround - only do this when num comp nodes > 1
        
        !Set up pressure incremented condition, if it exists
        NULLIFY(boundaryConditionsPressureIncremented)
        IF(boundaryConditionsVariable%dofCounts(BOUNDARY_CONDITION_PRESSURE_INCREMENTED)>0) THEN
          CALL BoundaryConditionsVariable_PressureIncrementedInitialise(boundaryConditionsVariable,err,error,*999)
          CALL BoundaryConditionsVariable_PressureIncConditionsGet(boundaryConditionsVariable, &
            & boundaryConditionsPressureIncremented,err,error,*999)
        ENDIF

        !Set up Neumann condition information if there are any Neumann conditions
        NULLIFY(boundaryConditionsNeumann)        
        IF(boundaryConditionsVariable%dofCounts(BOUNDARY_CONDITION_NEUMANN_POINT)>0.OR. &
          & boundaryConditionsVariable%dofCounts(BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)>0) THEN
          CALL BoundaryConditionsVariable_NeumannInitialise(boundaryConditionsVariable,err,error,*999)
          CALL BoundaryConditionsVariable_NeumannConditionsGet(boundaryConditionsVariable,boundaryConditionsNeumann,err,error,*999)
        ENDIF

        NULLIFY(boundaryConditionsRowVariable)
        NULLIFY(boundaryConditionsRHSRowVariable)
        IF(boundaryConditionsVariable%dofCounts(BOUNDARY_CONDITION_PRESSURE_INCREMENTED)>0.OR. &
          & boundaryConditionsVariable%dofCounts(BOUNDARY_CONDITION_NEUMANN_POINT)>0.OR. &
          & boundaryConditionsVariable%dofCounts(BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)>0) THEN
          CALL BoundaryConditions_RHSVariableExists(boundaryConditions,fieldVariable,boundaryConditionsRHSRowVariable, &
            & err,error,*999)
          CALL BoundaryConditions_LHSVariableExists(boundaryConditions,fieldVariable,boundaryConditionsRowVariable, &
            & err,error,*999)
        ENDIF
        
        !Loop over all global DOFs, keeping track of the dof indices of specific BC types where required
        pressureIdx=1
        neumannIdx=1
        DO dofIdx=1,fieldVariable%numberOfGlobalDofs
          IF(boundaryConditionsVariable%conditionTypes(dofIdx)==BOUNDARY_CONDITION_PRESSURE_INCREMENTED) THEN
            boundaryConditionsPressureIncremented%pressureIncrementedDOFIndices(pressureIdx)=dofIdx
            pressureIdx=pressureIdx+1
            IF(ASSOCIATED(boundaryConditionsRHSRowVariable)) THEN
              CALL DomainMapping_DomainNumberFromGlobalGet(domainMapping,dofIdx,1,domainNumber,err,error,*999)
              IF(domainNumber==myGroupComputationNodeNumber) THEN
                !On my rank
                CALL DomainMapping_LocalNumberFromGlobalGet(domainMapping,dofIdx,1,localDOFIdx,err,error,*999)
                CALL BoundaryConditionsRowVariable_RowConditionTypeGet(boundaryConditionsRHSRowVariable,localDOFIdx, &
                  & currentRowCondition,err,error,*999)
                IF(currentRowCondition==BOUNDARY_CONDITION_FREE_ROW) & 
                  & CALL BoundaryConditionsRowVariable_RowConditionTypeSet(boundaryConditionsRHSRowVariable,localDOFIdx, &
                  & BOUNDARY_CONDITION_NEUMANN_ROW,err,error,*999)
              ENDIF
            ENDIF
            IF(ASSOCIATED(boundaryConditionsRowVariable)) THEN
              CALL DomainMapping_DomainNumberFromGlobalGet(domainMapping,dofIdx,1,domainNumber,err,error,*999)
              IF(domainNumber==myGroupComputationNodeNumber) THEN
                !On my rank
                CALL DomainMapping_LocalNumberFromGlobalGet(domainMapping,dofIdx,1,localDOFIdx,err,error,*999)
                CALL BoundaryConditionsRowVariable_RowConditionTypeGet(boundaryConditionsRHSRowVariable,localDOFIdx, &
                  & currentRowCondition,err,error,*999)
                IF(currentRowCondition==BOUNDARY_CONDITION_FREE_ROW) & 
                  CALL BoundaryConditionsRowVariable_RowConditionTypeSet(boundaryConditionsRowVariable,localDOFIdx, &
                  & BOUNDARY_CONDITION_NEUMANN_ROW,err,error,*999)
              ENDIF
            ENDIF
          ELSE IF(boundaryConditionsVariable%conditionTypes(dofIdx)==BOUNDARY_CONDITION_NEUMANN_POINT.OR. &
            & boundaryConditionsVariable%conditionTypes(dofIdx)==BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED) THEN
            boundaryConditionsNeumann%setDofs(neumannIdx)=dofIdx
            neumannIdx=neumannIdx+1
            IF(ASSOCIATED(boundaryConditionsRHSRowVariable)) THEN
              CALL DomainMapping_DomainNumberFromGlobalGet(domainMapping,dofIdx,1,domainNumber,err,error,*999)
              IF(domainNumber==myGroupComputationNodeNumber) THEN
                !On my rank
                CALL DomainMapping_LocalNumberFromGlobalGet(domainMapping,dofIdx,1,localDOFIdx,err,error,*999)
                CALL BoundaryConditionsRowVariable_RowConditionTypeSet(boundaryConditionsRHSRowVariable,localDOFIdx, &
                  & BOUNDARY_CONDITION_NEUMANN_ROW,err,error,*999)
              ENDIF
            ENDIF
            IF(ASSOCIATED(boundaryConditionsRowVariable)) THEN
              CALL DomainMapping_DomainNumberFromGlobalGet(domainMapping,dofIdx,1,domainNumber,err,error,*999)
              IF(domainNumber==myGroupComputationNodeNumber) THEN
                !On my rank
                CALL DomainMapping_LocalNumberFromGlobalGet(domainMapping,dofIdx,1,localDOFIdx,err,error,*999)
                CALL BoundaryConditionsRowVariable_RowConditionTypeSet(boundaryConditionsRowVariable,localDOFIdx, &
                  & BOUNDARY_CONDITION_NEUMANN_ROW,err,error,*999)
              ENDIF
            ENDIF
          ENDIF
        ENDDO !dofIdx

        !Now that we know where Neumann point DOFs are, we can calculate matrix structure
        IF(boundaryConditionsVariable%dofCounts(BOUNDARY_CONDITION_NEUMANN_POINT)>0.OR. &
          & boundaryConditionsVariable%dofCounts(BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)>0) THEN
          CALL BoundaryConditionsVariable_NeumannMatricesInitialise(boundaryConditionsVariable,err,error,*999)
        ENDIF
        
        !Check that there is at least one dirichlet condition
        NULLIFY(boundaryConditionsDirichlet)
        NULLIFY(boundaryConditionsRowVariable)
        IF(boundaryConditionsVariable%numberOfDirichletConditions>0) THEN
          CALL BoundaryConditionsVariable_DirichletInitialise(boundaryConditionsVariable,err,error,*999)
          CALL BoundaryConditionsVariable_DirichletConditionsGet(boundaryConditionsVariable,boundaryConditionsDirichlet, &
            & err,error,*999)
          CALL BoundaryConditions_RowVariableExists(boundaryConditions,fieldVariable,boundaryConditionsRowVariable,err,error,*999)
          ! Find dirichlet conditions
          dirichletIdx=1
          DO dofIdx=1,fieldVariable%numberOfGlobalDofs
            IF(boundaryConditionsVariable%DOFTypes(dofIdx)==BOUNDARY_CONDITION_DOF_FIXED) THEN
              boundaryConditionsDirichlet%dirichletDOFIndices(dirichletIdx)=dofIdx
              dirichletIdx=dirichletIdx+1
              IF(ASSOCIATED(boundaryConditionsRowVariable)) THEN
                CALL DomainMapping_DomainNumberFromGlobalGet(domainMapping,dofIdx,1,domainNumber,err,error,*999)
                IF(domainNumber==myGroupComputationNodeNumber) THEN
                  !On my rank
                  CALL DomainMapping_LocalNumberFromGlobalGet(domainMapping,dofIdx,1,localDOFIdx,err,error,*999)
                  CALL BoundaryConditionsRowVariable_RowConditionTypeSet(boundaryConditionsRowVariable,localDOFIdx, &
                    & BOUNDARY_CONDITION_DIRICHLET_ROW,err,error,*999)
                ENDIF
              ENDIF
            ENDIF
          ENDDO !dofIdx
            
          !Store Dirichlet dof indices
          NULLIFY(solverEquations)
          CALL BoundaryConditions_SolverEquationsGet(boundaryConditions,solverEquations,err,error,*999)
          NULLIFY(solverMapping)
          CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
          CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
          DO equationsSetIdx=1,numberOfEquationsSets
            NULLIFY(equationsSet)
            CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
            NULLIFY(equations)
            CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
            NULLIFY(vectorEquations)
            CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
            NULLIFY(vectorMatrices)
            CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
            CALL EquationsMatricesVector_TotalNumberOfRowsGet(vectorMatrices,numberOfRows,err,error,*999)
            NULLIFY(vectorMapping)
            CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
            NULLIFY(linearMapping)
            CALL EquationsMappingVector_LinearMappingExists(vectorMapping,linearMapping,err,error,*999)
            IF(ASSOCIATED(linearMapping)) THEN
              NULLIFY(linearMatrices)
              CALL EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*999)
              !Iterate through equations matrices
              CALL EquationsMatricesLinear_NumberOfLinearMatricesGet(linearMatrices,numberOfLinearMatrices,err,error,*999)
              DO equationsMatrixIdx=1,numberOfLinearMatrices
                NULLIFY(equationsMatrix)
                CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,equationsMatrixIdx,equationsMatrix,err,error,*999)
                NULLIFY(matrixVariable)
                CALL EquationsMappingLinear_LinearMatrixVariableGet(linearMapping,equationsMatrixIdx,matrixVariable,err,error,*999)
                NULLIFY(distributedMatrix)
                CALL EquationsMatrix_DistributedMatrixGet(equationsMatrix,distributedMatrix,err,error,*999)
                IF(ASSOCIATED(matrixVariable,fieldVariable)) THEN
                  IF(boundaryConditionsVariable%numberOfDirichletConditions>0) THEN
                    CALL DistributedMatrix_TransposeRowsColumnsSet(distributedMatrix,boundaryConditionsDirichlet% &
                      & dirichletDOFIndices(1:boundaryConditionsVariable%numberOfDirichletConditions),err,error,*999)
                  ENDIF
                  CALL DistributedMatrix_StorageTypeGet(distributedMatrix,storageType,err,error,*999)
                  SELECT CASE(storageType)
                  CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
                    !Do nothing
                  CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
                    !Do nothing
                  CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
                    CALL FlagError("Not implemented for column major storage.",err,error,*999)
                  CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
                    CALL FlagError("Not implemented for row major storage.",err,error,*999)
                  CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
                    !Get Sparsity pattern, number of non zeros, number of rows
                    !CALL DistributedMatrix_StorageLocationsGet(distributedMatrix,rowIndices,columnIndices,err,error,*999)
                    CALL DistributedMatrix_NumberOfNonZerosGet(distributedMatrix,numberOfNonZeros,err,error,*999)
                    !Get the matrix stored as a linked list
                    NULLIFY(list)
                    CALL DistributedMatrix_LinkListGet(distributedMatrix,list,err,error,*999)
                    !Initialise sparsity indices arrays
                    CALL BoundaryConditionsSparsityIndices_Initialise(boundaryConditionsDirichlet% &
                      & linearSparsityIndices(equationsSetIdx,equationsMatrixIdx)%ptr, &
                      & boundaryConditionsVariable%numberOfDirichletConditions,err,error,*999)
                    !Find dirichlet columns and store the non zero indices (with respect to the 1D storage array)
                    NULLIFY(sparsityIndices)
                    CALL BoundaryConditionsDirichlet_LinearSparsityIndicesGet(boundaryConditionsDirichlet,equationsMatrixIdx, &
                      & equationsSetIdx,sparsityIndices,err,error,*999)
                    !Setup list for storing dirichlet non zero indices
                    NULLIFY(sparseIndices)
                    CALL List_CreateStart(sparseIndices,err,error,*999)
                    CALL List_DataTypeSet(sparseIndices,LIST_INTG_TYPE,err,error,*999)
                    CALL List_InitialSizeSet(sparseIndices,boundaryConditionsVariable%numberOfDirichletConditions* &
                      & numberOfNonZeros/numberOfRows,err,error,*999)
                    CALL List_CreateFinish(sparseIndices,err,error,*999)
                    count=0
                    sparsityIndices%sparseColumnIndices(1)=1
                    last=1
                    DO dirichletIdx=1,boundaryConditionsVariable%numberOfDirichletConditions
                      dirichletDOF=boundaryConditionsDirichlet%dirichletDOFIndices(dirichletIdx)
                      CALL LinkedList_to_Array(list(dirichletDOF),columnArray,err,error,*999)
                      DO rowIdx=1,SIZE(columnArray)
                        CALL List_ItemAdd(sparseIndices,columnArray(rowIdx),err,error,*999)
                        count=count+1
                        last=rowIdx+17
                      ENDDO !rowIdx
                      sparsityIndices%sparseColumnIndices(dirichletIdx+1)=count+1
                    ENDDO !dirichletIdx
                    CALL List_DetachAndDestroy(sparseIndices,dummy,sparsityIndices%sparseRowIndices,err,error,*999)
                    DO colIdx =1,numberOfRows
                      CALL LinkedList_Destroy(list(colIdx),err,error,*999)
                    ENDDO !colIdx
                  CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
                    CALL FlagError("Not implemented for compressed column storage.",err,error,*999)
                  CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
                    CALL FlagError("Not implemented for row column storage.",err,error,*999)
                  CASE DEFAULT
                    localError="The storage type of "//TRIM(NumberToVString(storageType,"*",err,error))//" is invalid."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                ENDIF
              ENDDO !equationMatrixIdx
            ENDIF

            NULLIFY(dynamicMapping)
            CALL EquationsMappingVector_DynamicMappingExists(vectorMapping,dynamicMapping,err,error,*999)
            IF(ASSOCIATED(dynamicMapping)) THEN
              NULLIFY(dynamicMatrices)
              CALL EquationsMatricesVector_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*999)
              !Iterate through equations matrices
              CALL EquationsMatricesDynamic_NumberOfDynamicMatricesGet(dynamicMatrices,numberOfDynamicMatrices,err,error,*999)
              DO equationsMatrixIdx=1,numberOfDynamicMatrices
                NULLIFY(equationsMatrix)
                CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,equationsMatrixIdx,equationsMatrix,err,error,*999)
                NULLIFY(matrixVariable)
                CALL EquationsMappingDynamic_DynamicVariableGet(dynamicMapping,matrixVariable,err,error,*999)
                NULLIFY(distributedMatrix)
                CALL EquationsMatrix_DistributedMatrixGet(equationsMatrix,distributedMatrix,err,error,*999)
                IF(ASSOCIATED(matrixVariable,fieldVariable)) THEN
                  IF(boundaryConditionsVariable%numberOfDirichletConditions>0) THEN
                    CALL DistributedMatrix_TransposeRowsColumnsSet(distributedMatrix,boundaryConditionsDirichlet% &
                      & dirichletDOFIndices(1:boundaryConditionsVariable%numberOfDirichletConditions),err,error,*999)
                  ENDIF
                  CALL DistributedMatrix_StorageTypeGet(distributedMatrix,storageType,err,error,*999)
                  SELECT CASE(storageType)
                  CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
                    !Do nothing
                  CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
                    !Do nothing
                  CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
                    CALL FlagError("Not implemented for column major storage.",err,error,*999)
                  CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
                    CALL FlagError("Not implemented for row major storage.",err,error,*999)
                  CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
                    !Get Sparsity pattern, number of non zeros, number of rows
                    !CALL DistributedMatrix_StorageLocationsGet(distributedMatrix,rowIndices,columnIndices,err,error,*999)
                    CALL DistributedMatrix_NumberOfNonZerosGet(distributedMatrix,numberOfNonZeros,err,error,*999)
                    !Sparse matrix in a list
                    CALL DistributedMatrix_LinkListGet(distributedMatrix,list,err,error,*999)
                    !Intialise sparsity indices arrays
                    CALL BoundaryConditionsSparsityIndices_Initialise(boundaryConditionsDirichlet% &
                      & dynamicSparsityIndices(equationsSetIdx,equationsMatrixIdx)%ptr, &
                      & boundaryConditionsVariable%numberOfDirichletConditions,err,error,*999)
                    !Find dirichlet columns and store the non zero indices (with respect to the 1D storage array)
                    NULLIFY(sparsityIndices)
                    CALL BoundaryConditionsDirichlet_DynamicSparsityIndicesGet(boundaryConditionsDirichlet,equationsMatrixIdx, &
                      & equationsSetIdx,sparsityIndices,err,error,*999)
                    ! Setup list for storing dirichlet non zero indices
                    NULLIFY(sparseIndices)
                    CALL List_CreateStart(sparseIndices,err,error,*999)
                    CALL List_DataTypeSet(sparseIndices,LIST_INTG_TYPE,err,error,*999)
                    CALL List_InitialSizeSet(sparseIndices,boundaryConditionsVariable%numberOfDirichletConditions* &
                      & numberOfNonZeros/numberOfRows,err,error,*999)
                    CALL List_CreateFinish(sparseIndices,err,error,*999)
                    count=0
                    sparsityIndices%sparseColumnIndices(1)=1
                    last=1
                    DO dirichletIdx=1,boundaryConditionsVariable%numberOfDirichletConditions
                      !Dirichlet columns
                      dirichletDOF=boundaryConditionsDirichlet%dirichletDOFIndices(dirichletIdx)
                      CALL LinkedList_to_Array(list(dirichletDOF),columnArray,err,error,*999)
                      !The row indices
                      DO rowIdx=1,SIZE(columnArray)
                        CALL List_ItemAdd(sparseIndices,columnArray(rowIdx),err,error,*999)
                        count=count+1
                        last=rowIdx+1
                      ENDDO !rowIdx
                      sparsityIndices%sparseColumnIndices(dirichletIdx+1)=count+1
                    ENDDO !dirichletIdx
                    CALL List_DetachAndDestroy(sparseIndices,dummy,sparsityIndices%sparseRowIndices,err,error,*999)
                    DO colIdx =1,numberOfRows
                      CALL LinkedList_Destroy(list(colIdx),err,error,*999)
                    ENDDO !colIdx
                  CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
                    CALL FlagError("Not implemented for compressed column storage.",err,error,*999)
                  CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
                    CALL FlagError("Not implemented for row column storage.",err,error,*999)
                  CASE DEFAULT
                    localError="The storage type of "//TRIM(NumberToVString(storageType,"*",err,error)) &
                      //" is invalid."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                ENDIF
              ENDDO !equationsMatrixIdx
            ENDIF
          ENDDO !equationsSetIdx
         
          !!\todo Update interface sparsity structure calculate first then update code below.
          !!Loop over interface conditions. Note that only linear interface matrices implemented so far.
          !CALL SolverMapping_NumberOfIterfaceConditionsGet(solverMapping,numberOfInterfaceConditions,err,error,*999)
          !DO interfaceConditionIdx=1,numberOfInterfaceConditions
          !  NULLIFY(interfaceCondition)
          !  CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*999)
          !  NULLIFY(interfaceEquations)
          !  CALL InterfaceCondition_InterfaceEquationsGet(interfaceCondition,interfaceEquations,err,error,*999)
          !  NULLIFY(interfaceMatrices)
          !  CALL InterfaceEquations_InterfaceMatricesGet(interfaceEquations,interfaceMatrices,err,error,*999)
          !  !Iterate through equations matrices
          !  CALL InterfaceMatrices_NumberOfInterfaceMatricesGet(interfaceMatrices,numberOfInterfaceMatrices,err,error,*999)
          !  DO interfaceMatrixIdx=1,numberOfInterfaceMatrices
          !    NULLIFY(interfaceMatrix)
          !    CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,interfaceMatrixIdx,interfaceMatrix,err,error,*999)
          !    NULLIFY(distributedMatrix)
          !    CALL InterfaceMatrix_DistributedMatrixGet(interfaceMatrix,distributedMatrix,err,error,*999)
          !    CALL DistributedMatrix_StorageTypeGet(distributedMatrix,storageType,err,error,*999)
          !    SELECT CASE(storageType)
          !    CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
          !      !Do nothing
          !    CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
          !      !Do nothing
          !    CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
          !      CALL FlagError("Not implemented for column major storage.",err,error,*999)
          !    CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
          !      CALL FlagError("Not implemented for row major storage.",err,error,*999)
          !    CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
          !      !Get Sparsity pattern, number of non zeros, number of rows
          !      CALL DistributedMatrix_StorageLocationsGet(distributedMatrix,rowIndices,columnIndices,err,error,*999)
          !      CALL DistributedMatrix_NumberOfNonZerosGet(distributedMatrix,numberOfNonZeros,err,error,*999)
          !      !Get the matrix stored as a linked list
          !      CALL DistributedMatrix_LinkListGet(distributedMatrix,list,err,error,*999)
          !      CALL InterfaceMatrix_TotalNumberOfRowsGet(interfaceMatrix,numberOfRows,err,error,*999)
          !      !Initialise sparsity indices arrays
          !      CALL BoundaryConditions_SparsityIndicesInitialise(boundaryConditionsDirichlet% &
          !        & linearSparsityIndices(interfaceConditionIdx,interfaceMatrixIdx)%ptr, &
          !        & boundaryConditionsVariable%numberOfDirichletConditions,err,error,*999)
          !      !Find dirichlet columns and store the non zero indices (with respect to the 1D storage array
          !      NULLIFY(sparsityIndices)
          !      CALL BoundaryConditionsDirichlet_LinearSparsityIndicesGet(boundaryConditionsDirichlet,interfaceMatrixIdx, &
          !        & interfaceConditionIdx,sparsityIndices,err,error,*999)
          !      !Setup list for storing dirichlet non zero indices
          !      NULLIFY(sparseIndices)
          !      CALL List_CreateStart(sparseIndices,err,error,*999)
          !      CALL List_DataTypeSet(sparseIndices,LIST_INTG_TYPE,err,error,*999)
          !      CALL List_InitialSizeSet(sparseIndices,numberOfDirichletConditions*(numberOfNonZeros/numberOfRows),err,error,*999)
          !      CALL List_CreateFinish(sparseIndices,err,error,*999)
          !      count=0
          !      sparsityIndices%sparseColumnIndices(1)=1
          !      last=1
          !      DO dirichletIdx=1,boundaryConditionsVariable%numberOfDirichletConditions
          !        dirichletDOF=boundaryConditionsDirichlet%dirichletDOFIndices(dirichletIdx)
          !        CALL LinkedList_to_Array(list(dirichletDOF),columnArray)
          !        DO rowIdx=1,SIZE(columnArray)
          !          CALL List_ItemAdd(sparseIndices,columnArray(rowIdx),err,error,*999)
          !          count=count+1
          !          last=rowIdx+1
          !        ENDDO !rowIdx
          !        sparsityIndices%sparseColumnIndices(dirichletIdx+1)=count+1
          !      ENDDO !dirichletIdx
          !      CALL List_DetachAndDestroy(sparseIndices,dummy,sparsityIndices%sparseRowIndices,err,error,*999)
          !      DO colIdx =1,numberOfRows
          !        CALL LinkedList_Destroy(list(colIdx))
          !      ENDDO !colIdx
          !    CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
          !      CALL FlagError("Not implemented for compressed column storage.",err,error,*999)
          !    CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
          !      CALL FlagError("Not implemented for row column storage.",err,error,*999)
          !    CASE DEFAULT
          !      localError="The storage type of "//TRIM(NumberToVString(storageType,"*",err,error))//" is invalid."
          !      CALL FlagError(localError,err,error,*999)
          !    END SELECT
          !  ENDDO !interfaceMatrixIdx
          !ENDDO !interfaceConditionIdx
          
        ENDIF
        ! Finish field update
        DO parameterSetIdx=1,FIELD_NUMBER_OF_SET_TYPES
          IF(boundaryConditionsVariable%parameterSetRequired(parameterSetIdx)) THEN
            CALL FieldVariable_ParameterSetUpdateFinish(fieldVariable,parameterSetIdx,err,error,*999)
          ENDIF
        ENDDO !parameterSetIdx
        
        !Finish creating the boundary conditions DOF constraints
        CALL BoundaryConditionsVariable_DOFConstraintsCreateFinish(boundaryConditionsVariable,err,error,*999)
      ENDDO ! variableIdx
      
    ENDIF
    !Set the finished flag
    boundaryConditions%boundaryConditionsFinished=.TRUE.
    
    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Boundary conditions:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of boundary conditions row variables = ", &
        & numberOfBoundaryConditionsRowVariables,err,error,*999)
      DO variableIdx=1,numberOfBoundaryConditionsRowVariables
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Row variable index : ",variableIdx,err,error,*999)
        NULLIFY(boundaryConditionsRowVariable)
        CALL BoundaryConditions_RowVariableIndexGet(boundaryConditions,variableIdx,boundaryConditionsRowVariable,err,error,*999)
        NULLIFY(lhsVariable)
        CALL BoundaryConditionsRowVariable_LHSVariableGet(boundaryConditionsRowVariable,lhsVariable,err,error,*999)
        CALL FieldVariable_VariableTypeGet(lhsVariable,lhsVariableType,err,error,*999)
        CALL FieldVariable_TotalNumberOfDOFsGet(lhsVariable,totalNumberOfLocal,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Variable type = ",lhsVariableType,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Total number of local DOFs = ",totalNumberOfLocal,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,totalNumberOfLocal,8,8,boundaryConditionsRowVariable%rowConditionTypes, &
          & '("    Row BCs:",8(X,I8))','(12X,8(X,I8))',err,error,*999)
      ENDDO !variableIdx
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of boundary conditions variables = ", &
        & numberOfBoundaryConditionsVariables,err,error,*999)
      DO variableIdx=1,numberOfBoundaryConditionsVariables
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Variable index : ",variableIdx,err,error,*999)
        NULLIFY(boundaryConditionsVariable)
        CALL BoundaryConditions_VariableIndexGet(boundaryConditions,variableIdx,boundaryConditionsVariable,err,error,*999)
        CALL BoundaryConditionsVariable_FieldVariableTypeGet(boundaryConditionsVariable,boundaryConditionsVarType,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Variable type = ",boundaryConditionsVarType,err,error,*999)
        NULLIFY(fieldVariable)
        CALL BoundaryConditionsVariable_FieldVariableGet(boundaryConditionsVariable,fieldVariable,err,error,*999)
        NULLIFY(domainMapping)
        CALL FieldVariable_DomainMappingGet(fieldVariable,domainMapping,err,error,*999)
        CALL DomainMapping_NumberOfGlobalGet(domainMapping,numberOfGlobal,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of global dofs = ",numberOfGlobal,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfGlobal,8,8,boundaryConditionsVariable%conditionTypes, &
          & '("    Global BCs:",8(X,I8))','(15X,8(X,I8))',err,error,*999)
      ENDDO !variableIdx
    ENDIF
    
    EXITS("BoundaryConditions_CreateFinish")
    RETURN
999 ERRORSEXITS("BoundaryConditions_CreateFinish",err,error)
    RETURN 1

  END SUBROUTINE BoundaryConditions_CreateFinish

  !
  !================================================================================================================================
  !

  !>Start the creation of boundary conditions for the equation set.
  SUBROUTINE BoundaryConditions_CreateStart(solverEquations,boundaryConditions,err,error,*)

    !Argument variables
    TYPE(SolverEquationsType), POINTER :: solverEquations !<A pointer to the solver equations to create boundary conditions for
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<On exit, a pointer to the created boundary conditions. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("BoundaryConditions_CreateStart",err,error,*998)

#ifdef WITH_PRECHECKS
    IF(ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is not associated.",err,error,*999)
#endif
    
    !Initialise the boundary conditions
    CALL BoundaryConditions_Initialise(solverEquations,err,error,*999)
    !Return the pointer
    boundaryConditions=>solverEquations%boundaryConditions

    EXITS("BoundaryConditions_CreateStart")
    RETURN
999 NULLIFY(boundaryConditions)
998 ERRORSEXITS("BoundaryConditions_CreateStart",err,error)
    RETURN 1

  END SUBROUTINE BoundaryConditions_CreateStart

  !
  !================================================================================================================================
  !

  !>Destroys boundary conditions
  SUBROUTINE BoundaryConditions_Destroy(boundaryConditions,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("BoundaryConditions_Destroy",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions is not associated.",err,error,*999)
#endif
    
    CALL BoundaryConditions_Finalise(boundaryConditions,err,error,*999)

    EXITS("BoundaryConditions_Destroy")
    RETURN
999 ERRORSEXITS("BoundaryConditions_Destroy",err,error)
    RETURN 1

  END SUBROUTINE BoundaryConditions_Destroy

  !
  !================================================================================================================================
  !

  !>Constrain a DOF to be a linear combination of other DOFs.
  SUBROUTINE BoundaryConditions_DOFConstraintSet(boundaryConditions,fieldVariable,globalDOF,dofs,coefficients,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER, INTENT(IN) :: boundaryConditions !<The boundary conditions for the solver equations in which to constrain the DOF.
    TYPE(FieldVariableType), POINTER, INTENT(IN) :: fieldVariable !<A pointer to the field variable containing the DOFs.
    INTEGER(INTG), INTENT(IN) :: globalDOF !<The global DOF to set the constraint on.
    INTEGER(INTG), INTENT(IN) :: dofs(:) !<The global DOFs that this DOF is constrained to depend on.
    REAL(DP), INTENT(IN) :: coefficients(:) !<The coefficient values in the DOF constraint.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message.
    !Local variables
    INTEGER(INTG) :: numberOfDOFs,dofIdx,dofIdx2
    TYPE(BoundaryConditionsDofConstraintPtrType), ALLOCATABLE :: newConstraints(:)
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable
    TYPE(BoundaryConditionsDofConstraintsType), POINTER :: dofConstraints
    TYPE(BoundaryConditionsDofConstraintType), POINTER :: dofConstraint
    TYPE(VARYING_STRING) :: localError

    ENTERS("BoundaryConditions_DOFConstraintSet",err,error,*998)

    !Check pointers for association
    CALL BoundaryConditions_AssertNotFinished(boundaryConditions,err,error,*999)
    NULLIFY(boundaryConditionsVariable)
    CALL BoundaryConditions_VariableExists(boundaryConditions,fieldVariable,boundaryConditionsVariable,err,error,*998)
    NULLIFY(dofConstraints)
    CALL BoundaryConditionsVariable_DOFConstraintsGet(boundaryConditionsVariable,dofConstraints,err,error,*999)

    numberOfDOFs=SIZE(dofs,1)
    IF(numberOfDOFs==0) THEN
      CALL FlagError("Empty DOFs list.",err,error,*998)
    ELSE IF(numberOfDOFs/=SIZE(coefficients,1)) THEN
      CALL FlagError("Length of coefficients does not match length of DOFs array.",err,error,*998)
    ELSE IF(numberOfDOFs>1) THEN
      CALL FlagError("Support for constraining an equations DOF to be depended on multiple "// &
        & "other DOFs is not yet implemented.",err,error,*998)
    ENDIF

    !Check for duplicate DOFs
    DO dofIdx=1,numberOfDOFs
      DO dofIdx2=dofIdx+1,numberOfDOFs
        IF(dofs(dofIdx)==dofs(dofIdx2)) THEN
          localError="DOF number "//TRIM(NumberToVstring(dofs(dofIdx),"*",err,error))// &
            & " is duplicated in the DOF constraint."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDDO !dofIdx2
    ENDDO !dofIdx

    !Check DOFs are free
    DO dofIdx=1,numberOfDOFs
      IF(boundaryConditionsVariable%DOFTypes(dofs(dofIdx))/=BOUNDARY_CONDITION_DOF_FREE) THEN
        CALL FlagError("DOF number "//TRIM(NumberToVstring(dofs(dofIdx),"*",err,error))// &
          & " is not free in the boundary conditions.",err,error,*998)
      ENDIF
    ENDDO !dofIdx

    !Allocate new DOF constraints and copy over old constraints
    ALLOCATE(newConstraints(dofConstraints%numberOfConstraints+1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new DOF constraints array.",err,error,*998)
    IF(dofConstraints%numberOfConstraints>0) THEN
      newConstraints(1:dofConstraints%numberOfConstraints)= &
        & dofConstraints%constraints(1:dofConstraints%numberOfConstraints)
    ENDIF

    !Set the new DOF constraint
    ALLOCATE(dofConstraint,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new DOF constraint.",err,error,*999)
    ALLOCATE(dofConstraint%dofs(numberOfDofs),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate constraint DOFs array.",err,error,*999)
    ALLOCATE(dofConstraint%coefficients(numberOfDofs),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate constraint coefficients array.",err,error,*999)
    dofConstraint%globalDof=globalDOF
    dofConstraint%numberOfDofs=numberOfDOFs
    dofConstraint%dofs(1:numberOfDofs)=dofs(1:numberOfDofs)
    dofConstraint%coefficients(1:numberOfDofs)=coefficients(1:numberOfDofs)

    !Add new DOF constraint to new array
    newConstraints(dofConstraints%numberOfConstraints+1)%ptr=>dofConstraint
    !Replace old DOF constraints with new ones
    CALL MOVE_ALLOC(newConstraints,dofConstraints%constraints)
    dofConstraints%numberOfConstraints=dofConstraints%numberOfConstraints+1

    !Set the DOF type and BC type of the constrained DOF
    CALL BoundaryConditionsVariable_ConditionTypeSet(boundaryConditionsVariable,globalDOF,BOUNDARY_CONDITION_LINEAR_CONSTRAINT, &
      & err,error,*999)

    EXITS("BoundaryConditions_DOFConstraintSet")
    RETURN
999 IF(ASSOCIATED(dofConstraint)) THEN
      IF(ALLOCATED(dofConstraint%dofs)) DEALLOCATE(dofConstraint%dofs)
      IF(ALLOCATED(dofConstraint%coefficients)) DEALLOCATE(dofConstraint%coefficients)
      DEALLOCATE(dofConstraint)
    ENDIF
    IF(ALLOCATED(newConstraints)) DEALLOCATE(newConstraints)
998 ERRORSEXITS("BoundaryConditions_DOFConstraintSet",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_DOFConstraintSet

  !
  !================================================================================================================================
  !

  !>Finalise the boundary conditions and deallocate all memory.
  SUBROUTINE BoundaryConditions_Finalise(boundaryConditions,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: variableIdx

    ENTERS("BoundaryConditions_Finalise",err,error,*999)

    IF(ASSOCIATED(boundaryConditions)) THEN
      IF(ALLOCATED(boundaryConditions%boundaryConditionsRowVariables)) THEN
        DO variableIdx=1,SIZE(boundaryConditions%boundaryConditionsRowVariables,1)
          CALL BoundaryConditionsRowVariable_Finalise(boundaryConditions%boundaryConditionsRowVariables(variableIdx)%ptr, &
            & err,error,*999)
        ENDDO !variableIdx
        DEALLOCATE(boundaryConditions%boundaryConditionsRowVariables)
      ENDIF
      IF(ALLOCATED(boundaryConditions%boundaryConditionsVariables)) THEN
        DO variableIdx=1,SIZE(boundaryConditions%boundaryConditionsVariables,1)
          CALL BoundaryConditionsVariable_Finalise(boundaryConditions%boundaryConditionsVariables(variableIdx)%ptr,err,error,*999)
        ENDDO !variableIdx        
        DEALLOCATE(boundaryConditions%boundaryConditionsVariables)
      ENDIF
      IF(ASSOCIATED(boundaryConditions%solverEquations)) THEN
        IF(ASSOCIATED(boundaryConditions%solverEquations%solver)) &
          & NULLIFY(boundaryConditions%solverEquations%solver%solverEquations)
      ENDIF
      DEALLOCATE(boundaryConditions)
    ENDIF

    EXITS("BoundaryConditions_Finalise")
    RETURN
999 ERRORSEXITS("BoundaryConditions_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the boundary conditions for an equations set.
  SUBROUTINE BoundaryConditions_Initialise(solverEquations,err,error,*)

    !Argument variables
    TYPE(SolverEquationsType), POINTER :: solverEquations !<A pointer to the solver equations to initialise the boundary conditions for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,equationsSetIdx,interfaceConditionIdx,linearityType,numberOfEquationsMatrices, &
      & numberOfEquationsSets,numberOfInterfaceConditions,numberOfInterfaceMatrices,numberOfLinearVariables, &
      & numberOfResiduals,numberOfResidualVariables,numberOfVariables,residualIdx,timeDependenceType,variableIdx,variableType
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(FieldVariableType), POINTER :: dynamicVariable,lagrangeVariable,lhsVariable,linearVariable,residualVariable,rhsVariable
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping
    TYPE(InterfaceMappingRHSType), POINTER :: interfaceRHSMapping
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("BoundaryConditions_Initialise",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is not associated",err,error,*998)
    IF(ASSOCIATED(solverEquations%boundaryConditions)) &
      & CALL FlagError("Boundary conditions is already associated for these solver equations.",err,error,*998)
#endif    
    
    ALLOCATE(solverEquations%boundaryConditions,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate boundary conditions.",err,error,*999)
    solverEquations%boundaryConditions%boundaryConditionsFinished=.FALSE.
    solverEquations%boundaryConditions%numberOfBoundaryConditionsRowVariables=0
    solverEquations%boundaryConditions%numberOfBoundaryConditionsVariables=0
    solverEquations%boundaryConditions%solverEquations=>solverEquations
    solverEquations%boundaryConditions%neumannMatrixSparsity=BOUNDARY_CONDITION_SPARSE_MATRICES
    NULLIFY(solverMapping)
    CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
    CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
    DO equationsSetIdx=1,numberOfEquationsSets
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
      NULLIFY(equations)
      CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
      CALL Equations_AssertIsFinished(equations,err,error,*999)
      NULLIFY(vectorEquations)
      CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
      NULLIFY(vectorMapping)
      CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
      CALL EquationsMappingVector_AssertIsFinished(vectorMapping,err,error,*999)
!!TODO: We shouldn't need the equations set to know about boundary conditions???
      equationsSet%boundaryConditions=>solverEquations%boundaryConditions
      CALL Equations_TimeDependenceTypeGet(equations,timeDependenceType,err,error,*999)
      CALL Equations_LinearityTypeGet(equations,linearityType,err,error,*999)
      SELECT CASE(timeDependenceType)
      CASE(EQUATIONS_STATIC,EQUATIONS_QUASISTATIC)
        SELECT CASE(linearityType)
        CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
          NULLIFY(linearMapping)
          CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
          CALL EquationsMappingLinear_NumberOfLinearVariablesGet(linearMapping,numberOfLinearVariables,err,error,*999)
          DO variableIdx=1,numberOfLinearVariables
            CALL EquationsMappingLinear_VariableNumberOfMatricesGet(linearMapping,variableIdx,numberOfEquationsMatrices, &
              & err,error,*999)
            IF(numberOfEquationsMatrices>0) THEN                
              NULLIFY(linearVariable)
              CALL EquationsMappingLinear_LinearVariableGet(linearMapping,variableIdx,linearVariable,err,error,*999)
              CALL BoundaryConditions_VariableAdd(solverEquations%boundaryConditions,linearVariable,err,error,*999)
            ENDIF
          ENDDO !variableIdx
          NULLIFY(rhsMapping)
          NULLIFY(rhsVariable)
          CALL EquationsMappingVector_RHSMappingExists(vectorMapping,rhsMapping,err,error,*999)
          IF(ASSOCIATED(rhsMapping)) THEN
            CALL EquationsMappingRHS_RHSVariableGet(rhsMapping,rhsVariable,err,error,*999)
            CALL BoundaryConditions_VariableAdd(solverEquations%boundaryConditions,rhsVariable,err,error,*999)
          ENDIF
        CASE(EQUATIONS_NONLINEAR)
          NULLIFY(nonlinearMapping)
          CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
          CALL EquationsMappingNonlinear_NumberOfResidualsGet(nonlinearMapping,numberOfResiduals,err,error,*999)
          DO residualIdx=1,numberOfResiduals
            NULLIFY(residualMapping)
            CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,residualIdx,residualMapping,err,error,*999)
            CALL EquationsMappingResidual_NumberOfResidualVariablesGet(residualMapping,numberOfResidualVariables,err,error,*999)
            DO variableIdx=1,numberOfResidualVariables
              NULLIFY(residualVariable)
              CALL EquationsMappingResidual_VariableGet(residualMapping,variableIdx,residualVariable,err,error,*999)
              CALL BoundaryConditions_VariableAdd(solverEquations%boundaryConditions,residualVariable,err,error,*999)
            ENDDO !variableIdx
          ENDDO !residualIdx
          NULLIFY(linearMapping)
          CALL EquationsMappingVector_LinearMappingExists(vectorMapping,linearMapping,err,error,*999)
          IF(ASSOCIATED(linearMapping)) THEN
            CALL EquationsMappingLinear_NumberOfLinearVariablesGet(linearMapping,numberOfLinearVariables,err,error,*999)
            DO variableIdx=1,numberOfLinearVariables
              CALL EquationsMappingLinear_VariableNumberOfMatricesGet(linearMapping,variableIdx,numberOfEquationsMatrices, &
                & err,error,*999)
              IF(numberOfEquationsMatrices>0) THEN                
                NULLIFY(linearVariable)
                CALL EquationsMappingLinear_LinearVariableGet(linearMapping,variableIdx,linearVariable,err,error,*999)
                CALL BoundaryConditions_VariableAdd(solverEquations%boundaryConditions,linearVariable,err,error,*999)
              ENDIF
            ENDDO !variableIdx
          ENDIF
          NULLIFY(rhsMapping)
          NULLIFY(rhsVariable)
          CALL EquationsMappingVector_RHSMappingExists(vectorMapping,rhsMapping,err,error,*999)
          IF(ASSOCIATED(rhsMapping)) THEN
            CALL EquationsMappingRHS_RHSVariableGet(rhsMapping,rhsVariable,err,error,*999)
            CALL BoundaryConditions_VariableAdd(solverEquations%boundaryConditions,rhsVariable,err,error,*999)
          ENDIF
        CASE DEFAULT
          localError="The equations linearity type of "//TRIM(NumberToVString(linearityType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_FIRST_ORDER_DYNAMIC,EQUATIONS_SECOND_ORDER_DYNAMIC)
        SELECT CASE(linearityType)
        CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
          NULLIFY(dynamicMapping)
          CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
          NULLIFY(dynamicVariable)
          CALL EquationsMappingDynamic_DynamicVariableGet(dynamicMapping,dynamicVariable,err,error,*999)
          CALL BoundaryConditions_VariableAdd(solverEquations%boundaryConditions,dynamicVariable,err,error,*999)
          NULLIFY(linearMapping)
          CALL EquationsMappingVector_LinearMappingExists(vectorMapping,linearMapping,err,error,*999)
          IF(ASSOCIATED(linearMapping)) THEN
            CALL EquationsMappingLinear_NumberOfLinearVariablesGet(linearMapping,numberOfLinearVariables,err,error,*999)
            DO variableIdx=1,numberOfLinearVariables
              CALL EquationsMappingLinear_VariableNumberOfMatricesGet(linearMapping,variableIdx,numberOfEquationsMatrices, &
                & err,error,*999)
              IF(numberOfEquationsMatrices>0) THEN                
                NULLIFY(linearVariable)
                CALL EquationsMappingLinear_LinearVariableGet(linearMapping,variableIdx,linearVariable,err,error,*999)
                CALL BoundaryConditions_VariableAdd(solverEquations%boundaryConditions,linearVariable,err,error,*999)
              ENDIF
            ENDDO !variableIdx
          ENDIF
          NULLIFY(rhsMapping)
          NULLIFY(rhsVariable)
          CALL EquationsMappingVector_RHSMappingExists(vectorMapping,rhsMapping,err,error,*999)
          IF(ASSOCIATED(rhsMapping)) THEN
            CALL EquationsMappingRHS_RHSVariableGet(rhsMapping,rhsVariable,err,error,*999)
            CALL BoundaryConditions_VariableAdd(solverEquations%boundaryConditions,rhsVariable,err,error,*999)
          ENDIF
        CASE(EQUATIONS_NONLINEAR)
          NULLIFY(dynamicMapping)
          CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
          NULLIFY(dynamicVariable)
          CALL EquationsMappingDynamic_DynamicVariableGet(dynamicMapping,dynamicVariable,err,error,*999)
          CALL BoundaryConditions_VariableAdd(solverEquations%boundaryConditions,dynamicVariable,err,error,*999)
          NULLIFY(nonlinearMapping)
          CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
          CALL EquationsMappingNonlinear_NumberOfResidualsGet(nonlinearMapping,numberOfResiduals,err,error,*999)
          DO residualIdx=1,numberOfResiduals
            NULLIFY(residualMapping)
            CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,residualIdx,residualMapping,err,error,*999)
            CALL EquationsMappingResidual_NumberOfResidualVariablesGet(residualMapping,numberOfResidualVariables,err,error,*999)
            DO variableIdx=1,numberOfResidualVariables
              NULLIFY(residualVariable)
              CALL EquationsMappingResidual_VariableGet(residualMapping,variableIdx,residualVariable,err,error,*999)
              CALL BoundaryConditions_VariableAdd(solverEquations%boundaryConditions,residualVariable,err,error,*999)
            ENDDO !variableIdx
          ENDDO !residualIdx
          NULLIFY(linearMapping)
          CALL EquationsMappingVector_LinearMappingExists(vectorMapping,linearMapping,err,error,*999)
          IF(ASSOCIATED(linearMapping)) THEN
            CALL EquationsMappingLinear_NumberOfLinearVariablesGet(linearMapping,numberOfLinearVariables,err,error,*999)
            DO variableIdx=1,numberOfLinearVariables
              CALL EquationsMappingLinear_VariableNumberOfMatricesGet(linearMapping,variableIdx,numberOfEquationsMatrices, &
                & err,error,*999)
              IF(numberOfEquationsMatrices>0) THEN                
                NULLIFY(linearVariable)
                CALL EquationsMappingLinear_LinearVariableGet(linearMapping,variableIdx,linearVariable,err,error,*999)
                CALL BoundaryConditions_VariableAdd(solverEquations%boundaryConditions,linearVariable,err,error,*999)
              ENDIF
            ENDDO !variableIdx
          ENDIF
          NULLIFY(rhsMapping)
          NULLIFY(rhsVariable)
          CALL EquationsMappingVector_RHSMappingExists(vectorMapping,rhsMapping,err,error,*999)
          IF(ASSOCIATED(rhsMapping)) THEN
            CALL EquationsMappingRHS_RHSVariableGet(rhsMapping,rhsVariable,err,error,*999)
            CALL BoundaryConditions_VariableAdd(solverEquations%boundaryConditions,rhsVariable,err,error,*999)
          ENDIF
        CASE DEFAULT
          localError="The equations linearity type of "//TRIM(NumberToVString(linearityType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The equations time dependence type of "//TRIM(NumberToVString(timeDependenceType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Add the row variable
      NULLIFY(lhsMapping)
      CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
      NULLIFY(lhsVariable)
      CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,lhsVariable,err,error,*999)
      CALL BoundaryConditions_RowVariableAdd(solverEquations%boundaryConditions,lhsVariable,rhsVariable,err,error,*999)      
    ENDDO !equationsSetIdx
    CALL SolverMapping_NumberOfInterfaceConditionsGet(solverMapping,numberOfInterfaceConditions,err,error,*999)
    DO interfaceConditionIdx=1,numberOfInterfaceConditions
      NULLIFY(interfaceCondition)
      CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*999)
      NULLIFY(interfaceEquations)
      CALL InterfaceCondition_InterfaceEquationsGet(interfaceCondition,interfaceEquations,err,error,*999)
      CALL InterfaceEquations_AssertIsFinished(interfaceEquations,err,error,*999)
      NULLIFY(interfaceMapping)
      CALL InterfaceEquations_InterfaceMappingGet(interfaceEquations,interfaceMapping,err,error,*999)
      CALL InterfaceMapping_AssertIsFinished(interfaceMapping,err,error,*999)
!!TODO: We shouldn't need the interface condition to know about boundary conditions???
      interfaceCondition%boundaryConditions=>solverEquations%boundaryConditions
      CALL InterfaceEquations_TimeDependenceTypeGet(interfaceEquations,timeDependenceType,err,error,*999)
      CALL InterfaceEquations_LinearityTypeGet(interfaceEquations,linearityType,err,error,*999)
      !Only linear interface equations implemented at the moment
      SELECT CASE(timeDependenceType)
      CASE(INTERFACE_EQUATIONS_STATIC,INTERFACE_EQUATIONS_QUASISTATIC)
        SELECT CASE(linearityType)
        CASE(INTERFACE_EQUATIONS_LINEAR)
          NULLIFY(interfaceRHSMapping)
          NULLIFY(rhsVariable)
          CALL InterfaceMapping_RHSMappingExists(interfaceMapping,interfaceRHSMapping,err,error,*999)
          IF(ASSOCIATED(interfaceRHSMapping)) THEN
            CALL InterfaceMappingRHS_RHSVariableGet(interfaceRHSMapping,rhsVariable,err,error,*999)
            CALL BoundaryConditions_VariableAdd(solverEquations%boundaryConditions,rhsVariable,err,error,*999)
          ENDIF
          CALL InterfaceMapping_NumberOfInterfaceMatricesGet(interfaceMapping,numberOfInterfaceMatrices,err,error,*999)
          IF(numberOfInterfaceMatrices>0) THEN
            NULLIFY(lagrangeVariable)
            CALL InterfaceMapping_LagrangeVariableGet(interfaceMapping,lagrangeVariable,err,error,*999)
            CALL BoundaryConditions_VariableAdd(solverEquations%boundaryConditions,lagrangeVariable,err,error,*999)
            CALL BoundaryConditions_RowVariableAdd(solverEquations%boundaryConditions,lagrangeVariable,rhsVariable,err,error,*999)
          ENDIF
        CASE DEFAULT
          localError="The interface equations linearity type of "//TRIM(NumberToVString(linearityType,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The interface equations time dependence type of "//TRIM(NumberToVString(timeDependenceType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDDO !interfaceConditionIdx
    
    EXITS("BoundaryConditions_Initialise")
    RETURN
999 CALL BoundaryConditions_Finalise(solverEquations%boundaryConditions,dummyErr,dummyError,*998)
998 ERRORSEXITS("BoundaryConditions_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_Initialise

  !
  !================================================================================================================================
  !

  !>Find the boundary conditions row variable for a given field LHS variable if it exists
  SUBROUTINE BoundaryConditions_LHSVariableExists(boundaryConditions,lhsVariable,boundaryConditionsRowVariable,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to see if the boundary conditions LHS variable exists for
    TYPE(FieldVariableType), POINTER :: lhsVariable !<A pointer to the LHS variable to check the boundary conditions row variable for.
    TYPE(BoundaryConditionsRowVariableType), POINTER, INTENT(OUT) :: boundaryConditionsRowVariable !<On return, a pointer to the boundary conditions row variable, or NULL if it wasn't found. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: lhsInterfaceUserNumber,lhsVarType,lhsVariableFieldUserNumber,lhsVariableRegionUserNumber, &
      & matchInterfaceUserNumber,matchVariableType,matchVariableFieldUserNumber,matchVariableRegionUserNumber, &
      & numberOfBoundaryConditionsRowVariables,rowVariableIdx
    LOGICAL :: lhsVariableIsRegion,matchVariableIsInterface,matchVariableIsRegion,rowVariableFound
    TYPE(FieldType), POINTER :: lhsVariableField,matchVariableField
    TYPE(FieldVariableType), POINTER :: matchVariable
    TYPE(InterfaceType), POINTER :: lhsVariableInterface,matchVariableInterface
    TYPE(RegionType), POINTER :: lhsInterfaceParentRegion,lhsVariableRegion,matchInterfaceParentRegion,matchVariableRegion
    TYPE(VARYING_STRING) :: localError

    ENTERS("BoundaryConditions_LHSVariableExists",err,error,*998)

    IF(ASSOCIATED(boundaryConditionsRowVariable)) &
      & CALL FlagError("Boundary conditions row variable is already associated.",err,error,*998)

    CALL FieldVariable_VariableTypeGet(lhsVariable,lhsVarType,err,error,*999)
    NULLIFY(lhsVariableField)
    CALL FieldVariable_FieldGet(lhsVariable,lhsVariableField,err,error,*999)
    CALL Field_UserNumberGet(lhsVariableField,lhsVariableFieldUserNumber,err,error,*999)
    CALL Field_IsRegionField(lhsVariableField,lhsVariableIsRegion,err,error,*999)
    NULLIFY(lhsVariableRegion)
    NULLIFY(lhsVariableInterface)
    IF(lhsVariableIsRegion) THEN
      CALL Field_RegionGet(lhsVariableField,lhsVariableRegion,err,error,*999)
      CALL Region_UserNumberGet(lhsVariableRegion,lhsVariableRegionUserNumber,err,error,*999)
    ELSE
      CALL Field_InterfaceGet(lhsVariableField,lhsVariableInterface,err,error,*999)
      CALL Interface_UserNumberGet(lhsVariableInterface,lhsInterfaceUserNumber,err,error,*999)
      NULLIFY(lhsInterfaceParentRegion)
      CALL Interface_ParentRegionGet(lhsVariableInterface,lhsInterfaceParentRegion,err,error,*999)
      CALL Region_UserNumberGet(lhsInterfaceParentRegion,lhsVariableRegionUserNumber,err,error,*999)      
    ENDIF
        
    CALL BoundaryConditions_NumberOfRowVariablesGet(boundaryConditions,numberOfBoundaryConditionsRowVariables,err,error,*999)
    rowVariableFound=.FALSE.
    rowVariableIdx=1
    DO WHILE(rowVariableIdx<=numberOfBoundaryConditionsRowVariables.AND..NOT.rowVariableFound)
      NULLIFY(boundaryConditionsRowVariable)
      CALL BoundaryConditions_RowVariableIndexGet(boundaryConditions,rowVariableIdx,boundaryConditionsRowVariable,err,error,*999)
      NULLIFY(matchVariable)
      CALL BoundaryConditionsRowVariable_LHSVariableExists(boundaryConditionsRowVariable,matchVariable,err,error,*999)
      IF(ASSOCIATED(matchVariable)) THEN
        CALL FieldVariable_VariableTypeGet(matchVariable,matchVariableType,err,error,*999)
        NULLIFY(matchVariableField)
        CALL FieldVariable_FieldGet(matchVariable,matchVariableField,err,error,*999)
        CALL Field_UserNumberGet(matchVariableField,matchVariableFieldUserNumber,err,error,*999)
        CALL Field_IsRegionField(matchVariableField,matchVariableIsRegion,err,error,*999)
        IF(lhsVariableIsRegion) THEN
          IF(matchVariableIsRegion) THEN
            NULLIFY(matchVariableRegion)
            CALL Field_RegionGet(matchVariableField,matchVariableRegion,err,error,*999)
            CALL Region_UserNumberGet(matchVariableRegion,matchVariableRegionUserNumber,err,error,*999)
            IF(lhsVariableRegionUserNumber==matchVariableRegionUserNumber.AND. &
              & lhsVariableFieldUserNumber==matchVariableFieldUserNumber.AND. &
              & lhsVarType==matchVariableType) rowVariableFound=.TRUE.
          ENDIF
        ELSE
          CALL Field_IsInterfaceField(matchVariableField,matchVariableIsInterface,err,error,*999)
          IF(matchVariableIsInterface) THEN
            NULLIFY(matchVariableInterface)
            CALL Field_InterfaceGet(matchVariableField,matchVariableInterface,err,error,*999)
            CALL Interface_UserNumberGet(matchVariableInterface,matchInterfaceUserNumber,err,error,*999)
            NULLIFY(matchInterfaceParentRegion)
            CALL Interface_ParentRegionGet(matchVariableInterface,matchInterfaceParentRegion,err,error,*999)
            CALL Region_UserNumberGet(matchInterfaceParentRegion,matchVariableRegionUserNumber,err,error,*999)
            IF(lhsInterfaceUserNumber==matchInterfaceUserNumber.AND. &
              & lhsVariableRegionUserNumber==matchVariableRegionUserNumber.AND. &
              & lhsVariableFieldUserNumber==matchVariableFieldUserNumber.AND. &
              & lhsVarType==matchVariableType) rowVariableFound=.TRUE.
          ENDIF
        ENDIF
      ENDIF
      rowVariableIdx=rowVariableIdx+1
    ENDDO
    IF(.NOT.rowVariableFound) NULLIFY(boundaryConditionsRowVariable)

    EXITS("BoundaryConditions_LHSVariableExists")
    RETURN
999 NULLIFY(boundaryConditionsRowVariable)
998 ERRORSEXITS("BoundaryConditions_LHSVariableExists",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_LHSVariableExists

  !
  !================================================================================================================================
  !

  !>Find the boundary conditions row variable for a given LHS variable
  SUBROUTINE BoundaryConditions_LHSVariableGet(boundaryConditions,lhsVariable,boundaryConditionsRowVariable,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to get the boundary conditions row variable for.
    TYPE(FieldVariableType), POINTER :: lhsVariable !<A pointer to the field variable to get the boundary conditions row variable for.
    TYPE(BoundaryConditionsRowVariableType), POINTER, INTENT(OUT) :: boundaryConditionsRowVariable !<On return, a pointer to the boundary conditions row variable. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("BoundaryConditions_LHSVariableGet",err,error,*998)

    IF(ASSOCIATED(boundaryConditionsRowVariable)) &
      & CALL FlagError("Boundary conditions row variable is already associated.",err,error,*998)
    
    CALL BoundaryConditions_LHSVariableExists(boundaryConditions,lhsVariable,boundaryConditionsRowVariable,err,error,*999)    

    IF(.NOT.ASSOCIATED(boundaryConditionsRowVariable)) THEN
      localError="Boundary condition row variable is not associated"
      IF(ASSOCIATED(lhsVariable)) THEN
        localError=localError//" for LHS variable type "//TRIM(NumberToVString(lhsVariable%variableType,"*",err,error))
        IF(ASSOCIATED(lhsVariable%field)) THEN
          localError=localError//" of field number "//TRIM(NumberToVString(lhsVariable%field%userNumber,"*",err,error))
          IF(ASSOCIATED(lhsVariable%field%region)) THEN
            localError=localError//" of region number "//TRIM(NumberToVString(lhsVariable%field%region%userNumber,"*",err,error))
          ELSE IF(ASSOCIATED(lhsVariable%field%interface)) THEN
            localError=localError//" of interface number "// &
              & TRIM(NumberToVString(lhsVariable%field%interface%userNumber,"*",err,error))
            IF(ASSOCIATED(lhsVariable%field%interface%parentRegion)) &
              & localError=localError//" of parent region number "// &
              & TRIM(NumberToVString(lhsVariable%field%interface%parentRegion%userNumber,"*",err,error))
          ENDIF
        ENDIF
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    EXITS("BoundaryConditions_LHSVariableGet")
    RETURN
999 NULLIFY(boundaryConditionsRowVariable)
998 ERRORSEXITS("BoundaryConditions_LHSVariableGet",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_LHSVariableGet

  !
  !================================================================================================================================
  !

  !>Sets/changes the sparsity type for the Neumann integration matrices
  SUBROUTINE BoundaryConditions_NeumannSparsityTypeSet(boundaryConditions,sparsityType,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: sparsityType !<The matrix sparsity type to be set \see SolverRoutines_EquationsSparsityTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions

    ENTERS("BoundaryConditions_NeumannSparsityTypeSet",err,error,*999)

    IF(.NOT.ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions are not associated.",err,error,*999)
    
    SELECT CASE(sparsityType)
    CASE(BOUNDARY_CONDITION_SPARSE_MATRICES)
      boundaryConditions%neumannMatrixSparsity=BOUNDARY_CONDITION_SPARSE_MATRICES
    CASE(BOUNDARY_CONDITION_FULL_MATRICES)
      boundaryConditions%neumannMatrixSparsity=BOUNDARY_CONDITION_FULL_MATRICES
    CASE DEFAULT
      CALL FlagError("The specified Neumann integration matrix sparsity type of "// &
        & TRIM(NumberToVString(sparsityType,"*",err,error))//" is invalid.",err,error,*999)
    END SELECT

    EXITS("BoundaryConditions_NeumannSparsityTypeSet")
    RETURN
999 ERRORSEXITS("BoundaryConditions_NeumannSparsityTypeSet",err,error)
    RETURN 1

  END SUBROUTINE BoundaryConditions_NeumannSparsityTypeSet

  !
  !================================================================================================================================
  !

  !>Find the boundary conditions row variable for a given field RHS variable if it exists
  SUBROUTINE BoundaryConditions_RHSVariableExists(boundaryConditions,rhsVariable,boundaryConditionsRowVariable,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to see if the boundary conditions RHS variable exists for
    TYPE(FieldVariableType), POINTER :: rhsVariable !<A pointer to the RHS variable to check the boundary conditions row variable for.
    TYPE(BoundaryConditionsRowVariableType), POINTER, INTENT(OUT) :: boundaryConditionsRowVariable !<On return, a pointer to the boundary conditions row variable, or NULL if it wasn't found. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: rhsInterfaceUserNumber,rhsVarType,rhsVariableFieldUserNumber,rhsVariableRegionUserNumber, &
      & matchInterfaceUserNumber,matchVariableType,matchVariableFieldUserNumber,matchVariableRegionUserNumber, &
      & numberOfBoundaryConditionsRowVariables,rowVariableIdx
    LOGICAL :: rhsVariableIsRegion,matchVariableIsInterface,matchVariableIsRegion,rowVariableFound
    TYPE(FieldType), POINTER :: rhsVariableField,matchVariableField
    TYPE(FieldVariableType), POINTER :: matchVariable
    TYPE(InterfaceType), POINTER :: rhsVariableInterface,matchVariableInterface
    TYPE(RegionType), POINTER :: rhsInterfaceParentRegion,rhsVariableRegion,matchInterfaceParentRegion,matchVariableRegion
    TYPE(VARYING_STRING) :: localError

    ENTERS("BoundaryConditions_RHSVariableExists",err,error,*998)

    IF(ASSOCIATED(boundaryConditionsRowVariable)) &
      & CALL FlagError("Boundary conditions row variable is already associated.",err,error,*998)

    CALL FieldVariable_VariableTypeGet(rhsVariable,rhsVarType,err,error,*999)
    NULLIFY(rhsVariableField)
    CALL FieldVariable_FieldGet(rhsVariable,rhsVariableField,err,error,*999)
    CALL Field_UserNumberGet(rhsVariableField,rhsVariableFieldUserNumber,err,error,*999)
    CALL Field_IsRegionField(rhsVariableField,rhsVariableIsRegion,err,error,*999)
    NULLIFY(rhsVariableRegion)
    NULLIFY(rhsVariableInterface)
    IF(rhsVariableIsRegion) THEN
      CALL Field_RegionGet(rhsVariableField,rhsVariableRegion,err,error,*999)
      CALL Region_UserNumberGet(rhsVariableRegion,rhsVariableRegionUserNumber,err,error,*999)
    ELSE
      CALL Field_InterfaceGet(rhsVariableField,rhsVariableInterface,err,error,*999)
      CALL Interface_UserNumberGet(rhsVariableInterface,rhsInterfaceUserNumber,err,error,*999)
      NULLIFY(rhsInterfaceParentRegion)
      CALL Interface_ParentRegionGet(rhsVariableInterface,rhsInterfaceParentRegion,err,error,*999)
      CALL Region_UserNumberGet(rhsInterfaceParentRegion,rhsVariableRegionUserNumber,err,error,*999)      
    ENDIF
        
    CALL BoundaryConditions_NumberOfRowVariablesGet(boundaryConditions,numberOfBoundaryConditionsRowVariables,err,error,*999)
    rowVariableFound=.FALSE.
    rowVariableIdx=1
    DO WHILE(rowVariableIdx<=numberOfBoundaryConditionsRowVariables.AND..NOT.rowVariableFound)
      NULLIFY(boundaryConditionsRowVariable)
      CALL BoundaryConditions_RowVariableIndexGet(boundaryConditions,rowVariableIdx,boundaryConditionsRowVariable,err,error,*999)
      NULLIFY(matchVariable)
      CALL BoundaryConditionsRowVariable_RHSVariableExists(boundaryConditionsRowVariable,matchVariable,err,error,*999)
      IF(ASSOCIATED(matchVariable)) THEN
        CALL FieldVariable_VariableTypeGet(matchVariable,matchVariableType,err,error,*999)
        NULLIFY(matchVariableField)
        CALL FieldVariable_FieldGet(matchVariable,matchVariableField,err,error,*999)
        CALL Field_UserNumberGet(matchVariableField,matchVariableFieldUserNumber,err,error,*999)
        CALL Field_IsRegionField(matchVariableField,matchVariableIsRegion,err,error,*999)
        IF(rhsVariableIsRegion) THEN
          IF(matchVariableIsRegion) THEN
            NULLIFY(matchVariableRegion)
            CALL Field_RegionGet(matchVariableField,matchVariableRegion,err,error,*999)
            CALL Region_UserNumberGet(matchVariableRegion,matchVariableRegionUserNumber,err,error,*999)
            IF(rhsVariableRegionUserNumber==matchVariableRegionUserNumber.AND. &
              & rhsVariableFieldUserNumber==matchVariableFieldUserNumber.AND. &
              & rhsVarType==matchVariableType) rowVariableFound=.TRUE.
          ENDIF
        ELSE
          CALL Field_IsInterfaceField(matchVariableField,matchVariableIsInterface,err,error,*999)
          IF(matchVariableIsInterface) THEN
            NULLIFY(matchVariableInterface)
            CALL Field_InterfaceGet(matchVariableField,matchVariableInterface,err,error,*999)
            CALL Interface_UserNumberGet(matchVariableInterface,matchInterfaceUserNumber,err,error,*999)
            NULLIFY(matchInterfaceParentRegion)
            CALL Interface_ParentRegionGet(matchVariableInterface,matchInterfaceParentRegion,err,error,*999)
            CALL Region_UserNumberGet(matchInterfaceParentRegion,matchVariableRegionUserNumber,err,error,*999)
            IF(rhsInterfaceUserNumber==matchInterfaceUserNumber.AND. &
              & rhsVariableRegionUserNumber==matchVariableRegionUserNumber.AND. &
              & rhsVariableFieldUserNumber==matchVariableFieldUserNumber.AND. &
              & rhsVarType==matchVariableType) rowVariableFound=.TRUE.
          ENDIF
        ENDIF
      ENDIF
      rowVariableIdx=rowVariableIdx+1
    ENDDO
    IF(.NOT.rowVariableFound) NULLIFY(boundaryConditionsRowVariable)

    EXITS("BoundaryConditions_RHSVariableExists")
    RETURN
999 NULLIFY(boundaryConditionsRowVariable)
998 ERRORSEXITS("BoundaryConditions_RHSVariableExists",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_RHSVariableExists

  !
  !================================================================================================================================
  !

  !>Find the boundary conditions row variable for a given RHS variable
  SUBROUTINE BoundaryConditions_RHSVariableGet(boundaryConditions,rhsVariable,boundaryConditionsRowVariable,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to get the boundary conditions row variable for.
    TYPE(FieldVariableType), POINTER :: rhsVariable !<A pointer to the field variable to get the boundary conditions row variable for.
    TYPE(BoundaryConditionsRowVariableType), POINTER, INTENT(OUT) :: boundaryConditionsRowVariable !<On return, a pointer to the boundary conditions row variable. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("BoundaryConditions_RHSVariableGet",err,error,*998)

    IF(ASSOCIATED(boundaryConditionsRowVariable)) &
      & CALL FlagError("Boundary conditions row variable is already associated.",err,error,*998)
    
    CALL BoundaryConditions_RHSVariableExists(boundaryConditions,rhsVariable,boundaryConditionsRowVariable,err,error,*999)    

    IF(.NOT.ASSOCIATED(boundaryConditionsRowVariable)) THEN
      localError="Boundary condition row variable is not associated"
      IF(ASSOCIATED(rhsVariable)) THEN
        localError=localError//" for RHS variable type "//TRIM(NumberToVString(rhsVariable%variableType,"*",err,error))
        IF(ASSOCIATED(rhsVariable%field)) THEN
          localError=localError//" of field number "//TRIM(NumberToVString(rhsVariable%field%userNumber,"*",err,error))
          IF(ASSOCIATED(rhsVariable%field%region)) THEN
            localError=localError//" of region number "//TRIM(NumberToVString(rhsVariable%field%region%userNumber,"*",err,error))
          ELSE IF(ASSOCIATED(rhsVariable%field%interface)) THEN
            localError=localError//" of interface number "// &
              & TRIM(NumberToVString(rhsVariable%field%interface%userNumber,"*",err,error))
            IF(ASSOCIATED(rhsVariable%field%interface%parentRegion)) &
              & localError=localError//" of parent region number "// &
              & TRIM(NumberToVString(rhsVariable%field%interface%parentRegion%userNumber,"*",err,error))
          ENDIF
        ENDIF
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    EXITS("BoundaryConditions_RHSVariableGet")
    RETURN
999 NULLIFY(boundaryConditionsRowVariable)
998 ERRORSEXITS("BoundaryConditions_RHSVariableGet",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_RHSVariableGet

  !
  !================================================================================================================================
  !

  !>Find the boundary conditions row variable for a given field row variable if it exists
  SUBROUTINE BoundaryConditions_RowVariableExists(boundaryConditions,fieldVariable,boundaryConditionsRowVariable,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to see if the boundary conditions row variable exists for
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to check the boundary conditions row variable for.
    TYPE(BoundaryConditionsRowVariableType), POINTER, INTENT(OUT) :: boundaryConditionsRowVariable !<On return, a pointer to the boundary conditions row variable, or NULL if it wasn't found. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: fieldInterfaceUserNumber,fieldVarType,fieldVariableFieldUserNumber,fieldVariableRegionUserNumber, &
      & matchInterfaceUserNumber,matchVariableType,matchVariableFieldUserNumber,matchVariableRegionUserNumber, &
      & numberOfBoundaryConditionsRowVariables,rowVariableIdx
    LOGICAL :: fieldVariableIsRegion,matchVariableIsInterface,matchVariableIsRegion,rowVariableFound
    TYPE(FieldType), POINTER :: fieldVariableField,matchVariableField
    TYPE(FieldVariableType), POINTER :: matchVariable
    TYPE(InterfaceType), POINTER :: fieldVariableInterface,matchVariableInterface
    TYPE(RegionType), POINTER :: fieldInterfaceParentRegion,fieldVariableRegion,matchInterfaceParentRegion,matchVariableRegion
    TYPE(VARYING_STRING) :: localError

    ENTERS("BoundaryConditions_RowVariableExists",err,error,*998)

    IF(ASSOCIATED(boundaryConditionsRowVariable)) &
      & CALL FlagError("Boundary conditions row variable is already associated.",err,error,*998)

    CALL FieldVariable_VariableTypeGet(fieldVariable,fieldVarType,err,error,*999)
    NULLIFY(fieldVariableField)
    CALL FieldVariable_FieldGet(fieldVariable,fieldVariableField,err,error,*999)
    CALL Field_UserNumberGet(fieldVariableField,fieldVariableFieldUserNumber,err,error,*999)
    CALL Field_IsRegionField(fieldVariableField,fieldVariableIsRegion,err,error,*999)
    NULLIFY(fieldVariableRegion)
    NULLIFY(fieldVariableInterface)
    IF(fieldVariableIsRegion) THEN
      CALL Field_RegionGet(fieldVariableField,fieldVariableRegion,err,error,*999)
      CALL Region_UserNumberGet(fieldVariableRegion,fieldVariableRegionUserNumber,err,error,*999)
    ELSE
      CALL Field_InterfaceGet(fieldVariableField,fieldVariableInterface,err,error,*999)
      CALL Interface_UserNumberGet(fieldVariableInterface,fieldInterfaceUserNumber,err,error,*999)
      NULLIFY(fieldInterfaceParentRegion)
      CALL Interface_ParentRegionGet(fieldVariableInterface,fieldInterfaceParentRegion,err,error,*999)
      CALL Region_UserNumberGet(fieldInterfaceParentRegion,fieldVariableRegionUserNumber,err,error,*999)      
    ENDIF
        
    CALL BoundaryConditions_NumberOfRowVariablesGet(boundaryConditions,numberOfBoundaryConditionsRowVariables,err,error,*999)
    rowVariableFound=.FALSE.
    rowVariableIdx=1
    DO WHILE(rowVariableIdx<=numberOfBoundaryConditionsRowVariables.AND..NOT.rowVariableFound)
      NULLIFY(boundaryConditionsRowVariable)
      CALL BoundaryConditions_RowVariableIndexGet(boundaryConditions,rowVariableIdx,boundaryConditionsRowVariable,err,error,*999)
      NULLIFY(matchVariable)
      CALL BoundaryConditionsRowVariable_LHSVariableGet(boundaryConditionsRowVariable,matchVariable,err,error,*999)
      CALL FieldVariable_VariableTypeGet(matchVariable,matchVariableType,err,error,*999)
      NULLIFY(matchVariableField)
      CALL FieldVariable_FieldGet(matchVariable,matchVariableField,err,error,*999)
      CALL Field_UserNumberGet(matchVariableField,matchVariableFieldUserNumber,err,error,*999)
      CALL Field_IsRegionField(matchVariableField,matchVariableIsRegion,err,error,*999)
      IF(fieldVariableIsRegion) THEN
        IF(matchVariableIsRegion) THEN
          NULLIFY(matchVariableRegion)
          CALL Field_RegionGet(matchVariableField,matchVariableRegion,err,error,*999)
          CALL Region_UserNumberGet(matchVariableRegion,matchVariableRegionUserNumber,err,error,*999)
          IF(fieldVariableRegionUserNumber==matchVariableRegionUserNumber.AND. &
            & fieldVariableFieldUserNumber==matchVariableFieldUserNumber.AND. &
            & fieldVarType==matchVariableType) rowVariableFound=.TRUE.
        ENDIF
      ELSE
        CALL Field_IsInterfaceField(matchVariableField,matchVariableIsInterface,err,error,*999)
        IF(matchVariableIsInterface) THEN
          NULLIFY(matchVariableInterface)
          CALL Field_InterfaceGet(matchVariableField,matchVariableInterface,err,error,*999)
          CALL Interface_UserNumberGet(matchVariableInterface,matchInterfaceUserNumber,err,error,*999)
          NULLIFY(matchInterfaceParentRegion)
          CALL Interface_ParentRegionGet(matchVariableInterface,matchInterfaceParentRegion,err,error,*999)
          CALL Region_UserNumberGet(matchInterfaceParentRegion,matchVariableRegionUserNumber,err,error,*999)
          IF(fieldInterfaceUserNumber==matchInterfaceUserNumber.AND. &
            & fieldVariableRegionUserNumber==matchVariableRegionUserNumber.AND. &
            & fieldVariableFieldUserNumber==matchVariableFieldUserNumber.AND. &
            & fieldVarType==matchVariableType) rowVariableFound=.TRUE.
        ENDIF
      ENDIF
      rowVariableIdx=rowVariableIdx+1
    ENDDO
    IF(.NOT.rowVariableFound) NULLIFY(boundaryConditionsRowVariable)

    EXITS("BoundaryConditions_RowVariableExists")
    RETURN
999 NULLIFY(boundaryConditionsRowVariable)
998 ERRORSEXITS("BoundaryConditions_RowVariableExists",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_RowVariableExists

  !
  !================================================================================================================================
  !

  !>Find the boundary conditions row variable for a given field variable
  SUBROUTINE BoundaryConditions_RowVariableGet(boundaryConditions,fieldVariable,boundaryConditionsRowVariable,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to get the boundary conditions row variable for.
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to get the boundary conditions row variable for.
    TYPE(BoundaryConditionsRowVariableType), POINTER, INTENT(OUT) :: boundaryConditionsRowVariable !<On return, a pointer to the boundary conditions row variable. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("BoundaryConditions_RowVariableGet",err,error,*998)

    IF(ASSOCIATED(boundaryConditionsRowVariable)) &
      & CALL FlagError("Boundary conditions row variable is already associated.",err,error,*998)
    
    CALL BoundaryConditions_RowVariableExists(boundaryConditions,fieldVariable,boundaryConditionsRowVariable,err,error,*999)    

    IF(.NOT.ASSOCIATED(boundaryConditionsRowVariable)) THEN
      localError="Boundary condition row variable is not associated"
      IF(ASSOCIATED(fieldVariable)) THEN
        localError=localError//" for variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
        IF(ASSOCIATED(fieldVariable%field)) THEN
          localError=localError//" of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
          IF(ASSOCIATED(fieldVariable%field%region)) THEN
            localError=localError//" of region number "//TRIM(NumberToVString(fieldVariable%field%region%userNumber,"*",err,error))
          ELSE IF(ASSOCIATED(fieldVariable%field%interface)) THEN
            localError=localError//" of interface number "// &
              & TRIM(NumberToVString(fieldVariable%field%interface%userNumber,"*",err,error))
            IF(ASSOCIATED(fieldVariable%field%interface%parentRegion)) &
              & localError=localError//" of parent region number "// &
              & TRIM(NumberToVString(fieldVariable%field%interface%parentRegion%userNumber,"*",err,error))
          ENDIF
        ENDIF
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    EXITS("BoundaryConditions_RowVariableGet")
    RETURN
999 NULLIFY(boundaryConditionsRowVariable)
998 ERRORSEXITS("BoundaryConditions_RowVariableGet",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_RowVariableGet

  !
  !================================================================================================================================
  !

  !>Adds a boundary conditions row variable if that variable has not already been added, otherwise do nothing.
  SUBROUTINE BoundaryConditions_RowVariableAdd(boundaryConditions,lhsVariable,rhsVariable,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to add row variable for.
    TYPE(FieldVariableType), POINTER :: lhsVariable !<A pointer to the row variable to add to the boundary conditions for.
    TYPE(FieldVariableType), POINTER :: rhsVariable !<A pointer to the row variable to add to the boundary conditions for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,numberOfDOFs,rowVariableIdx,totalNumberOfDOFs
    TYPE(BoundaryConditionsRowVariablePtrType), ALLOCATABLE :: newBoundaryConditionsRowVariables(:)
    TYPE(BoundaryConditionsRowVariableType), POINTER :: boundaryConditionsRowVariable
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("BoundaryConditions_RowVariableAdd",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions is not associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(lhsVariable)) CALL FlagError("Row variable is not associated.",err,error,*998)
#endif    
    
    !Check if boundary conditions variable has already been added, if so then we don't do anything as different equations
    !sets can have the same dependent field variables and will both want to add the variable
    NULLIFY(boundaryConditionsRowVariable)
    CALL BoundaryConditions_RowVariableExists(boundaryConditions,lhsVariable,boundaryConditionsRowVariable,err,error,*999)
    IF(.NOT.ASSOCIATED(boundaryConditionsRowVariable)) THEN
      !The row variable has not been added so add it to the list
      ALLOCATE(newBoundaryConditionsRowVariables(boundaryConditions%numberOfBoundaryConditionsRowVariables+1),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new boundary conditions row variables array.",err,error,*998)
      IF(ALLOCATED(boundaryConditions%boundaryConditionsRowVariables)) THEN
        DO rowVariableIdx=1,boundaryConditions%numberOfBoundaryConditionsRowVariables
          newBoundaryConditionsRowVariables(rowVariableIdx)%ptr=> &
            & boundaryConditions%boundaryConditionsRowVariables(RowVariableIdx)%ptr
        ENDDO !rowVariableIdx
      ENDIF

      NULLIFY(newBoundaryConditionsRowVariables(boundaryConditions%numberOfBoundaryConditionsRowVariables+1)%ptr)
      CALL BoundaryConditionsRowVariable_Initialise(newBoundaryConditionsRowVariables(boundaryConditions% &
        & numberOfBoundaryConditionsRowVariables+1)%ptr,err,error,*999)
      boundaryConditionsRowVariable=>newBoundaryConditionsRowVariables(boundaryConditions% &
        & numberOfBoundaryConditionsRowVariables+1)%ptr
      boundaryConditionsRowVariable%boundaryConditions=>boundaryConditions
      boundaryConditionsRowVariable%lhsVariable=>lhsVariable
      boundaryConditionsRowVariable%rhsVariable=>rhsVariable
      CALL FieldVariable_NumberOfDOFsGet(lhsVariable,boundaryConditionsRowVariable%numberOfDOFs,err,error,*999)
      CALL FieldVariable_TotalNumberOfDOFsGet(lhsVariable,boundaryConditionsRowVariable%totalNumberOfDOFs,err,error,*999)
      CALL FieldVariable_NumberOfGlobalDOFsGet(lhsVariable,boundaryConditionsRowVariable%numberOfGLobalDOFs,err,error,*999)
      !TEMP WHILST WE SWITCH FROM GLOBAL DOFS
      !ALLOCATE(boundaryConditionsRowVariable%rowConditionTypes(boundaryConditionsRowVariable%totalNumberOfDOFs),STAT=err)
      ALLOCATE(boundaryConditionsRowVariable%rowConditionTypes(boundaryConditionsRowVariable%numberOfGlobalDOFs),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate row condition types.",err,error,*999)
      boundaryConditionsRowVariable%rowConditionTypes=BOUNDARY_CONDITION_FREE_ROW
      
      CALL MOVE_ALLOC(newBoundaryConditionsRowVariables,boundaryConditions%boundaryConditionsRowVariables)
      boundaryConditions%numberOfBoundaryConditionsRowVariables=boundaryConditions%numberOfBoundaryConditionsRowVariables+1

    ENDIF

    EXITS("BoundaryConditions_RowVariableAdd")
    RETURN
999 CALL BoundaryConditionsRowVariable_Finalise(boundaryConditionsRowVariable,dummyErr,dummyError,*998)
    IF(ALLOCATED(newBoundaryConditionsRowVariables)) DEALLOCATE(newBoundaryConditionsRowVariables)
998 ERRORSEXITS("BoundaryConditions_RowVariableAdd",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_RowVariableAdd

  !
  !================================================================================================================================
  !

  !>Sets the analytic value as a boundary condition on the specified node.
  SUBROUTINE BoundaryConditions_SetAnalyticBoundaryNode(boundaryConditions,numberOfDimensions,dependentVariable,componentNumber, &
    & domainNodes,localNodenumber,boundaryNode,tangents,normal,analyticValue,gradientAnalyticValue,hessianAnalyticValue, &
    & setVelocity,analyticVelocity,setAcceleration,analyticAcceleration,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to set the boundary condition for
    INTEGER(INTG), INTENT(IN) :: numberOfDimensions !<The number of dimensions for the problem.
    TYPE(FieldVariableType), POINTER :: dependentVariable !<The dependent field variable to set the boundary condition on.    
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The component number to set the boundary condition at
    TYPE(DomainNodesType), POINTER :: domainNodes !<A pointer to the domain nodes for the local node to set the boundary condition for
    INTEGER(INTG), INTENT(IN) :: localNodeNumber !<The local node number to set the boundary condition for
    LOGICAL, INTENT(IN) :: boundaryNode !<The boundary node flag for the node
    REAL(DP), INTENT(IN) :: tangents(:,:) !<tangents(componentIdx,tangentIdx). The tangent vectors for the node.
    REAL(DP), INTENT(IN) :: normal(:) !<normal(componentIdx). The normal vector at the node.
    REAL(DP), INTENT(IN) :: analyticValue !<The value of the analytic value at the node
    REAL(DP), INTENT(IN) :: gradientAnalyticValue(:) !<gradientAnalyticValue(componentIdx). The gradient of the analytic value at the node.
    REAL(DP), INTENT(IN) :: hessianAnalyticValue(:,:) !<hessianAnalyticValue(componentIdx1,componentIdx2). The Hessian of the analyticc value at the node.
    LOGICAL, INTENT(IN) :: setVelocity !<Is .TRUE. if velocity values are to be set, .FALSE. if not.
    REAL(DP), INTENT(IN) :: analyticVelocity !<The analytic velocity to set
    LOGICAL, INTENT(IN) :: setAcceleration !<Is .TRUE. if acceleration values are to be set, .FALSE. if not.
    REAL(DP), INTENT(IN) :: analyticAcceleration !<The analytic acceleration to set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx1,componentIdx2,derivativeGlobalIndex,derivativeIdx,localDOFIdx,numberOfNodeDerivatives, &
      & variableType
    REAL(DP) :: dofValue
    TYPE(VARYING_STRING) :: localError

    ENTERS("BoundaryConditions_SetAnalyticBoundaryNode",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions is not associated.",err,error,*999)
    CALL FieldVariable_AssertComponentNumberOK(dependentVariable,componentNumber,err,error,*999)
    CALL DomainNodes_LocalNumberCheck(domainNodes,localNodeNumber,err,error,*999)
#endif    
    CALL FieldVariable_VariableTypeGet(dependentVariable,variableType,err,error,*999)

    CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,localNodeNumber,numberOfNodeDerivatives,err,error,*999)
    !Loop over the derivatives
    DO derivativeIdx=1,numberOfNodeDerivatives
      CALL DomainNodes_DerivativeGlobalIndexGet(domainNodes,derivativeIdx,localNodeNumber,derivativeGlobalIndex,err,error,*999)
      !Default to version 1 of each node derivative
      CALL FieldVariable_LocalNodeDOFGet(dependentVariable,1,derivativeIdx,localNodeNumber,componentNumber,localDOFIdx, &
        & err,error,*999)
      dofValue=0.0_DP
      SELECT CASE(variableType)
      CASE(FIELD_U_VARIABLE_TYPE) 
        SELECT CASE(derivativeGlobalIndex)
        CASE(NO_GLOBAL_DERIV)
          dofValue=analyticValue
        CASE(GLOBAL_DERIV_S1)
          !Compute grad(u).s1
          DO componentIdx1=1,numberOfDimensions
            dofValue=dofValue+gradientAnalyticValue(componentIdx1)*tangents(componentIdx1,1)
          ENDDO !componentIdx1
        CASE(GLOBAL_DERIV_S2)
          !Compute grad(u).s2
          DO componentIdx1=1,numberOfDimensions
            dofValue=dofValue+gradientAnalyticValue(componentIdx1)*tangents(componentIdx1,2)
          ENDDO !componentIdx1
        CASE(GLOBAL_DERIV_S3)
          !Compute grad(u).s3
          DO componentIdx1=1,numberOfDimensions
            dofValue=dofValue+gradientAnalyticValue(componentIdx1)*tangents(componentIdx1,3)
          ENDDO !componentIdx1
        CASE(GLOBAL_DERIV_S1_S2)
          !Compute grad s1^t.H(u).s2
          DO componentIdx2=1,numberOfDimensions
            DO componentIdx1=1,numberOfDimensions
              dofValue=dofValue+ &
                & tangents(componentIdx1,1)*hessianAnalyticValue(componentIdx1,componentIdx2)*tangents(componentIdx2,2)
            ENDDO !componentIdx1
          ENDDO !componentIdx2
        CASE(GLOBAL_DERIV_S1_S3)
          !Compute grad s1^t.H(u).s3
          DO componentIdx2=1,numberOfDimensions
            DO componentIdx1=1,numberOfDimensions
              dofValue=dofValue+ &
                & tangents(componentIdx1,1)*hessianAnalyticValue(componentIdx1,componentIdx2)*tangents(componentIdx2,3)
            ENDDO !componentIdx1
          ENDDO !componentIdx2
        CASE(GLOBAL_DERIV_S2_S3)
          !Compute grad s2^t.H(u).s3
          DO componentIdx2=1,numberOfDimensions
            DO componentIdx1=1,numberOfDimensions
              dofValue=dofValue+ &
                & tangents(componentIdx1,2)*hessianAnalyticValue(componentIdx1,componentIdx2)*tangents(componentIdx2,3)
            ENDDO !componentIdx1
          ENDDO !componentIdx2
        CASE(GLOBAL_DERIV_S1_S2_S3)
!!TODO: not implemented
          dofValue=0.0_DP
        CASE DEFAULT
          localError="The global derivative index of "//TRIM(NumberToVString(derivativeGlobalIndex,"*",err,error))// &
            & " for derivative index "//TRIM(NumberToVString(derivativeIdx,"*",err,error))//" of local node number "// &
            & TRIM(NumberToVString(localNodeNumber,"*",err,error))//" of component number "// &
            & TRIM(NumberToVString(componentNumber,"*",err,error))//" of variable type "// &
            & TRIM(NumberToVString(variableType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(FIELD_DELUDELN_VARIABLE_TYPE)
        SELECT CASE(derivativeGlobalIndex)
        CASE(NO_GLOBAL_DERIV)
          !Compute grad(u).n
          DO componentIdx1=1,numberOfDimensions
            dofValue=dofValue+gradientAnalyticValue(componentIdx1)*normal(componentIdx1)
          ENDDO !componentIdx1
          dofValue=analyticValue
        CASE(GLOBAL_DERIV_S1)
          !Compute grad s1^t.H(u).n
          DO componentIdx2=1,numberOfDimensions
            DO componentIdx1=1,numberOfDimensions
              dofValue=dofValue+ &
                & tangents(componentIdx1,1)*hessianAnalyticValue(componentIdx1,componentIdx2)*normal(componentIdx2)
            ENDDO !componentIdx1
          ENDDO !componentIdx2
        CASE(GLOBAL_DERIV_S2)
          !Compute grad s2^t.H(u).n
          DO componentIdx2=1,numberOfDimensions
            DO componentIdx1=1,numberOfDimensions
              dofValue=dofValue+ &
                & tangents(componentIdx1,2)*hessianAnalyticValue(componentIdx1,componentIdx2)*normal(componentIdx2)
            ENDDO !componentIdx1
          ENDDO !componentIdx2
        CASE(GLOBAL_DERIV_S3)
          !Compute grad s3^t.H(u).n
          DO componentIdx2=1,numberOfDimensions
            DO componentIdx1=1,numberOfDimensions
              dofValue=dofValue+ &
                & tangents(componentIdx1,3)*hessianAnalyticValue(componentIdx1,componentIdx2)*normal(componentIdx2)
            ENDDO !componentIdx1
          ENDDO !componentIdx2
        CASE(GLOBAL_DERIV_S1_S2)
!!TODO: not implemented
          dofValue=0.0_DP
        CASE(GLOBAL_DERIV_S1_S3)
!!TODO: not implemented
          dofValue=0.0_DP
        CASE(GLOBAL_DERIV_S2_S3)
!!TODO: not implemented
          dofValue=0.0_DP
        CASE(GLOBAL_DERIV_S1_S2_S3)
!!TODO: not implemented
          dofValue=0.0_DP
        CASE DEFAULT
          localError="The global derivative index of "//TRIM(NumberToVString(derivativeGlobalIndex,"*",err,error))// &
            & " for derivative index "//TRIM(NumberToVString(derivativeIdx,"*",err,error))//" of local node number "// &
            & TRIM(NumberToVString(localNodeNumber,"*",err,error))//" of component number "// &
            & TRIM(NumberToVString(componentNumber,"*",err,error))//" of variable type "// &
            & TRIM(NumberToVString(variableType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        !Do nothing or error???
      END SELECT
      IF(boundaryNode) THEN
        SELECT CASE(variableType)
        CASE(FIELD_U_VARIABLE_TYPE) !Dirichlet
          CALL BoundaryConditions_SetLocalDOF(boundaryConditions,dependentVariable,localDOFIdx,BOUNDARY_CONDITION_FIXED, &
            & dofValue,err,error,*999)
          IF(setVelocity) THEN
            CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_VELOCITY_VALUES_SET_TYPE,localDOFIdx, &
              & analyticVelocity,err,error,*999)
            CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_ANALYTIC_VELOCITY_VALUES_SET_TYPE,localDOFIdx, &
              & analyticVelocity,err,error,*999)
            IF(setAcceleration) THEN
              CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_ACCELERATION_VALUES_SET_TYPE, &
                & localDOFIdx,analyticAcceleration,err,error,*999)
              CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_ANALYTIC_ACCELERATION_VALUES_SET_TYPE, &
                & localDOFIdx,analyticAcceleration,err,error,*999)
            ENDIF
          ENDIF
        CASE(FIELD_DELUDELN_VARIABLE_TYPE) !Neummann
          CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,localDOFIdx, &
            & dofValue,err,error,*999)
        CASE DEFAULT
          !Do nothing
        END SELECT
      ELSE
        CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,localDOFIdx, &
          & dofValue,err,error,*999)
        IF(setVelocity) THEN
          CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_VELOCITY_VALUES_SET_TYPE,localDOFIdx, &
            & analyticVelocity,err,error,*999)
          CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_ANALYTIC_VELOCITY_VALUES_SET_TYPE,localDOFIdx, &
            & analyticVelocity,err,error,*999)
          IF(setAcceleration) THEN
            CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_ACCELERATION_VALUES_SET_TYPE, &
              & localDOFIdx,analyticAcceleration,err,error,*999)
            CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_ANALYTIC_ACCELERATION_VALUES_SET_TYPE, &
              & localDOFIdx,analyticAcceleration,err,error,*999)
          ENDIF
        ENDIF
      ENDIF
    ENDDO !derivativeIdx
    
    EXITS("BoundaryConditions_SetAnalyticBoundaryNode")
    RETURN
999 ERRORSEXITS("BoundaryConditions_SetAnalyticBoundaryNode",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_SetAnalyticBoundaryNode

  !
  !================================================================================================================================
  !

  !>Sets a boundary condition on the specified constant. \see OpenCMISS::Iron::cmfe_BoundaryConditions_SetConstant
  SUBROUTINE BoundaryConditions_SetConstantField(boundaryConditions,field,variableType,componentNumber,condition,bcValue, &
    & err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FieldType), POINTER :: field !<The dependent field to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: variableType !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The component number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: condition !<The boundary condition type to set \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
    REAL(DP), INTENT(IN) :: bcValue !<The value of the boundary condition to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(FieldVariableType), POINTER :: fieldVariable

    ENTERS("BoundaryConditions_SetConstantField",err,error,*999)

    NULLIFY(fieldVariable)
    CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
    CALL BoundaryConditions_SetConstantVariable(boundaryConditions,fieldVariable,componentNumber,condition,bcValue,err,error,*999)
    
    EXITS("BoundaryConditions_SetConstantField")
    RETURN
999 ERRORSEXITS("BoundaryConditions_SetConstantField",err,error)
    RETURN 1

  END SUBROUTINE BoundaryConditions_SetConstantField

  !
  !================================================================================================================================
  !

  !>Sets a boundary condition on the specified constant. 
  SUBROUTINE BoundaryConditions_SetConstantVariable(boundaryConditions,fieldVariable,componentNumber,condition,bcValue,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FieldVariableType), POINTER :: fieldVariable !<The dependent field variable to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The component number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: condition !<The boundary condition type to set \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
    REAL(DP), INTENT(IN) :: bcValue !<The value of the boundary condition to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: localDOF,globalDOF

    ENTERS("BoundaryConditions_SetConstantVariable",err,error,*999)

    !Note: This routine is for constant interpolation
    CALL BoundaryConditions_AssertNotFinished(boundaryConditions,err,error,*999)
    CALL FieldVariable_ConstantDOFGet(fieldVariable,componentNumber,localDOF,err,error,*999)
    CALL BoundaryConditions_CheckInterpolationType(condition,fieldVariable,componentNumber,err,error,*999)
    CALL BoundaryConditions_SetLocalDOF(boundaryConditions,fieldVariable,localDOF,condition,bcValue,err,error,*999)
    
    EXITS("BoundaryConditions_SetConstantVariable")
    RETURN
999 ERRORSEXITS("BoundaryConditions_SetConstantVariable",err,error)
    RETURN 1

  END SUBROUTINE BoundaryConditions_SetConstantVariable

  !
  !================================================================================================================================
  !

  !>Sets a boundary condition on the specified user element. \see OpenCMISS::Iron::cmfe_BoundaryConditions_SetElement
  SUBROUTINE BoundaryConditions_SetElementField(boundaryConditions,field,variableType,userElementNumber,componentNumber, &
    & condition,bcValue,err,error,*)
    
    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FieldType), POINTER :: field !<The dependent field to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: variableType !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The component number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: condition !<The boundary condition type to set \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
    REAL(DP), INTENT(IN) :: bcValue !<The value of the boundary condition to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(FieldVariableType), POINTER :: fieldVariable
 
    ENTERS("BoundaryConditions_SetElementField",err,error,*999)

    !Note: this routine is for element based interpolation
    NULLIFY(fieldVariable)
    CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
    CALL BoundaryConditions_SetElementVariable(boundaryConditions,fieldVariable,userElementNumber,componentNumber, &
      & condition,bcValue,err,error,*999)
 
    EXITS("BoundaryConditions_SetElementField")
    RETURN
999 ERRORSEXITS("BoundaryConditions_SetElementField",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_SetElementField

  !
  !================================================================================================================================
  !

  !>Sets a boundary condition on the specified user element. 
  SUBROUTINE BoundaryConditions_SetElementVariable(boundaryConditions,fieldVariable,userElementNumber,componentNumber, &
    & condition,bcValue,err,error,*)
    
    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FieldVariableType), POINTER :: fieldVariable !<The dependent field variable to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The component number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: condition !<The boundary condition type to set \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
    REAL(DP), INTENT(IN) :: bcValue !<The value of the boundary condition to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: localDOF
    LOGICAL :: ghostDOF

    ENTERS("BoundaryConditions_SetElementVariable",err,error,*999)

    !Note: this routine is for element based interpolation
    CALL BoundaryConditions_AssertNotFinished(boundaryConditions,err,error,*999)
    CALL FieldVariable_UserElementDOFGet(fieldVariable,userElementNumber,componentNumber,localDOF,ghostDOF,err,error,*999)
    CALL BoundaryConditions_CheckInterpolationType(condition,fieldVariable,componentNumber,err,error,*999)
    CALL BoundaryConditions_SetLocalDOF(boundaryConditions,fieldVariable,localDOF,condition,bcValue,err,error,*999)
 
    EXITS("BoundaryConditions_SetElementVariable")
    RETURN
999 ERRORSEXITS("BoundaryConditions_SetElementVariable",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_SetElementVariable

  !
  !================================================================================================================================
  !

  !>Sets a boundary condition on the specified DOF.
  SUBROUTINE BoundaryConditions_SetLocalDOF0(boundaryConditions,fieldVariable,dofIndex,condition,bcValue,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FieldVariableType), POINTER :: fieldVariable !<The dependent field variable to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: dofIndex !<The local dof index to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: condition !<The boundary condition type to set \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
    REAL(DP), INTENT(IN) :: bcValue !<The value of the boundary condition to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("BoundaryConditions_SetLocalDOF0",err,error,*999)

    CALL BoundaryConditions_SetLocalDOF1(boundaryConditions,fieldVariable,[dofIndex],[condition],[bcValue],err,error,*999)

    EXITS("BoundaryConditions_SetLocalDOF0")
    RETURN
999 ERRORSEXITS("BoundaryConditions_SetLocalDOF0",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_SetLocalDOF0

  !
  !================================================================================================================================
  !

  !>Sets a boundary condition on the specified DOFs.
  SUBROUTINE BoundaryConditions_SetLocalDOF1(boundaryConditions,fieldVariable,dofIndices,conditions,bcValues,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FieldVariableType), POINTER :: fieldVariable !<The dependent field to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: dofIndices(:) !<dofIndices(:). The local dof index for the i'th dof to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: conditions(:) !<conditions(:). The boundary condition type to set for the i'th dof \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
    REAL(DP), INTENT(IN) :: bcValues(:) !<bcValues(:). The value of the boundary condition for the i'th dof to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dofIdx,globalDOF,localDOF,numberOfLocalDOFs
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable
    TYPE(BoundaryConditionsRowVariableType), POINTER :: boundaryConditionsRowVariable
    TYPE(DomainMappingType), POINTER :: domainMapping
    TYPE(VARYING_STRING) :: localError

    ENTERS("BoundaryConditions_SetLocalDOF1",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(SIZE(dofIndices,1)/=SIZE(conditions,1)) THEN
      localError="The size of the DOF indices array of "//TRIM(NumberToVString(SIZE(dofIndices,1),"*",err,error))// &
        & " does not match the size of the fixed conditions array of "//TRIM(NumberToVString(SIZE(conditions,1),"*",err,error))// &
        & "."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(dofIndices,1)/=SIZE(bcValues,1)) THEN
      localError="The size of the DOF indices array of "//TRIM(NumberToVString(SIZE(dofIndices,1),"*",err,error))// &
        & " does not match the size of the values array of "//TRIM(NumberToVString(SIZE(bcValues,1),"*",err,error))//")."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    CALL BoundaryConditions_AssertNotFinished(boundaryConditions,err,error,*999)
    NULLIFY(boundaryConditionsVariable)
    CALL BoundaryConditions_VariableGet(boundaryConditions,fieldVariable,boundaryConditionsVariable,err,error,*999)
    NULLIFY(boundaryConditionsRowVariable)
    CALL BoundaryConditions_RowVariableExists(boundaryConditions,fieldVariable,boundaryConditionsRowVariable,err,error,*999)
    
    NULLIFY(domainMapping)
    CALL FieldVariable_DomainMappingGet(fieldVariable,domainMapping,err,error,*999)
    CALL DomainMapping_NumberOfLocalGet(domainMapping,numberOfLocalDOFs,err,error,*999)    
    DO dofIdx=1,SIZE(dofIndices,1)
      localDOF=dofIndices(dofIdx)
      IF(localDOF<1.OR.localDOF>numberOfLocalDOFs) THEN
        localError="The local DOF of  "//TRIM(NumberToVString(localDOF,"*",err,error))//" at DOF index "// &
          & TRIM(NumberToVString(dofIdx,"*",err,error))//" is invalid. The DOF should be >= 1 and <= "// &
          & TRIM(NumberToVString(numberOfLocalDOFs,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      CALL DomainMapping_LocalToGlobalGet(domainMapping,localDOF,globalDOF,err,error,*999)
      ! Set boundary condition and dof type
      CALL BoundaryConditionsVariable_ConditionTypeSet(boundaryConditionsVariable,globalDOF,conditions(dofIdx),err,error,*999)
      ! Update field sets with boundary condition value                  
      SELECT CASE(conditions(dofIdx))
      CASE(BOUNDARY_CONDITION_NONE)
        ! No field update
      CASE(BOUNDARY_CONDITION_FIXED)
        CALL FieldVariable_ParameterSetUpdateLocalDOF(fieldVariable,FIELD_VALUES_SET_TYPE,localDOF,bcValues(dofIdx),err,error,*999)
      CASE(BOUNDARY_CONDITION_FIXED_INLET)
        CALL FieldVariable_ParameterSetUpdateLocalDOF(fieldVariable,FIELD_VALUES_SET_TYPE,localDOF,bcValues(dofIdx),err,error,*999)
      CASE(BOUNDARY_CONDITION_FIXED_OUTLET)
        CALL FieldVariable_ParameterSetUpdateLocalDOF(fieldVariable,FIELD_VALUES_SET_TYPE,localDOF,bcValues(dofIdx),err,error,*999)
      CASE(BOUNDARY_CONDITION_FIXED_WALL)
        CALL FieldVariable_ParameterSetUpdateLocalDOF(fieldVariable,FIELD_VALUES_SET_TYPE,localDOF,bcValues(dofIdx),err,error,*999)
      CASE(BOUNDARY_CONDITION_MOVED_WALL)
        CALL FieldVariable_ParameterSetUpdateLocalDOF(fieldVariable,FIELD_VALUES_SET_TYPE,localDOF,bcValues(dofIdx),err,error,*999)
      CASE(BOUNDARY_CONDITION_FREE_WALL)
        CALL FieldVariable_ParameterSetUpdateLocalDOF(fieldVariable,FIELD_VALUES_SET_TYPE,localDOF,bcValues(dofIdx),err,error,*999)
      CASE(BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED)
        CALL FieldVariable_ParameterSetUpdateLocalDOF(fieldVariable,FIELD_VALUES_SET_TYPE,localDOF,bcValues(dofIdx), &
          & err,error,*999)
        CALL FieldVariable_ParameterSetUpdateLocalDOF(fieldVariable,FIELD_BOUNDARY_CONDITIONS_SET_TYPE,localDOF,bcValues(dofIdx), &
          & err,error,*999)
      CASE(BOUNDARY_CONDITION_FIXED_INCREMENTED) !For load increment loops
        CALL FieldVariable_ParameterSetUpdateLocalDOF(fieldVariable,FIELD_BOUNDARY_CONDITIONS_SET_TYPE,localDOF,bcValues(dofIdx), &
          & err,error,*999)
      CASE(BOUNDARY_CONDITION_PRESSURE)
        CALL FieldVariable_ParameterSetUpdateLocalDOF(fieldVariable,FIELD_PRESSURE_VALUES_SET_TYPE,localDOF,bcValues(dofIdx), &
          & err,error,*999)
      CASE(BOUNDARY_CONDITION_FIXED_PRESSURE)
        CALL FieldVariable_ParameterSetUpdateLocalDOF(fieldVariable,FIELD_VALUES_SET_TYPE,localDOF,bcValues(dofIdx),err,error,*999)
      CASE(BOUNDARY_CONDITION_PRESSURE_INCREMENTED)
        CALL FieldVariable_ParameterSetUpdateLocalDOF(fieldVariable,FIELD_PRESSURE_VALUES_SET_TYPE,localDOF,bcValues(dofIdx), &
          & err,error,*999)
      CASE(BOUNDARY_CONDITION_CORRECTION_MASS_INCREASE)
        ! No field update
      CASE(BOUNDARY_CONDITION_NEUMANN_POINT)
        CALL FieldVariable_ParameterSetUpdateLocalDOF(fieldVariable,FIELD_BOUNDARY_CONDITIONS_SET_TYPE,localDOF,bcValues(dofIdx), &
          & err,error,*999)
      CASE(BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)
        CALL FieldVariable_ParameterSetUpdateLocalDOF(fieldVariable,FIELD_BOUNDARY_CONDITIONS_SET_TYPE,localDOF,bcValues(dofIdx), &
          & err,error,*999)
      CASE(BOUNDARY_CONDITION_NEUMANN_INTEGRATED,BOUNDARY_CONDITION_NEUMANN_INTEGRATED_ONLY)
        CALL FieldVariable_ParameterSetUpdateLocalDOF(fieldVariable,FIELD_VALUES_SET_TYPE,localDOF,bcValues(dofIdx),err,error,*999)
      CASE(BOUNDARY_CONDITION_IMPERMEABLE_WALL)
        CALL FieldVariable_ParameterSetUpdateLocalDOF(fieldVariable,FIELD_IMPERMEABLE_FLAG_VALUES_SET_TYPE,localDOF, &
          & bcValues(dofIdx),err,error,*999)
      CASE(BOUNDARY_CONDITION_FIXED_FITTED,BOUNDARY_CONDITION_FIXED_NONREFLECTING, &
        & BOUNDARY_CONDITION_FIXED_CELLML,BOUNDARY_CONDITION_FIXED_STREE)
        CALL FieldVariable_ParameterSetUpdateLocalDOF(fieldVariable,FIELD_VALUES_SET_TYPE,localDOF,bcValues(dofIdx),err,error,*999)
      CASE DEFAULT
        localError="The specified boundary condition type for DOF index "//TRIM(NumberToVString(dofIdx,"*",err,error))//" of "// &
          & TRIM(NumberToVString(conditions(dofIdx),"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDDO !dofIdx

    EXITS("BoundaryConditions_SetLocalDOF1")
    RETURN
999 ERRORSEXITS("BoundaryConditions_SetLocalDOF1",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_SetLocalDOF1

  !
  !================================================================================================================================
  !

  !>Sets a boundary condition on the specified user node. \see OpenCMISS::Iron:cmfe_BoundaryConditions_SetNode
  SUBROUTINE BoundaryConditions_SetNodeField(boundaryConditions,field,variableType,versionNumber,derivativeNumber, &
    & userNodeNumber,componentNumber,condition,bcValue,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FieldType), POINTER :: field !<The dependent field to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: variableType !<The variable type to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: versionNumber !<The derivative version to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: derivativeNumber !<The derivative to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: userNodeNumber !<The user node number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The component number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: condition !<The boundary condition type to set \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
    REAL(DP), INTENT(IN) :: bcValue !<The value of the boundary condition to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: localDOF,globalDOF
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("BoundaryConditions_SetNodeField",err,error,*999)

    NULLIFY(fieldVariable)
    CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
    CALL BoundaryConditions_SetNodeVariable(boundaryConditions,fieldVariable,versionNumber,derivativeNumber, &
      & userNodeNumber,componentNumber,condition,bcValue,err,error,*999)
    
    EXITS("BoundaryConditions_SetNodeField")
    RETURN
999 ERRORSEXITS("BoundaryConditions_SetNodeField",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_SetNodeField

  !
  !================================================================================================================================
  !

  !>Sets a boundary condition on the specified user node. 
  SUBROUTINE BoundaryConditions_SetNodeVariable(boundaryConditions,fieldVariable,versionNumber,derivativeNumber, &
    & userNodeNumber,componentNumber,condition,bcValue,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to set the boundary condition for
    TYPE(FieldVariableType), POINTER :: fieldVariable !<The dependent field variable to set the boundary condition on.
    INTEGER(INTG), INTENT(IN) :: versionNumber !<The derivative version to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: derivativeNumber !<The derivative to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: userNodeNumber !<The user node number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The component number to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: condition !<The boundary condition type to set \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
    REAL(DP), INTENT(IN) :: bcValue !<The value of the boundary condition to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: localDOF,globalDOF
    TYPE(VARYING_STRING) :: localError

    ENTERS("BoundaryConditions_SetNodeVariable",err,error,*999)

    CALL BoundaryConditions_AssertNotFinished(boundaryConditions,err,error,*999)
    CALL FieldVariable_ComponentDOFGetUserNode(fieldVariable,versionNumber,derivativeNumber,userNodeNumber,componentNumber, &
      & localDOF,globalDOF,err,error,*999)
    CALL BoundaryConditions_CheckInterpolationType(condition,fieldVariable,componentNumber,err,error,*999)
    CALL BoundaryConditions_SetLocalDOF(boundaryConditions,fieldVariable,localDOF,condition,bcValue,err,error,*999)
    
    EXITS("BoundaryConditions_SetNodeVariable")
    RETURN
999 ERRORSEXITS("BoundaryConditions_SetNodeVariable",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_SetNodeVariable

  !
  !================================================================================================================================
  !

  !>Updates the analytic value as a boundary condition on the specified node.
  SUBROUTINE BoundaryConditions_UpdateAnalyticNode(boundaryConditions,numberOfDimensions,dependentVariable, &
    & componentNumber,domainNodes,localNodenumber,tangents,normal,analyticValue,gradientAnalyticValue, &
    & hessianAnalyticValue,setVelocity,analyticVelocity,setAcceleration,analyticAcceleration,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to update the boundary condition for
    INTEGER(INTG), INTENT(IN) :: numberOfDimensions !<The number of dimensions for the problem.
    TYPE(FieldVariableType), POINTER :: dependentVariable !<The dependent field variable to update the boundary condition on.    
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The component number to update the boundary condition at
    TYPE(DomainNodesType), POINTER :: domainNodes !<A pointer to the domain nodes for the local node to update the boundary condition for
    INTEGER(INTG), INTENT(IN) :: localNodeNumber !<The local node number to update the boundary condition for
    REAL(DP), INTENT(IN) :: tangents(:,:) !<tangents(componentIdx,tangentIdx). The tangent vectors for the node.
    REAL(DP), INTENT(IN) :: normal(:) !<normal(componentIdx). The normal vector at the node.
    REAL(DP), INTENT(IN) :: analyticValue !<The value of the analytic value at the node
    REAL(DP), INTENT(IN) :: gradientAnalyticValue(:) !<gradientAnalyticValue(componentIdx). The gradient of the analytic value at the node.
    REAL(DP), INTENT(IN) :: hessianAnalyticValue(:,:) !<hessianAnalyticValue(componentIdx1,componentIdx2). The Hessian of the analyticc value at the node.
    LOGICAL, INTENT(IN) :: setVelocity !<Is .TRUE. if velocity values are to be set, .FALSE. if not.
    REAL(DP), INTENT(IN) :: analyticVelocity !<The analytic velocity to set
    LOGICAL, INTENT(IN) :: setAcceleration !<Is .TRUE. if acceleration values are to be set, .FALSE. if not.
    REAL(DP), INTENT(IN) :: analyticAcceleration !<The analytic acceleration to set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx1,componentIdx2,conditionType,derivativeGlobalIndex,derivativeIdx,localDOFIdx, &
      & numberOfNodeDerivatives,variableType
    REAL(DP) :: dofValue
    TYPE(VARYING_STRING) :: localError
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable

    ENTERS("BoundaryConditions_UpdateAnalyticNode",err,error,*999)

#ifdef WITH_PRECHECKS    
    CALL FieldVariable_AssertComponentNumberOK(dependentVariable,componentNumber,err,error,*999)
    CALL DomainNodes_LocalNumberCheck(domainNodes,localNodeNumber,err,error,*999)
#endif    
    NULLIFY(boundaryConditionsVariable)
    CALL BoundaryConditions_VariableGet(boundaryConditions,dependentVariable,boundaryConditionsVariable,err,error,*999)
    CALL FieldVariable_VariableTypeGet(dependentVariable,variableType,err,error,*999)
 
    CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,localNodeNumber,numberOfNodeDerivatives,err,error,*999)
    !Loop over the derivatives
    DO derivativeIdx=1,numberOfNodeDerivatives
      CALL DomainNodes_DerivativeGlobalIndexGet(domainNodes,derivativeIdx,localNodeNumber,derivativeGlobalIndex,err,error,*999)
      !Default to version 1 of each node derivative
      CALL FieldVariable_LocalNodeDOFGet(dependentVariable,1,derivativeIdx,localNodeNumber,componentNumber,localDOFIdx, &
        & err,error,*999)
      dofValue=0.0_DP
      SELECT CASE(variableType)
      CASE(FIELD_U_VARIABLE_TYPE) 
        SELECT CASE(derivativeGlobalIndex)
        CASE(NO_GLOBAL_DERIV)
          dofValue=analyticValue
        CASE(GLOBAL_DERIV_S1)
          !Compute grad(u).s1
          DO componentIdx1=1,numberOfDimensions
            dofValue=dofValue+gradientAnalyticValue(componentIdx1)*tangents(componentIdx1,1)
          ENDDO !componentIdx1
        CASE(GLOBAL_DERIV_S2)
          !Compute grad(u).s2
          DO componentIdx1=1,numberOfDimensions
            dofValue=dofValue+gradientAnalyticValue(componentIdx1)*tangents(componentIdx1,2)
          ENDDO !componentIdx1
        CASE(GLOBAL_DERIV_S3)
          !Compute grad(u).s3
          DO componentIdx1=1,numberOfDimensions
            dofValue=dofValue+gradientAnalyticValue(componentIdx1)*tangents(componentIdx1,3)
          ENDDO !componentIdx1
        CASE(GLOBAL_DERIV_S1_S2)
          !Compute grad s1^t.H(u).s2
          DO componentIdx2=1,numberOfDimensions
            DO componentIdx1=1,numberOfDimensions
              dofValue=dofValue+ &
                & tangents(componentIdx1,1)*hessianAnalyticValue(componentIdx1,componentIdx2)*tangents(componentIdx2,2)
            ENDDO !componentIdx1
          ENDDO !componentIdx2
        CASE(GLOBAL_DERIV_S1_S3)
          !Compute grad s1^t.H(u).s3
          DO componentIdx2=1,numberOfDimensions
            DO componentIdx1=1,numberOfDimensions
              dofValue=dofValue+ &
                & tangents(componentIdx1,1)*hessianAnalyticValue(componentIdx1,componentIdx2)*tangents(componentIdx2,3)
            ENDDO !componentIdx1
          ENDDO !componentIdx2
        CASE(GLOBAL_DERIV_S2_S3)
          !Compute grad s2^t.H(u).s3
          DO componentIdx2=1,numberOfDimensions
            DO componentIdx1=1,numberOfDimensions
              dofValue=dofValue+ &
                & tangents(componentIdx1,2)*hessianAnalyticValue(componentIdx1,componentIdx2)*tangents(componentIdx2,3)
            ENDDO !componentIdx1
          ENDDO !componentIdx2
        CASE(GLOBAL_DERIV_S1_S2_S3)
!!TODO: not implemented
          dofValue=0.0_DP
        CASE DEFAULT
          localError="The global derivative index of "//TRIM(NumberToVString(derivativeGlobalIndex,"*",err,error))// &
            & " for derivative index "//TRIM(NumberToVString(derivativeIdx,"*",err,error))//" of local node number "// &
            & TRIM(NumberToVString(localNodeNumber,"*",err,error))//" of component number "// &
            & TRIM(NumberToVString(componentNumber,"*",err,error))//" of variable type "// &
            & TRIM(NumberToVString(variableType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(FIELD_DELUDELN_VARIABLE_TYPE)
        SELECT CASE(derivativeGlobalIndex)
        CASE(NO_GLOBAL_DERIV)
          !Compute grad(u).n
          DO componentIdx1=1,numberOfDimensions
            dofValue=dofValue+gradientAnalyticValue(componentIdx1)*normal(componentIdx1)
          ENDDO !componentIdx1
          dofValue=analyticValue
        CASE(GLOBAL_DERIV_S1)
            !Compute grad s1^t.H(u).n
          DO componentIdx2=1,numberOfDimensions
            DO componentIdx1=1,numberOfDimensions
              dofValue=dofValue+ &
                & tangents(componentIdx1,1)*hessianAnalyticValue(componentIdx1,componentIdx2)*normal(componentIdx2)
            ENDDO !componentIdx1
          ENDDO !componentIdx2
        CASE(GLOBAL_DERIV_S2)
          !Compute grad s2^t.H(u).n
          DO componentIdx2=1,numberOfDimensions
            DO componentIdx1=1,numberOfDimensions
              dofValue=dofValue+ &
                & tangents(componentIdx1,2)*hessianAnalyticValue(componentIdx1,componentIdx2)*normal(componentIdx2)
            ENDDO !componentIdx1
          ENDDO !componentIdx2
        CASE(GLOBAL_DERIV_S3)
          !Compute grad s3^t.H(u).n
          DO componentIdx2=1,numberOfDimensions
            DO componentIdx1=1,numberOfDimensions
              dofValue=dofValue+ &
                  & tangents(componentIdx1,3)*hessianAnalyticValue(componentIdx1,componentIdx2)*normal(componentIdx2)
            ENDDO !componentIdx1
          ENDDO !componentIdx2
        CASE(GLOBAL_DERIV_S1_S2)
!!TODO: not implemented
          dofValue=0.0_DP
        CASE(GLOBAL_DERIV_S1_S3)
!!TODO: not implemented
          dofValue=0.0_DP
        CASE(GLOBAL_DERIV_S2_S3)
!!TODO: not implemented
          dofValue=0.0_DP
        CASE(GLOBAL_DERIV_S1_S2_S3)
!!TODO: not implemented
          dofValue=0.0_DP
        CASE DEFAULT
          localError="The global derivative index of "//TRIM(NumberToVString(derivativeGlobalIndex,"*",err,error))// &
            & " for derivative index "//TRIM(NumberToVString(derivativeIdx,"*",err,error))//" of local node number "// &
            & TRIM(NumberToVString(localNodeNumber,"*",err,error))//" of component number "// &
            & TRIM(NumberToVString(componentNumber,"*",err,error))//" of variable type "// &
            & TRIM(NumberToVString(variableType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        !Do nothing or error???
      END SELECT
      CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,localDOFIdx,dofValue, &
        & err,error,*999)
      IF(setVelocity) THEN
        CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_ANALYTIC_VELOCITY_VALUES_SET_TYPE,localDOFIdx, &
          & analyticVelocity,err,error,*999)
        IF(setAcceleration) THEN
          CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_ANALYTIC_ACCELERATION_VALUES_SET_TYPE, &
            & localDOFIdx,analyticAcceleration,err,error,*999)
        ENDIF
      ENDIF
      CALL BoundaryConditionsVariable_ConditionTypeGet(boundaryConditionsVariable,localDOFIdx,conditionType,err,error,*999)
      IF(conditionType==BOUNDARY_CONDITION_FIXED) THEN
        CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_VALUES_SET_TYPE,localDOFIdx,dofValue,err,error,*999)
        IF(setVelocity) THEN
          CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_VELOCITY_VALUES_SET_TYPE,localDOFIdx, &
            & analyticVelocity,err,error,*999)
          IF(setAcceleration) THEN
            CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_ACCELERATION_VALUES_SET_TYPE, &
              & localDOFIdx,analyticAcceleration,err,error,*999)
          ENDIF
        ENDIF
      ENDIF
    ENDDO !derivativeIdx
    
    EXITS("BoundaryConditions_UpdateAnalyticNode")
    RETURN
999 ERRORSEXITS("BoundaryConditions_UpdateAnalyticNode",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_UpdateAnalyticNode

  !
  !================================================================================================================================
  !

  !>Finalise the boundary conditions variable and deallocate all memory.
  SUBROUTINE BoundaryCondition_VariableFinalise(boundaryConditionsVariable,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("BoundaryCondition_VariableFinalise",err,error,*999)

    IF(ASSOCIATED(boundaryConditionsVariable)) THEN
      IF(ALLOCATED(boundaryConditionsVariable%conditionTypes)) DEALLOCATE(boundaryConditionsVariable%conditionTypes)
      IF(ALLOCATED(boundaryConditionsVariable%DOFTypes)) DEALLOCATE(boundaryConditionsVariable%DOFTypes)
      IF(ASSOCIATED(boundaryConditionsVariable%dirichletBoundaryConditions)) THEN
        CALL BoundaryConditionsSparsityIndices_ArrayFinalise(boundaryConditionsVariable%dirichletBoundaryConditions% &
          & linearSparsityIndices,err,error,*999)
        CALL BoundaryConditionsSparsityIndices_ArrayFinalise(boundaryConditionsVariable%dirichletBoundaryConditions% &
          & dynamicSparsityIndices,err,error,*999)
        IF(ALLOCATED(boundaryConditionsVariable%dirichletBoundaryConditions%dirichletDOFIndices)) &
          & DEALLOCATE(boundaryConditionsVariable%dirichletBoundaryConditions%dirichletDOFIndices)
        DEALLOCATE(boundaryConditionsVariable%dirichletBoundaryConditions)
      ENDIF
      CALL BoundaryConditionsVariable_NeumannFinalise(boundaryConditionsVariable%neumannBoundaryConditions,err,error,*999)
      IF(ASSOCIATED(boundaryConditionsVariable%pressureIncrementedBoundaryConditions)) &
        & DEALLOCATE(boundaryConditionsVariable%pressureIncrementedBoundaryConditions)
      IF(ASSOCIATED(boundaryConditionsVariable%dofConstraints)) THEN
        CALL BoundaryConditionsDOFConstraints_Finalise(boundaryConditionsVariable%dofConstraints,err,error,*999)
        DEALLOCATE(boundaryConditionsVariable%dofConstraints)
      END IF
      DEALLOCATE(boundaryConditionsVariable)
    ENDIF

    EXITS("BoundaryCondition_VariableFinalise")
    RETURN
999 ERRORSEXITS("BoundaryCondition_VariableFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryCondition_VariableFinalise


  !
  !================================================================================================================================
  !

  !>Adds a boundary conditions variable if that variable has not already been added, otherwise do nothing.
  SUBROUTINE BoundaryConditions_VariableAdd(boundaryConditions,fieldVariable,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to add a variable for.
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to add to  the boundary conditions variables.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,numberOfGlobalDOFs,variableIdx
    TYPE(BoundaryConditionsVariablePtrType), ALLOCATABLE :: newBoundaryConditionsVariables(:)
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("BoundaryConditions_VariableAdd",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions is not associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*998)
#endif
    
    !Check if boundary conditions variable has already been added, if so then we don't do anything as different equations
    !sets can have the same dependent field variables and will both want to add the variable
    NULLIFY(boundaryConditionsVariable)
    CALL BoundaryConditions_VariableExists(boundaryConditions,fieldVariable,boundaryConditionsVariable,err,error,*999)
    IF(.NOT.ASSOCIATED(boundaryConditionsVariable)) THEN
      ALLOCATE(newBoundaryConditionsVariables(boundaryConditions%numberOfBoundaryConditionsVariables+1),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new boundary conditions variables array.",err,error,*998)
      IF(ALLOCATED(boundaryConditions%boundaryConditionsVariables)) THEN
        DO variableIdx=1,boundaryConditions%numberOfBoundaryConditionsVariables
          newBoundaryConditionsVariables(variableIdx)%ptr=>boundaryConditions%boundaryConditionsVariables(variableIdx)%ptr
        ENDDO !variableIdx
      ENDIF

      NULLIFY(newBoundaryConditionsVariables(boundaryConditions%numberOfBoundaryConditionsVariables+1)%ptr)
      CALL BoundaryConditionsVariable_Initialise(newBoundaryConditionsVariables(boundaryConditions% &
        & numberOfBoundaryConditionsVariables+1)%ptr,err,error,*999)
      boundaryConditionsVariable=>newBoundaryConditionsVariables(boundaryConditions%numberOfBoundaryConditionsVariables+1)%ptr
      boundaryConditionsVariable%boundaryConditions=>boundaryConditions
      boundaryConditionsVariable%variableType=fieldVariable%variableType
      boundaryConditionsVariable%variable=>fieldVariable
      CALL FieldVariable_NumberOfGlobalDOFsGet(fieldVariable,numberOfGlobalDOFs,err,error,*999)
      ALLOCATE(boundaryConditionsVariable%conditionTypes(numberOfGlobalDOFs),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate global boundary condition types.",err,error,*999)
      ALLOCATE(boundaryConditionsVariable%DOFTypes(numberOfGlobalDOFs),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate global boundary condition dof types.",err,error,*999)
      boundaryConditionsVariable%conditionTypes=BOUNDARY_CONDITION_NONE
      boundaryConditionsVariable%DOFTypes=BOUNDARY_CONDITION_DOF_FREE
      ALLOCATE(boundaryConditionsVariable%dofCounts(MAX_BOUNDARY_CONDITION_NUMBER),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate boundary condition DOF counts array.",err,error,*999)
      boundaryConditionsVariable%dofCounts=0
      NULLIFY(boundaryConditionsVariable%dirichletBoundaryConditions)
      boundaryConditionsVariable%numberOfDirichletConditions=0
      NULLIFY(boundaryConditionsVariable%neumannBoundaryConditions)
      NULLIFY(boundaryConditionsVariable%pressureIncrementedBoundaryConditions)
      ALLOCATE(boundaryConditionsVariable%parameterSetRequired(FIELD_NUMBER_OF_SET_TYPES),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate boundary condition parameter set required array.",err,error,*999)
      boundaryConditionsVariable%parameterSetRequired=.FALSE.
      boundaryConditionsVariable%parameterSetRequired(FIELD_VALUES_SET_TYPE)=.TRUE.
      
      CALL MOVE_ALLOC(newBoundaryConditionsVariables,boundaryConditions%boundaryConditionsVariables)
      boundaryConditions%numberOfBoundaryConditionsVariables=boundaryConditions%numberOfBoundaryConditionsVariables+1

      NULLIFY(boundaryConditionsVariable%DOFConstraints)
      CALL BoundaryConditionsDOFConstraints_Initialise(boundaryConditionsVariable%DOFconstraints,err,error,*999)

    ENDIF

    EXITS("BoundaryConditions_VariableAdd")
    RETURN
999 CALL BoundaryCondition_VariableFinalise(boundaryConditionsVariable,dummyErr,dummyError,*998)
    IF(ALLOCATED(newBoundaryConditionsVariables)) DEALLOCATE(newBoundaryConditionsVariables)
998 ERRORSEXITS("BoundaryConditions_VariableAdd",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_VariableAdd

  !
  !================================================================================================================================
  !

  !>Find the boundary conditions variable for a given field variable if it exists
  SUBROUTINE BoundaryConditions_VariableExists(boundaryConditions,fieldVariable,boundaryConditionsVariable,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to see if the boundary conditions variable exists for
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to check the boundary conditions variable for.
    TYPE(BoundaryConditionsVariableType), POINTER, INTENT(OUT) :: boundaryConditionsVariable !<On return, a pointer to the boundary conditions variable, or NULL if it wasn't found. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: fieldInterfaceUserNumber,fieldVarType,fieldVariableFieldUserNumber,fieldVariableRegionUserNumber, &
      & matchInterfaceUserNumber,matchVariableType,matchVariableFieldUserNumber,matchVariableRegionUserNumber, &
      & numberOfBoundaryConditionsVariables,variableIdx
    LOGICAL :: fieldVariableIsRegion,matchVariableIsInterface,matchVariableIsRegion,variableFound
    TYPE(FieldType), POINTER :: fieldVariableField,matchVariableField
    TYPE(FieldVariableType), POINTER :: matchVariable
    TYPE(InterfaceType), POINTER :: fieldVariableInterface,matchVariableInterface
    TYPE(RegionType), POINTER :: fieldInterfaceParentRegion,fieldVariableRegion,matchInterfaceParentRegion,matchVariableRegion
    TYPE(VARYING_STRING) :: localError

    ENTERS("BoundaryConditions_VariableExists",err,error,*998)

    IF(ASSOCIATED(boundaryConditionsVariable)) CALL FlagError("Boundary conditions variable is already associated.",err,error,*998)

    CALL FieldVariable_VariableTypeGet(fieldVariable,fieldVarType,err,error,*999)
    NULLIFY(fieldVariableField)
    CALL FieldVariable_FieldGet(fieldVariable,fieldVariableField,err,error,*999)
    CALL Field_UserNumberGet(fieldVariableField,fieldVariableFieldUserNumber,err,error,*999)
    CALL Field_IsRegionField(fieldVariableField,fieldVariableIsRegion,err,error,*999)
    NULLIFY(fieldVariableRegion)
    NULLIFY(fieldVariableInterface)
    IF(fieldVariableIsRegion) THEN
      CALL Field_RegionGet(fieldVariableField,fieldVariableRegion,err,error,*999)
      CALL Region_UserNumberGet(fieldVariableRegion,fieldVariableRegionUserNumber,err,error,*999)
    ELSE
      CALL Field_InterfaceGet(fieldVariableField,fieldVariableInterface,err,error,*999)
      CALL Interface_UserNumberGet(fieldVariableInterface,fieldInterfaceUserNumber,err,error,*999)
      NULLIFY(fieldInterfaceParentRegion)
      CALL Interface_ParentRegionGet(fieldVariableInterface,fieldInterfaceParentRegion,err,error,*999)
      CALL Region_UserNumberGet(fieldInterfaceParentRegion,fieldVariableRegionUserNumber,err,error,*999)      
    ENDIF
        
    CALL BoundaryConditions_NumberOfVariablesGet(boundaryConditions,numberOfBoundaryConditionsVariables,err,error,*999)
    variableFound=.FALSE.
    variableIdx=1
    DO WHILE(variableIdx<=numberOfBoundaryConditionsVariables.AND..NOT.variableFound)
      NULLIFY(boundaryConditionsVariable)
      CALL BoundaryConditions_VariableIndexGet(boundaryConditions,variableIdx,boundaryConditionsVariable,err,error,*999)
      NULLIFY(matchVariable)
      CALL BoundaryConditionsVariable_FieldVariableGet(boundaryConditionsVariable,matchVariable,err,error,*999)
      CALL FieldVariable_VariableTypeGet(matchVariable,matchVariableType,err,error,*999)
      NULLIFY(matchVariableField)
      CALL FieldVariable_FieldGet(matchVariable,matchVariableField,err,error,*999)
      CALL Field_UserNumberGet(matchVariableField,matchVariableFieldUserNumber,err,error,*999)
      CALL Field_IsRegionField(matchVariableField,matchVariableIsRegion,err,error,*999)
      IF(fieldVariableIsRegion) THEN
        IF(matchVariableIsRegion) THEN
          NULLIFY(matchVariableRegion)
          CALL Field_RegionGet(matchVariableField,matchVariableRegion,err,error,*999)
          CALL Region_UserNumberGet(matchVariableRegion,matchVariableRegionUserNumber,err,error,*999)
          IF(fieldVariableRegionUserNumber==matchVariableRegionUserNumber.AND. &
            & fieldVariableFieldUserNumber==matchVariableFieldUserNumber.AND. &
            & fieldVarType==matchVariableType) variableFound=.TRUE.
        ENDIF
      ELSE
        CALL Field_IsInterfaceField(matchVariableField,matchVariableIsInterface,err,error,*999)
        IF(matchVariableIsInterface) THEN
          NULLIFY(matchVariableInterface)
          CALL Field_InterfaceGet(matchVariableField,matchVariableInterface,err,error,*999)
          CALL Interface_UserNumberGet(matchVariableInterface,matchInterfaceUserNumber,err,error,*999)
          NULLIFY(matchInterfaceParentRegion)
          CALL Interface_ParentRegionGet(matchVariableInterface,matchInterfaceParentRegion,err,error,*999)
          CALL Region_UserNumberGet(matchInterfaceParentRegion,matchVariableRegionUserNumber,err,error,*999)
          IF(fieldInterfaceUserNumber==matchInterfaceUserNumber.AND. &
            & fieldVariableRegionUserNumber==matchVariableRegionUserNumber.AND. &
            & fieldVariableFieldUserNumber==matchVariableFieldUserNumber.AND. &
            & fieldVarType==matchVariableType) variableFound=.TRUE.
        ENDIF
      ENDIF
      variableIdx=variableIdx+1
    ENDDO
    IF(.NOT.variableFound) NULLIFY(boundaryConditionsVariable)

    EXITS("BoundaryConditions_VariableExists")
    RETURN
999 NULLIFY(boundaryConditionsVariable)
998 ERRORSEXITS("BoundaryConditions_VariableExists",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_VariableExists

  !
  !================================================================================================================================
  !

  !>Find the boundary conditions variable for a given field variable
  SUBROUTINE BoundaryConditions_VariableGet(boundaryConditions,fieldVariable,boundaryConditionsVariable,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to get the boundary conditions variable for.
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to get the boundary conditions variable for.
    TYPE(BoundaryConditionsVariableType), POINTER, INTENT(OUT) :: boundaryConditionsVariable !<On return, a pointer to the boundary conditions variable. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("BoundaryConditions_VariableGet",err,error,*998)

    IF(ASSOCIATED(boundaryConditionsVariable)) CALL FlagError("Boundary conditions variable is already associated.",err,error,*998)
    
    CALL BoundaryConditions_VariableExists(boundaryConditions,fieldVariable,boundaryConditionsVariable,err,error,*999)    

    IF(.NOT.ASSOCIATED(boundaryConditionsVariable)) THEN
      localError="Boundary condition variable is not associated"
      IF(ASSOCIATED(fieldVariable)) THEN
        localError=localError//" for variable type "//TRIM(NumberToVString(fieldVariable%variableType,"*",err,error))
        IF(ASSOCIATED(fieldVariable%field)) THEN
          localError=localError//" of field number "//TRIM(NumberToVString(fieldVariable%field%userNumber,"*",err,error))
          IF(ASSOCIATED(fieldVariable%field%region)) THEN
            localError=localError//" of region number "//TRIM(NumberToVString(fieldVariable%field%region%userNumber,"*",err,error))
          ELSE IF(ASSOCIATED(fieldVariable%field%interface)) THEN
            localError=localError//" of interface number "// &
              & TRIM(NumberToVString(fieldVariable%field%interface%userNumber,"*",err,error))
            IF(ASSOCIATED(fieldVariable%field%interface%parentRegion)) &
              & localError=localError//" of parent region number "// &
              & TRIM(NumberToVString(fieldVariable%field%interface%parentRegion%userNumber,"*",err,error))
          ENDIF
        ENDIF
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    EXITS("BoundaryConditions_VariableGet")
    RETURN
999 NULLIFY(boundaryConditionsVariable)
998 ERRORSEXITS("BoundaryConditions_VariableGet",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_VariableGet

  !
  !================================================================================================================================
  !

  !>Finalise the coupled DOFs structure and deallocates all memory.
  SUBROUTINE BoundaryConditionsCoupledDOFs_Finalise(coupledDOFs,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsCoupledDOFsType), POINTER :: coupledDOFs !<A pointer to the coupled DOFs to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables

    ENTERS("BoundaryConditionsCoupledDOFs_Finalise",err,error,*999)

    IF(ASSOCIATED(coupledDOFs)) THEN
      IF(ALLOCATED(coupledDOFs%globalDOFs)) DEALLOCATE(coupledDOFs%globalDOFs)
      IF(ALLOCATED(coupledDOFs%localDOFs)) DEALLOCATE(coupledDOFs%localDOFs)
      IF(ALLOCATED(coupledDOFs%coefficients)) DEALLOCATE(coupledDOFs%coefficients)
      DEALLOCATE(coupledDOFs)
    ENDIF

    EXITS("BoundaryConditionsCoupledDOFs_Finalise")
    RETURN
999 ERRORSEXITS("BoundaryConditionsCoupledDOFs_Finalise",err,error)
    RETURN 1

  END SUBROUTINE BoundaryConditionsCoupledDOFs_Finalise

  !
  !================================================================================================================================
  !

  !>Allocates and initialises the coupled DOFs structure
  SUBROUTINE BoundaryConditionsCoupledDOFs_Initialise(coupledDOFs,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsCoupledDOFsType), POINTER :: coupledDOFs !<A pointer to the coupled DOFs to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("BoundaryConditionsCoupledDOFs_Initialise",err,error,*998)

    IF(ASSOCIATED(coupledDOFs)) CALL FlagError("The DOF constraint pointer is already associated.",err,error,*998)

    ALLOCATE(coupledDOFs,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate the coupled DOFs.",err,error,*999)
    coupledDOFs%numberOfDOFs=0

    EXITS("BoundaryConditionsCoupledDOFs_Initialise")
    RETURN
999 CALL BoundaryConditionsCoupledDOFs_Finalise(coupledDOFs,dummyErr,dummyError,*998)
998 ERRORSEXITS("BoundaryConditionsCoupledDOFs_Initialise",err,error)
    RETURN 1

  END SUBROUTINE BoundaryConditionsCoupledDOFs_Initialise

  !
  !================================================================================================================================
  !

  !>Finalise the DOF constraint structure and deallocates all memory.
  SUBROUTINE BoundaryConditionsDOFConstraint_Finalise(dofConstraint,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsDofConstraintType), POINTER :: dofConstraint !<A pointer to the dof constraint to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables

    ENTERS("BoundaryConditionsDOFConstraint_Finalise",err,error,*999)

    IF(ASSOCIATED(dofConstraint)) THEN
      IF(ALLOCATED(dofConstraint%dofs)) DEALLOCATE(dofConstraint%dofs)
      IF(ALLOCATED(dofConstraint%coefficients)) DEALLOCATE(dofConstraint%coefficients)
      DEALLOCATE(dofConstraint)
    ENDIF

    EXITS("BoundaryConditionsDOFConstraint_Finalise")
    RETURN
999 ERRORSEXITS("BoundaryConditionsDOFConstraint_Finalise",err,error)
    RETURN 1

  END SUBROUTINE BoundaryConditionsDOFConstraint_Finalise

  !
  !================================================================================================================================
  !

  !>Allocates and initialises the DOF constraint structure
  SUBROUTINE BoundaryConditionsDOFConstraint_Initialise(dofConstraint,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsDOFConstraintType), POINTER :: dofConstraint !<A pointer to the DOF constraint to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("BoundaryConditionsDOFConstraint_Initialise",err,error,*998)

    IF(ASSOCIATED(dofConstraint)) CALL FlagError("The DOF constraint pointer is already associated.",err,error,*998)

    ALLOCATE(dofConstraint,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate the DOF constraint.",err,error,*999)
    dofConstraint%numberOfDOFs=0
    dofConstraint%globalDOF=0

    EXITS("BoundaryConditionsDOFConstraint_Initialise")
    RETURN
999 CALL BoundaryConditionsDOFConstraint_Finalise(dofConstraint,dummyErr,dummyError,*998)
998 ERRORSEXITS("BoundaryConditionsDOFConstraint_Initialise",err,error)
    RETURN 1

  END SUBROUTINE BoundaryConditionsDOFConstraint_Initialise

  !
  !================================================================================================================================
  !

  !>Finalise the DOF constraints structure and deallocate all memory.
  SUBROUTINE BoundaryConditionsDOFConstraints_Finalise(dofConstraints,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsDofConstraintsType), POINTER :: dofConstraints !<A pointer to the dof constraints to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: constraintIdx,dofIdx

    ENTERS("BoundaryConditionsDOFConstraints_Finalise",err,error,*999)

    IF(ASSOCIATED(dofConstraints)) THEN
      IF(ALLOCATED(dofConstraints%constraints)) THEN
        DO constraintIdx=1,SIZE(dofConstraints%constraints,1)
          CALL BoundaryConditionsDOFConstraint_Finalise(dofConstraints%constraints(constraintIdx)%ptr,err,error,*999)
        ENDDO !constraintIdx
        DEALLOCATE(dofConstraints%constraints)
      ENDIF
      IF(ALLOCATED(dofConstraints%dofCouplings)) THEN
        DO dofIdx=1,SIZE(dofConstraints%dofCouplings,1)
          CALL BoundaryConditionsCoupledDOFs_Finalise(dofConstraints%dofCouplings(dofIdx)%ptr,err,error,*999)
        ENDDO !dofIdx
        DEALLOCATE(dofConstraints%dofCouplings)
      ENDIF
      DEALLOCATE(dofConstraints)
    ENDIF

    EXITS("BoundaryConditionsDOFConstraints_Finalise")
    RETURN
999 ERRORSEXITS("BoundaryConditionsDOFConstraints_Finalise",err,error)
    RETURN 1

  END SUBROUTINE BoundaryConditionsDOFConstraints_Finalise

  !
  !================================================================================================================================
  !

  !>Allocates and initialise the DOF constraints structure
  SUBROUTINE BoundaryConditionsDOFConstraints_Initialise(dofConstraints,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsDofConstraintsType), POINTER :: dofConstraints !<A pointer to the DOF constraints to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("BoundaryConditionsDOFConstraints_Initialise",err,error,*998)

    IF(ASSOCIATED(dofConstraints)) CALL FlagError("The DOF constraints pointer is already associated.",err,error,*998)

    ALLOCATE(dofConstraints,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate the DOF constraints.",err,error,*999)
    dofConstraints%numberOfConstraints=0
    dofConstraints%numberOfDOFs=0

    EXITS("BoundaryConditionsDOFConstraints_Initialise")
    RETURN
999 CALL BoundaryConditionsDOFConstraints_Finalise(dofConstraints,dummyErr,dummyError,*998)
998 ERRORSEXITS("BoundaryConditionsDOFConstraints_Initialise",err,error)
    RETURN 1

  END SUBROUTINE BoundaryConditionsDOFConstraints_Initialise

  !
  !================================================================================================================================
  !

  !Finalises a boundary conditions row variable and deallocates all memory.
  SUBROUTINE BoundaryConditionsRowVariable_Finalise(boundaryConditionsRowVariable,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsRowVariableType), POINTER :: boundaryConditionsRowVariable !<A pointer to the boundary conditions row variable to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables

    ENTERS("BoundaryConditionsRowVariable_Finalise",err,error,*999)

    IF(ASSOCIATED(boundaryConditionsRowVariable)) THEN
      IF(ALLOCATED(boundaryConditionsRowVariable%rowConditionTypes)) DEALLOCATE(boundaryConditionsRowVariable%rowConditionTypes)
      DEALLOCATE(boundaryConditionsRowVariable)
    ENDIF
    
    EXITS("BoundaryConditionsRowVariable_Finalise")
    RETURN
999 ERRORSEXITS("BoundaryConditionsRowVariable_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsRowVariable_Finalise

  !
  !================================================================================================================================
  !

  !Allocates and initialises a boundary conditions row variable.
  SUBROUTINE BoundaryConditionsRowVariable_Initialise(boundaryConditionsRowVariable,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsRowVariableType), POINTER :: boundaryConditionsRowVariable !<A pointer to the boundary conditions row variable to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("BoundaryConditionsRowVariable_Initialise",err,error,*998)

    IF(ASSOCIATED(boundaryConditionsRowVariable)) &
      & CALL FlagError("Boundary conditions row variable is already associated.",err,error,*998)

    ALLOCATE(boundaryConditionsRowVariable,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocated boundary conditions row variable.",err,error,*999)
    NULLIFY(boundaryConditionsRowVariable%boundaryConditions)
    NULLIFY(boundaryConditionsRowVariable%lhsVariable)
    NULLIFY(boundaryConditionsRowVariable%rhsVariable)
    boundaryConditionsRowVariable%numberOfDOFs=0
    boundaryConditionsRowVariable%totalNumberOfDOFs=0    
    
    EXITS("BoundaryConditionsRowVariable_Initialise")
    RETURN
999 CALL BoundaryConditionsRowVariable_Finalise(boundaryConditionsRowVariable,dummyErr,dummyError,*998)
998 ERRORSEXITS("BoundaryConditionsRowVariable_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsRowVariable_Initialise

  !
  !================================================================================================================================
  !

  !>Sets the boundary condition row type  for the row boundary conditions.
  SUBROUTINE BoundaryConditionsRowVariable_RowConditionTypeSet(boundaryConditionsRowVariable,dofIdx,rowCondition,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsRowVariableType), POINTER :: boundaryConditionsRowVariable !<A pointer to the boundary conditions row variable to set the row condition for
    INTEGER(INTG), INTENT(IN) :: dofIdx !<The DOF index to set the row boundary condition at
    INTEGER(INTG), INTENT(IN) :: rowCondition !<The row condition type to set \see BoundaryConditionsRoutines_RowTypes,BoundaryConditionsRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: currentRowCondition,rowConditionType,variableType
    TYPE(FieldVariableType), POINTER :: rowVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("BoundaryConditionsRowVariable_RowConditionTypeSet",err,error,*999)

    NULLIFY(rowVariable)
    CALL BoundaryConditionsRowVariable_LHSVariableGet(boundaryConditionsRowVariable,rowVariable,err,error,*999)
    CALL BoundaryConditionsRowVariable_RowConditionTypeGet(boundaryConditionsRowVariable,dofIdx,currentRowCondition,err,error,*999)

    IF(rowCondition/=currentRowCondition) THEN
      !Check that the row boundary condition has not already been set.
      IF(currentRowCondition/=BOUNDARY_CONDITION_FREE_ROW) THEN        
        localError="The row boundary condition has already been set. The current row boundary condition is "// &
          & TRIM(NumberToVString(currentRowCondition,"*",err,error))//" for DOF index "// &
          & TRIM(NumberToVString(dofIdx,"*",err,error))//" of variable type "// &
          & TRIM(NumberToVString(rowVariable%variableType,"*",err,error))
        IF(ASSOCIATED(rowVariable%field)) THEN
          localError=localError//" of field number "//TRIM(NumberToVString(rowVariable%field%userNumber,"*",err,error))
          IF(ASSOCIATED(rowVariable%field%region)) THEN
            localError=localError//" of region number "//TRIM(NumberToVString(rowVariable%field%region%userNumber,"*",err,error))
          ELSE IF(ASSOCIATED(rowVariable%field%INTERFACE)) THEN
            localError=localError//" of interface number "// &
              & TRIM(NumberToVString(rowVariable%field%interface%userNumber,"*",err,error))            
          ENDIF
        ENDIF
        localError=localError//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      !Check condition
      SELECT CASE(rowCondition)
      CASE(BOUNDARY_CONDITION_FREE_ROW)
        rowConditionType=BOUNDARY_CONDITION_FREE_ROW
      CASE(BOUNDARY_CONDITION_DIRICHLET_ROW)
        rowConditionType=BOUNDARY_CONDITION_DIRICHLET_ROW
      CASE(BOUNDARY_CONDITION_POINT_NEUMANN_ROW)
        rowConditionType=BOUNDARY_CONDITION_POINT_NEUMANN_ROW
      CASE(BOUNDARY_CONDITION_INTEGRATED_NEUMANN_ROW)
        rowConditionType=BOUNDARY_CONDITION_INTEGRATED_NEUMANN_ROW
      CASE(BOUNDARY_CONDITION_ROBIN_ROW)
        rowConditionType=BOUNDARY_CONDITION_ROBIN_ROW
      CASE(BOUNDARY_CONDITION_CAUCHY_ROW)
        rowConditionType=BOUNDARY_CONDITION_CAUCHY_ROW
      CASE(BOUNDARY_CONDITION_CONSTRAINED_ROW)
        rowConditionType=BOUNDARY_CONDITION_CONSTRAINED_ROW
      CASE DEFAULT
        localError="The specified boundary row condition type of "//TRIM(NumberToVString(rowCondition,"*",err,error))// &
          & "for DOF index "//TRIM(NumberToVString(dofIdx,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Set the row condition type
      boundaryConditionsRowVariable%rowConditionTypes(dofIdx)=rowConditionType
    ENDIF    
    
    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Row boundary condition being set:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  DOF index     = ", dofIdx,err,error,*999)
      CALL FieldVariable_VariableTypeGet(rowVariable,variableType,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Variable type = ",variableType,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Old condition = ",currentRowCondition,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  New condition = ",rowCondition,err,error,*999)
    ENDIF
    
    EXITS("BoundaryConditionsRowVariable_RowConditionTypeSet")
    RETURN
999 ERRORSEXITS("BoundaryConditionsRowVariable_RowConditionTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsRowVariable_RowConditionTypeSet

  !
  !================================================================================================================================
  !

  !>Finalise an array of sparcity indices and deallocate all memory.
  SUBROUTINE BoundaryConditionsSparsityIndices_ArrayFinalise(sparsityIndicesArray,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsSparsityIndicesPtrType), ALLOCATABLE :: sparsityIndicesArray(:,:)
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx, equationsMatrixIdx

    ENTERS("BoundaryConditionsSparsityIndices_ArrayFinalise",err,error,*999)

    IF(ALLOCATED(sparsityIndicesArray)) THEN
      DO equationsSetIdx=1,SIZE(sparsityIndicesArray,1)
        DO equationsMatrixIdx=1,SIZE(sparsityIndicesArray,2)
          CALL BoundaryConditionsSparsityIndices_Finalise(sparsityIndicesArray(equationsSetIdx,equationsMatrixIdx)%ptr, &
            & err,error,*999)
        ENDDO !equationsMatrixIdx
      ENDDO !equtionsSetIdx
      DEALLOCATE(sparsityIndicesArray)
    ENDIF

    EXITS("BoundaryConditionsSparsityIndices_ArrayFinalise")
    RETURN
999 ERRORS("BoundaryConditionsSparsityIndices_ArrayFinalise",err,error)
    EXITS("BoundaryConditionsSparsityIndices_ArrayFinalise")
    RETURN 1

  END SUBROUTINE BoundaryConditionsSparsityIndices_ArrayFinalise
  
  !
  !================================================================================================================================
  !

  !>Finalises the boundary condition sparsity indices and deallocates all memory
  SUBROUTINE BoundaryConditionsSparsityIndices_Finalise(sparsityIndices,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsSparsityIndicesType), POINTER :: sparsityIndices !<A pointer to the sparsity indices to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("BoundaryConditionsSparsityIndices_Finalise",err,error,*999)

    IF(ASSOCIATED(sparsityIndices)) THEN
      IF(ALLOCATED(sparsityIndices%sparseRowIndices)) DEALLOCATE(sparsityIndices%sparseRowIndices)
      IF(ALLOCATED(sparsityIndices%sparseColumnIndices)) DEALLOCATE(sparsityIndices%sparseColumnIndices)
      DEALLOCATE(sparsityIndices)
    ENDIF
    
    EXITS("BoundaryConditionsSparsityIndices_Finalise")
    RETURN
999 ERRORS("BoundaryConditionsSparsityIndices_Finalise",err,error)
    EXITS("BoundaryConditionsSparsityIndices_Finalise")
    RETURN 1

  END SUBROUTINE BoundaryConditionsSparsityIndices_Finalise

  !
  !================================================================================================================================
  !

  !>Initialise Sparsity Indices type
  SUBROUTINE BoundaryConditionsSparsityIndices_Initialise(sparsityIndices,numberOfDirichlet,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsSparsityIndicesType), POINTER :: sparsityIndices !<A pointer to the Sparsity Indices type tp initialise
    INTEGER(INTG), INTENT(IN) :: numberOfDirichlet !<The number of dirichlet conditions this sparsity indices type will hold
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("BoundaryConditionsSparsityIndices_Initialise",err,error,*998)

    IF(ASSOCIATED(sparsityIndices)) CALL FlagError("Sparsity Indices are already associated.",err,error,*998)
   
    ALLOCATE(sparsityIndices,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate sparsity indicies.",err,error,*999)
    ALLOCATE(sparsityIndices%sparseColumnIndices(numberOfDirichlet+1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate sparsity column indices array",err,error,*999)

    EXITS("BoundaryConditionsSparsityIndices_Initialise")
    RETURN
999 CALL BoundaryConditionsSparsityIndices_Finalise(sparsityIndices,dummyErr,dummyError,*998)
998 ERRORS("BoundaryConditionsSparsityIndices_Initialise",err,error)
    EXITS("BoundaryConditionsSparsityIndices_Initialise")
    RETURN 1

  END SUBROUTINE BoundaryConditionsSparsityIndices_Initialise

  !
  !================================================================================================================================
  !

  !>Checks that the applied boundary conditions are supported by the equations sets in the solver equations
  SUBROUTINE BoundaryConditionsVariable_CheckEquations(boundaryConditionsVariable,err,error,*)

    ! Argument variables
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to check
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    type(varying_string), intent(out) :: error !<The error string
    ! Local variables
    INTEGER(INTG) :: boundaryConditionType,equationsSetIdx,esSpecification(3),numberOfEquationsSets,specificationSize
    LOGICAL :: validEquationsSetFound
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(VARYING_STRING) :: localError

    ENTERS("BoundaryConditionsVariable_CheckEquations",err,error,*999)

    NULLIFY(boundaryConditions)
    CALL BoundaryConditionsVariable_BoundaryConditionsGet(boundaryConditionsVariable,boundaryConditions,err,error,*999)
    NULLIFY(solverEquations)
    CALL BoundaryConditions_SolverEquationsGet(boundaryConditions,solverEquations,err,error,*999)
    NULLIFY(solverMapping)
    CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)

    DO boundaryConditionType=1,MAX_BOUNDARY_CONDITION_NUMBER
      !Check if any DOFs have been set to this BC type
      IF(boundaryConditionsVariable%dofCounts(boundaryConditionType)>0) THEN
        validEquationsSetFound=.FALSE.
        CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
        DO equationsSetIdx=1,numberOfEquationsSets
          NULLIFY(equationsSet)
          CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
          CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)

          SELECT CASE(boundaryConditionType)
          CASE(BOUNDARY_CONDITION_NONE)
            !Valid for any equations set
            validEquationsSetFound=.TRUE.
          CASE(BOUNDARY_CONDITION_FIXED)
            validEquationsSetFound=.TRUE.
          CASE(BOUNDARY_CONDITION_LINEAR_CONSTRAINT)
            validEquationsSetFound=.TRUE.
          CASE(BOUNDARY_CONDITION_FIXED_INLET, &
            & BOUNDARY_CONDITION_FIXED_OUTLET)
            IF(esSpecification(1)==EQUATIONS_SET_FLUID_MECHANICS_CLASS.AND. &
              & (esSpecification(2)==EQUATIONS_SET_STOKES_EQUATION_TYPE.OR. &
              & esSpecification(2)==EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE.OR. &
              & esSpecification(2)==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE)) &
              & validEquationsSetFound=.TRUE.
          CASE(BOUNDARY_CONDITION_FIXED_WALL,BOUNDARY_CONDITION_MOVED_WALL, &
            & BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED,BOUNDARY_CONDITION_FREE_WALL)
            IF(esSpecification(1)==EQUATIONS_SET_FLUID_MECHANICS_CLASS.AND. &
              & (esSpecification(2)==EQUATIONS_SET_STOKES_EQUATION_TYPE.OR. &
              & esSpecification(2)==EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE.OR. &
              & esSpecification(2)==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE.OR. &
              & esSpecification(2)==EQUATIONS_SET_DARCY_EQUATION_TYPE)) &
              validEquationsSetFound=.TRUE.
            IF(esSpecification(1)==EQUATIONS_SET_CLASSICAL_FIELD_CLASS.AND. &
              & esSpecification(2)==EQUATIONS_SET_LAPLACE_EQUATION_TYPE.AND. &
              & esSpecification(3)==EQUATIONS_SET_MOVING_MESH_LAPLACE_SUBTYPE) &
              & validEquationsSetFound=.TRUE.
          CASE(BOUNDARY_CONDITION_FIXED_INCREMENTED)
            validEquationsSetFound=.TRUE.
          CASE(BOUNDARY_CONDITION_PRESSURE, &
              & BOUNDARY_CONDITION_PRESSURE_INCREMENTED)
            IF(esSpecification(1)==EQUATIONS_SET_ELASTICITY_CLASS.AND. &
              & esSpecification(2)==EQUATIONS_SET_FINITE_ELASTICITY_TYPE) &
              & validEquationsSetFound=.TRUE.
            IF(esSpecification(1)==EQUATIONS_SET_FLUID_MECHANICS_CLASS .AND. &
              & (esSpecification(2)==EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE.OR. &
              & esSpecification(2)==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE)) &
              &  validEquationsSetFound=.TRUE.
          CASE(BOUNDARY_CONDITION_CORRECTION_MASS_INCREASE)
            !Not actually used anywhere? So keep it as invalid, although maybe it should be removed?
            validEquationsSetFound=.FALSE.
          CASE(BOUNDARY_CONDITION_IMPERMEABLE_WALL)
            IF(esSpecification(1)==EQUATIONS_SET_ELASTICITY_CLASS.AND. &
              & esSpecification(2)==EQUATIONS_SET_FINITE_ELASTICITY_TYPE.AND. &
              & (esSpecification(3)==EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE.OR. &
              & esSpecification(3)==EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE.OR. &
              & esSpecification(3)==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)) &
              & validEquationsSetFound=.TRUE.
          CASE(BOUNDARY_CONDITION_NEUMANN_POINT,BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)
            validEquationsSetFound=.TRUE.
          CASE(BOUNDARY_CONDITION_NEUMANN_INTEGRATED,BOUNDARY_CONDITION_NEUMANN_INTEGRATED_ONLY)
            validEquationsSetFound=.TRUE.
          CASE(BOUNDARY_CONDITION_FIXED_FITTED)
            IF(esSpecification(1)==EQUATIONS_SET_FLUID_MECHANICS_CLASS.AND. &
              & (esSpecification(2)==EQUATIONS_SET_STOKES_EQUATION_TYPE.OR. &
              & esSpecification(2)==EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE.OR. &
              & esSpecification(2)==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE)) &
              & validEquationsSetFound=.TRUE.
          CASE(BOUNDARY_CONDITION_FIXED_NONREFLECTING,BOUNDARY_CONDITION_FIXED_CELLML,BOUNDARY_CONDITION_FIXED_STREE)
            IF(esSpecification(1)==EQUATIONS_SET_FLUID_MECHANICS_CLASS.AND. &
              & (esSpecification(2)==EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE.OR. &
              & esSpecification(2)==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE)) &
              validEquationsSetFound=.TRUE.
          CASE DEFAULT
            localError="The specified boundary condition type of "//TRIM(NumberToVString(boundaryConditionType,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ENDDO !equationsSetIdx

        IF(.NOT.validEquationsSetFound) THEN
          localError="The specified boundary condition type of "//TRIM(NumberToVString(boundaryConditionType,"*",err,error))// &
            & " is invalid for the equations sets in the solver equations."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDIF
    ENDDO !boundaryConditionType

    EXITS("BoundaryConditionsVariable_CheckEquations")
    RETURN
999 ERRORSEXITS("BoundaryConditionsVariable_CheckEquations",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsVariable_CheckEquations

  !
  !================================================================================================================================
  !

  !> Checks the boundary condition type and sets the boundary condition type and dof type for the boundary conditions.
  !> Makes sure any field parameter sets required are created, and sets the parameter set required array value.
  SUBROUTINE BoundaryConditionsVariable_ConditionTypeSet(boundaryConditionsVariable,globalDof,condition,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to set the boundary condition for
    INTEGER(INTG), INTENT(IN) :: globalDof !<The globalDof to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: condition !<The boundary condition type to set \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: dofType,numberOfGlobalDOFs,previousCondition,previousDOF,variableType
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("BoundaryConditionsVariable_ConditionTypeSet",err,error,*999)

    NULLIFY(fieldVariable)
    CALL BoundaryConditionsVariable_FieldVariableGet(boundaryConditionsVariable,fieldVariable,err,error,*999)
    CALL FieldVariable_NumberOfGlobalDOFsGet(fieldVariable,numberOfGlobalDOFs,err,error,*999)
    IF(globalDOF<1.OR.globalDOF>numberOfGlobalDOFs) THEN
      localError="The specified global DOF number of "//TRIM(NumberToVString(globalDOF,"*",err,error))// &
        & " is invalid. The global number should be >= 1 and <= "//TRIM(NumberToVString(numberOfGlobalDOFs,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(condition)
    CASE(BOUNDARY_CONDITION_NONE)
      dofType=BOUNDARY_CONDITION_DOF_FREE
    CASE(BOUNDARY_CONDITION_FIXED)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
    CASE(BOUNDARY_CONDITION_LINEAR_CONSTRAINT)
      dofType=BOUNDARY_CONDITION_DOF_CONSTRAINED
    CASE(BOUNDARY_CONDITION_FIXED_INLET)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
    CASE(BOUNDARY_CONDITION_FIXED_OUTLET)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
    CASE(BOUNDARY_CONDITION_FIXED_WALL)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
    CASE(BOUNDARY_CONDITION_MOVED_WALL)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
    CASE(BOUNDARY_CONDITION_FREE_WALL)
      dofType=BOUNDARY_CONDITION_DOF_FREE
    CASE(BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
      CALL FieldVariable_ParameterSetEnsureCreated(fieldVariable,FIELD_BOUNDARY_CONDITIONS_SET_TYPE,err,error,*999)
      boundaryConditionsVariable%parameterSetRequired(FIELD_BOUNDARY_CONDITIONS_SET_TYPE)=.TRUE.
    CASE(BOUNDARY_CONDITION_FIXED_INCREMENTED) !For load increment loops
      dofType=BOUNDARY_CONDITION_DOF_FIXED
      CALL FieldVariable_ParameterSetEnsureCreated(fieldVariable,FIELD_BOUNDARY_CONDITIONS_SET_TYPE,err,error,*999)
      boundaryConditionsVariable%parameterSetRequired(FIELD_BOUNDARY_CONDITIONS_SET_TYPE)=.TRUE.
    CASE(BOUNDARY_CONDITION_PRESSURE)
      ! Pressure boundary conditions leave the RHS dof as free, as the Neumann terms
      ! are calculated in finite elasticity routines when calculating the element residual
      dofType=BOUNDARY_CONDITION_DOF_FREE
      CALL FieldVariable_ParameterSetEnsureCreated(fieldVariable,FIELD_PRESSURE_VALUES_SET_TYPE,err,error,*999)
      boundaryConditionsVariable%parameterSetRequired(FIELD_PRESSURE_VALUES_SET_TYPE)=.TRUE.
    CASE(BOUNDARY_CONDITION_PRESSURE_INCREMENTED)
      dofType=BOUNDARY_CONDITION_DOF_FREE
      CALL FieldVariable_ParameterSetEnsureCreated(fieldVariable,FIELD_PRESSURE_VALUES_SET_TYPE,err,error,*999)
      boundaryConditionsVariable%parameterSetRequired(FIELD_PRESSURE_VALUES_SET_TYPE)=.TRUE.
      CALL FieldVariable_ParameterSetEnsureCreated(fieldVariable,FIELD_PREVIOUS_PRESSURE_SET_TYPE,err,error,*999)
      boundaryConditionsVariable%parameterSetRequired(FIELD_PREVIOUS_PRESSURE_SET_TYPE)=.TRUE.
    CASE(BOUNDARY_CONDITION_CORRECTION_MASS_INCREASE)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
    CASE(BOUNDARY_CONDITION_IMPERMEABLE_WALL)
      dofType=BOUNDARY_CONDITION_DOF_FREE
      CALL FieldVariable_ParameterSetEnsureCreated(fieldVariable,FIELD_IMPERMEABLE_FLAG_VALUES_SET_TYPE,err,error,*999)
      boundaryConditionsVariable%parameterSetRequired(FIELD_IMPERMEABLE_FLAG_VALUES_SET_TYPE)=.TRUE.
    CASE(BOUNDARY_CONDITION_NEUMANN_POINT,BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
      CALL FieldVariable_ParameterSetEnsureCreated(fieldVariable,FIELD_BOUNDARY_CONDITIONS_SET_TYPE,err,error,*999)
      boundaryConditionsVariable%parameterSetRequired(FIELD_BOUNDARY_CONDITIONS_SET_TYPE)=.TRUE.
      CALL FieldVariable_ParameterSetEnsureCreated(fieldVariable,FIELD_INTEGRATED_NEUMANN_SET_TYPE,err,error,*999)
      boundaryConditionsVariable%parameterSetRequired(FIELD_INTEGRATED_NEUMANN_SET_TYPE)=.TRUE.
    CASE(BOUNDARY_CONDITION_NEUMANN_INTEGRATED,BOUNDARY_CONDITION_NEUMANN_INTEGRATED_ONLY)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
    CASE(BOUNDARY_CONDITION_FIXED_FITTED,BOUNDARY_CONDITION_FIXED_NONREFLECTING,BOUNDARY_CONDITION_FIXED_CELLML, &
      & BOUNDARY_CONDITION_FIXED_STREE)
      dofType=BOUNDARY_CONDITION_DOF_FIXED
    CASE DEFAULT
      localError="The specified boundary condition type for DOF number "//TRIM(NumberToVString(globalDof,"*",err,error))// &
        & " of "//TRIM(NumberToVString(condition,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    !We have a valid boundary condition type
    !Update condition type counts
    previousCondition=boundaryConditionsVariable%conditionTypes(globalDof)
    IF(previousCondition/=condition) THEN
      ! dofCounts array doesn't include a count for BOUNDARY_CONDITION_NONE, which equals zero
      IF(previousCondition/=BOUNDARY_CONDITION_NONE) THEN
        boundaryConditionsVariable%dofCounts(previousCondition)= &
          & boundaryConditionsVariable%dofCounts(previousCondition)-1
      ENDIF
      IF(condition/=BOUNDARY_CONDITION_NONE) THEN
        boundaryConditionsVariable%dofCounts(condition)= &
          & boundaryConditionsVariable%dofCounts(condition)+1
      ENDIF
    ENDIF
    !Update Dirichlet DOF count
    previousDof=boundaryConditionsVariable%DOFTypes(globalDof)
    IF(dofType==BOUNDARY_CONDITION_DOF_FIXED.AND.previousDof/=BOUNDARY_CONDITION_DOF_FIXED) THEN
      boundaryConditionsVariable%numberOfDirichletConditions= &
        & boundaryConditionsVariable%numberOfDirichletConditions+1
    ELSE IF(dofType/=BOUNDARY_CONDITION_DOF_FIXED.AND.previousDof==BOUNDARY_CONDITION_DOF_FIXED) THEN
      boundaryConditionsVariable%numberOfDirichletConditions= &
        & boundaryConditionsVariable%numberOfDirichletConditions-1
    ENDIF

    !Set the boundary condition type and DOF type
    boundaryConditionsVariable%conditionTypes(globalDof)=condition
    boundaryConditionsVariable%DOFTypes(globalDof)=dofType
    
    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Boundary condition being set:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Global DOF = ", globalDOF,err,error,*999)
      CALL FieldVariable_VariableTypeGet(fieldVariable,variableType,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Variable Type = ",variableType,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  New Condition = ",condition,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  DOF type = ",dofType,err,error,*999)
    ENDIF
    
    EXITS("BoundaryConditionsVariable_ConditionTypeSet")
    RETURN
999 ERRORSEXITS("BoundaryConditionsVariable_ConditionTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsVariable_ConditionTypeSet

  !
  !================================================================================================================================
  !

  !>Finalise the Dirichlet boundary conditions for a boundary conditions variable and deallocate all memory.
  SUBROUTINE BoundaryConditionsVariable_DirichletFinalise(dirichletBoundaryConditions,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsDirichletType), POINTER :: dirichletBoundaryConditions !<A pointer to the Dirichlet boundary conditions  to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsMatrixIdx,equationsSetIdx

    ENTERS("BoundaryConditionsVariable_DirichletFinalise",err,error,*999)

    IF(ASSOCIATED(dirichletBoundaryConditions)) THEN
      IF(ALLOCATED(dirichletBoundaryConditions%dirichletDOFIndices)) DEALLOCATE(dirichletBoundaryConditions%dirichletDOFIndices)
      IF(ALLOCATED(dirichletBoundaryConditions%linearSparsityIndices)) THEN
        DO equationsMatrixIdx=1,SIZE(dirichletBoundaryConditions%linearSparsityIndices,2)
          DO equationsSetIdx=1,SIZE(dirichletBoundaryConditions%linearSparsityIndices,1)
            CALL BoundaryConditionsSparsityIndices_Finalise(dirichletBoundaryConditions% &
              & linearSparsityIndices(equationsSetIdx,equationsMatrixIdx)%ptr,err,error,*999)
          ENDDO !equationsSetIdx
        ENDDO !equationsMatrixIdx
        DEALLOCATE(dirichletBoundaryConditions%linearSparsityIndices)
      ENDIF
      IF(ALLOCATED(dirichletBoundaryConditions%dynamicSparsityIndices)) THEN
        DO equationsMatrixIdx=1,SIZE(dirichletBoundaryConditions%dynamicSparsityIndices,2)
          DO equationsSetIdx=1,SIZE(dirichletBoundaryConditions%dynamicSparsityIndices,1)
            CALL BoundaryConditionsSparsityIndices_Finalise(dirichletBoundaryConditions% &
              & dynamicSparsityIndices(equationsSetIdx,equationsMatrixIdx)%ptr,err,error,*999)
          ENDDO !equationsSetIdx
        ENDDO !equationsMatrixIdx
        DEALLOCATE(dirichletBoundaryConditions%dynamicSparsityIndices)
      ENDIF
      DEALLOCATE(dirichletBoundaryConditions)
    ENDIF
    
    EXITS("BoundaryConditionsVariable_DirichletFinalise")
    RETURN
999 ERRORSEXITS("BoundaryConditions_DirichletFinalise",err,error)
    RETURN 1

  END SUBROUTINE BoundaryConditionsVariable_DirichletFinalise

  !
  !================================================================================================================================
  !

  !>Initialise dirichlet boundary conditions for a boundary conditions.
  SUBROUTINE BoundaryConditionsVariable_DirichletInitialise(boundaryConditionsVariable,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to initialise a boundary conditions dirichlet type for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,equationsSetIdx,matrixIdx,maxNumberOfLinearMatrices,maxNumberOfDynamicMatrices, &
      & numberOfDirichletConditions,numberOfEquationsSets,numberOfLinearMatrices,numberOfDynamicMatrices
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions
    TYPE(BoundaryConditionsDirichletType), POINTER :: boundaryConditionsDirichlet
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("BoundaryConditionsVariable_DirichletInitialise",err,error,*999)

    IF(.NOT.ASSOCIATED(boundaryConditionsVariable)) &
      & CALL FlagError("Boundary conditions variable is not associated.",err,error,*999)
    IF(ASSOCIATED(boundaryConditionsVariable%dirichletBoundaryConditions)) &
      CALL FlagError("Dirichlet boundary conditions are already associated for this boundary conditions variable." &
      & ,err,error,*999)
      
    ALLOCATE(boundaryConditionsVariable%dirichletBoundaryConditions,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate Dirichlet Boundary Conditions",err,error,*999)
    boundaryConditionsDirichlet=>boundaryConditionsVariable%dirichletBoundaryConditions
    numberOfDirichletConditions=boundaryConditionsVariable%numberOfDirichletConditions
    ALLOCATE(boundaryConditionsDirichlet%dirichletDOFIndices(numberOfDirichletConditions),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate Dirichlet DOF indices array",err,error,*999)

    NULLIFY(boundaryConditions)
    CALL BoundaryConditionsVariable_BoundaryConditionsGet(boundaryConditionsVariable,boundaryConditions,err,error,*999)
    NULLIFY(solverEquations)
    CALL BoundaryConditions_SolverEquationsGet(boundaryConditions,solverEquations,err,error,*999)
    NULLIFY(solverMapping)
    CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
    
    maxNumberOfLinearMatrices=0
    maxNumberOfDynamicMatrices=0
    CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
    DO equationsSetIdx=1,numberOfEquationsSets
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
      NULLIFY(equations)
      CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
      NULLIFY(vectorEquations)
      CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
      NULLIFY(vectorMapping)
      CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
      NULLIFY(linearMapping)
      CALL EquationsMappingVector_LinearMappingExists(vectorMapping,linearMapping,err,error,*999)
      IF(ASSOCIATED(linearMapping)) THEN
        CALL EquationsMappingLinear_NumberOfLinearMatricesGet(linearMapping,numberOfLinearMatrices,err,error,*999)
        IF(numberOfLinearMatrices>maxNumberOfLinearMatrices) maxNumberOfLinearMatrices=numberOfLinearMatrices
      ENDIF
      NULLIFY(dynamicMapping)
      CALL EquationsMappingVector_DynamicMappingExists(vectorMapping,dynamicMapping,err,error,*999)
      IF(ASSOCIATED(dynamicMapping)) THEN
        CALL EquationsMappingDynamic_NumberOfDynamicMatricesGet(dynamicMapping,numberOfDynamicMatrices,err,error,*999)
        IF(numberOfDynamicMatrices>maxNumberOfDynamicMatrices) maxNumberOfDynamicMatrices=numberOfDynamicMatrices
      ENDIF
    ENDDO !equationsSetIdx
    ALLOCATE(boundaryConditionsDirichlet%linearSparsityIndices(numberOfEquationsSets, maxNumberOfLinearMatrices),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate Dirichlet linear sparsity indices array",err,error,*999)
    ALLOCATE(boundaryConditionsDirichlet%dynamicSparsityIndices(numberOfEquationsSets,maxNumberOfDynamicMatrices),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate Dirichlet dynamic sparsity indices array",err,error,*999)
    DO equationsSetIdx=1,numberOfEquationsSets
      DO matrixIdx=1,maxNumberOfLinearMatrices
        NULLIFY(boundaryConditionsDirichlet%linearSparsityIndices(equationsSetIdx,matrixIdx)%ptr)
      ENDDO !matrixIdx
      DO matrixIdx=1,maxNumberOfDynamicMatrices
        NULLIFY(boundaryConditionsDirichlet%dynamicSparsityIndices(equationsSetIdx,matrixIdx)%ptr)
      ENDDO !matrixIdx
    ENDDO !equationsSetIdx

    EXITS("BoundaryConditionsVariable_DirichletInitialise")
    RETURN
999 CALL BoundaryConditionsVariable_DirichletFinalise(boundaryConditionsVariable%dirichletBoundaryConditions, &
      & dummyErr,dummyError,*998)
998 ERRORSEXITS("BoundaryConditionsVariable_DirichletInitialise",err,error)
    RETURN 1

  END SUBROUTINE BoundaryConditionsVariable_DirichletInitialise

  !
  !================================================================================================================================
  !

  !>Finish the creation of the dof constraints for a boundary conditions variable
  SUBROUTINE BoundaryConditionsVariable_DOFConstraintsCreateFinish(boundaryConditionsVariable,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable !<A pointer to boundary conditions variable to finish the dof constraints for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: constraintIdx,dofIdx,globalDOF,globalDOF2,localDOF,localDOF2,numberOfCoupledDofs,numberOfGlobalDOFs, &
      & numberOfGroupNodes,otherDOFDomain,thisDOFDomain
    INTEGER(INTG), ALLOCATABLE :: newCoupledGlobalDofs(:),newCoupledLocalDofs(:)
    REAL(DP), ALLOCATABLE :: newCoefficients(:)
    TYPE(BoundaryConditionsDofConstraintsType), POINTER :: dofConstraints
    TYPE(BoundaryConditionsDofConstraintType), POINTER :: dofConstraint
    TYPE(BoundaryConditionsCoupledDofsType), POINTER :: dofCoupling
    TYPE(DomainMappingType), POINTER :: variableDomainMapping
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(VARYING_STRING) :: localError
    TYPE(WorkGroupType), POINTER :: workGroup

    ENTERS("BoundaryConditionsVariable_DOFConstraintsCreateFinish",err,error,*998)

    !We have a list of DOF constraints, which give the values for a field variable
    !DOF as a linear combination of other DOFs.
    !In order to be able to construct the solver matrices in the solver mapping routines,
    !we create a set of couplings, where a coupling is a set of field variable DOFs
    !mapped to a single solver row or column.

    IF(.NOT.ASSOCIATED(boundaryConditionsVariable)) &
      & CALL FlagError("Boundary conditions variable is not associated.",err,error,*999)
    NULLIFY(fieldVariable)
    CALL BoundaryConditionsVariable_FieldVariableGet(boundaryConditionsVariable,fieldVariable,err,error,*999)
    NULLIFY(dofConstraints)
    CALL BoundaryConditionsVariable_DOFConstraintsGet(boundaryConditionsVariable,dofConstraints,err,error,*999)

    NULLIFY(variableDomainMapping)
    CALL FieldVariable_DomainMappingGet(fieldVariable,variableDomainMapping,err,error,*999)
    NULLIFY(workGroup)
    CALL DomainMapping_WorkGroupGet(variableDomainMapping,workGroup,err,error,*999)
    CALL WorkGroup_NumberOfGroupNodesGet(workGroup,numberOfGroupNodes,err,error,*999)
    
    !Allocate an array of pointers to DOF couplings
    IF(dofConstraints%numberOfConstraints>0) THEN
      CALL FieldVariable_NumberOfGlobalDOFsGet(fieldVariable,numberOfGlobalDOFs,err,error,*999)
      ALLOCATE(dofConstraints%dofCouplings(numberOfGlobalDofs),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate dof constraints dof couplings array.",err,error,*998)
      dofConstraints%numberOfDofs=numberOfGlobalDofs
      DO dofIdx=1,numberOfGlobalDofs
        NULLIFY(dofConstraints%dofCouplings(dofIdx)%ptr)
      ENDDO !dofIdx
    ENDIF

    !Loop over all constraints
    DO constraintIdx=1,dofConstraints%numberOfConstraints
      NULLIFY(dofConstraint)
      CALL BoundaryConditionsVariableDOFConstraints_ConstraintGet(dofConstraints,constraintIdx,dofConstraint,err,error,*999)
      
      globalDOF=dofConstraint%globalDof
      CALL DomainMapping_LocalNumberFromGlobalGet(variableDomainMapping,globalDOF,1,localDOF,err,error,*999)
      CALL DomainMapping_DomainNumberFromGlobalGet(variableDomainMapping,globalDOF,1,thisDOFDomain,err,error,*999)
      
      !Check that the constrained DOFs are still set to be constrained, as
      !subsequently setting a boundary condition would change the DOF type but
      !not update the DOF constraints structure.
      IF(boundaryConditionsVariable%DOFTypes(globalDof)/=BOUNDARY_CONDITION_DOF_CONSTRAINED) THEN
        localError="Global DOF number "//TRIM(NumberToVstring(globalDOF,"*",err,error))// &
          & " is part of a linear constraint but the DOF type has been changed"// &
          & " by applying a boundary condition."
        CALL FlagError(localError,err,error,*999)
      ENDIF

      DO dofIdx=1,dofConstraint%numberOfDofs
        globalDOF2=dofConstraint%dofs(dofIdx)
        CALL DomainMapping_LocalNumberFromGlobalGet(variableDomainMapping,globalDOF2,1,localDOF2,err,error,*999)
        !Check a Dirichlet conditions hasn't also been set on this DOF
        IF(boundaryConditionsVariable%DOFTypes(globalDOF2)/=BOUNDARY_CONDITION_DOF_FREE) THEN
          localError="A Dirichlet boundary condition has been set on DOF number "// &
            & TRIM(NumberToVstring(globalDOF2,"*",err,error))// &
            & " which is part of a linear constraint."
          CALL FlagError(localError,err,error,*999)
        ENDIF

        !Check we don't have DOF constraints that are split over domains
        !\todo Implement support for DOF constraints that are split over domains
        IF(numberOfGroupNodes>1) THEN
          CALL DomainMapping_DomainNumberFromGlobalGet(variableDomainMapping,globalDOF2,1,otherDOFDomain,err,error,*999)
          IF(thisDOFDomain/=otherDOFDomain) THEN
            CALL FlagError("An equal DOF constraint is split over multiple domains, "// &
              & "support for this has not yet been implemented.",err,error,*999)
          ENDIF
        ENDIF

        !Add to the DOFs that are coupled with globalDOF2
        !This might be quite inefficient if there are many dofs mapped to a single row/column
        !due to the reallocation at each step
        IF(ASSOCIATED(dofConstraints%dofCouplings(globalDOF2)%ptr)) THEN
          numberOfCoupledDofs=dofConstraints%dofCouplings(globalDOF2)%ptr%numberOfDofs
          ALLOCATE(newCoupledGlobalDofs(numberOfCoupledDofs+1),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate new DOF coupling global DOFs.",err,error,*999)
          ALLOCATE(newCoupledLocalDofs(numberOfCoupledDofs+1),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate new DOF coupling local DOFs.",err,error,*999)
          ALLOCATE(newCoefficients(numberOfCoupledDofs+1),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate new DOF coupling values.",err,error,*999)
          newCoupledGlobalDofs(1:numberOfCoupledDofs)=dofCoupling%globalDofs(1:numberOfCoupledDofs)
          newCoupledLocalDofs(1:numberOfCoupledDofs)=dofCoupling%localDofs(1:numberOfCoupledDofs)
          newCoefficients(1:numberOfCoupledDofs)=dofCoupling%coefficients(1:numberOfCoupledDofs)
        ELSE
          !Set up a a new dofCoupling and set globalDof2 as the first DOF
          ALLOCATE(dofConstraints%dofCouplings(globalDof2)%ptr,STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate new DOF coupling type.",err,error,*999)
          ALLOCATE(newCoupledGlobalDofs(2),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate new DOF coupling global DOFs.",err,error,*999)
          ALLOCATE(newCoupledLocalDofs(2),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate new DOF coupling local DOFs.",err,error,*999)
          ALLOCATE(newCoefficients(2),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate new DOF coupling values.",err,error,*999)
          newCoupledGlobalDofs(1)=globalDOF2
          newCoupledLocalDofs(1)=localDOF2
          newCoefficients(1)=1.0_DP
          numberOfCoupledDofs=1
        ENDIF
        NULLIFY(dofCoupling)
        CALL BoundaryConditionsVariableDOFConstraints_DOFCouplingGet(dofConstraints,globalDOF2,dofCoupling,err,error,*999)
        newCoupledGlobalDofs(numberOfCoupledDofs+1)=globalDOF
        newCoupledLocalDofs(numberOfCoupledDofs+1)=localDOF
        newCoefficients(numberOfCoupledDofs+1)=dofConstraint%coefficients(dofIdx)
        CALL MOVE_ALLOC(newCoupledGlobalDofs,dofCoupling%globalDofs)
        CALL MOVE_ALLOC(newCoupledLocalDofs,dofCoupling%localDofs)
        CALL MOVE_ALLOC(newCoefficients,dofCoupling%coefficients)
        dofCoupling%numberOfDofs=numberOfCoupledDofs+1
      ENDDO !dofIdx
    ENDDO !constraintIdx

    EXITS("BoundaryConditionsVariable_DOFConstraintsCreateFinish")
    RETURN
999 IF(ALLOCATED(newCoupledGlobalDofs)) DEALLOCATE(newCoupledGlobalDofs)
    IF(ALLOCATED(newCoupledLocalDofs)) DEALLOCATE(newCoupledLocalDofs)
    IF(ALLOCATED(newCoefficients)) DEALLOCATE(newCoefficients)
    CALL BoundaryConditionsDOFConstraints_Finalise(dofConstraints,err,error,*998)
998 ERRORS("BoundaryConditionsVariable_DOFConstraintsCreateFinish",err,error)
    EXITS("BoundaryConditionsVariable_DOFConstraintsCreateFinish")
    RETURN 1

  END SUBROUTINE BoundaryConditionsVariable_DOFConstraintsCreateFinish

  !
  !================================================================================================================================
  !

  !>Sets the boundary condition DOF type for a boundary condition variable.
  SUBROUTINE BoundaryConditionsVariable_DOFTypeSet(boundaryConditionsVariable,globalDOF,dofType,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to set the boundary condition DOF typefor
    INTEGER(INTG), INTENT(IN) :: globalDOF !<The globalDof to set the boundary condition at
    INTEGER(INTG), INTENT(IN) :: dofType !<The boundary condition DOF type to set \see BoundaryConditionsRoutines_DOFTypes,BoundaryConditionsRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: bcDOFType,currentDOFType,numberOfGlobalDOFs,previousCondition,previousDOF,variableType
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("BoundaryConditionsVariable_DOFTypeSet",err,error,*999)

    NULLIFY(fieldVariable)
    CALL BoundaryConditionsVariable_FieldVariableGet(boundaryConditionsVariable,fieldVariable,err,error,*999)
    CALL FieldVariable_NumberOfGlobalDOFsGet(fieldVariable,numberOfGlobalDOFs,err,error,*999)
    IF(globalDOF<1.OR.globalDOF>numberOfGlobalDOFs) THEN
      localError="The specified global DOF number of "//TRIM(NumberToVString(globalDOF,"*",err,error))// &
        & " is invalid. The global number should be >= 1 and <= "//TRIM(NumberToVString(numberOfGlobalDOFs,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    CALL BoundaryConditionsVariable_DOFTypeGet(boundaryConditionsVariable,globalDOF,currentDOFType,err,error,*999)

    IF(dofType/=currentDOFType) THEN
      SELECT CASE(dofType)
      CASE(BOUNDARY_CONDITION_DOF_FREE)
        bcDOFType=BOUNDARY_CONDITION_DOF_FREE
      CASE(BOUNDARY_CONDITION_DOF_FIXED)
        bcDOFType=BOUNDARY_CONDITION_DOF_FIXED
      CASE(BOUNDARY_CONDITION_DOF_INCREMENTED)
        bcDOFType=BOUNDARY_CONDITION_DOF_INCREMENTED
      CASE(BOUNDARY_CONDITION_DOF_MIXED)
        bcDOFType=BOUNDARY_CONDITION_DOF_MIXED
      CASE(BOUNDARY_CONDITION_DOF_CONSTRAINED)
        bcDOFType=BOUNDARY_CONDITION_DOF_CONSTRAINED
      CASE DEFAULT
        localError="The specified boundary condition DOF type of "//TRIM(NumberToVString(dofType,"*",err,error))// &
          & " for DOF number "//TRIM(NumberToVString(globalDOF,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT

      !We have a valid boundary condition type
      boundaryConditionsVariable%DOFTypes(globalDOF)=bcDOFType
      
    ENDIF
  
    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Boundary condition being set:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Global DOF    = ", globalDOF,err,error,*999)
      CALL FieldVariable_VariableTypeGet(fieldVariable,variableType,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Variable type = ",variableType,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Old DOF type  = ",currentDOFType,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  New DOF type  = ",dofType,err,error,*999)
    ENDIF
    
    EXITS("BoundaryConditionsVariable_DOFTypeSet")
    RETURN
999 ERRORSEXITS("BoundaryConditionsVariable_DOFTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsVariable_DOFTypeSet

  !
  !================================================================================================================================
  !

  !Finalises a boundary conditions variable and deallocates all memory.
  SUBROUTINE BoundaryConditionsVariable_Finalise(boundaryConditionsVariable,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables

    ENTERS("BoundaryConditionsVariable_Finalise",err,error,*999)

    IF(ASSOCIATED(boundaryConditionsVariable)) THEN
      IF(ALLOCATED(boundaryConditionsVariable%dofTypes)) DEALLOCATE(boundaryConditionsVariable%dofTypes)
      IF(ALLOCATED(boundaryConditionsVariable%conditionTypes)) DEALLOCATE(boundaryConditionsVariable%conditionTypes)
      CALL BoundaryConditionsVariable_DirichletFinalise(boundaryConditionsVariable%dirichletBoundaryConditions,err,error,*999)
      CALL BoundaryConditionsVariable_NeumannFinalise(boundaryConditionsVariable%neumannBoundaryConditions,err,error,*999)
      CALL BoundaryConditionsVariable_PressureIncrementedFinalise(boundaryConditionsVariable% &
        & pressureIncrementedBoundaryConditions,err,error,*999)
      IF(ALLOCATED(boundaryConditionsVariable%dofCounts)) DEALLOCATE(boundaryConditionsVariable%dofCounts)
      IF(ALLOCATED(boundaryConditionsVariable%parameterSetRequired)) DEALLOCATE(boundaryConditionsVariable%parameterSetRequired)
      CALL BoundaryConditionsDOFConstraints_Finalise(boundaryConditionsVariable%dofConstraints,err,error,*999)
      DEALLOCATE(boundaryConditionsVariable)
    ENDIF
    
    EXITS("BoundaryConditionsVariable_Finalise")
    RETURN
999 ERRORSEXITS("BoundaryConditionsVariable_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsVariable_Finalise

  !
  !================================================================================================================================
  !

  !Allocates and initialises a boundary conditions variable.
  SUBROUTINE BoundaryConditionsVariable_Initialise(boundaryConditionsVariable,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("BoundaryConditionsVariable_Initialise",err,error,*998)

    IF(ASSOCIATED(boundaryConditionsVariable)) &
      & CALL FlagError("Boundary conditions variable is already associated.",err,error,*998)

    ALLOCATE(boundaryConditionsVariable,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocated boundary conditions variable.",err,error,*999)
    NULLIFY(boundaryConditionsVariable%boundaryConditions)
    boundaryConditionsVariable%variableType=0
    NULLIFY(boundaryConditionsVariable%variable)
    boundaryConditionsVariable%numberOfDirichletConditions=0
    NULLIFY(boundaryConditionsVariable%dirichletBoundaryConditions)
    NULLIFY(boundaryConditionsVariable%neumannBoundaryConditions)
    NULLIFY(boundaryConditionsVariable%pressureIncrementedBoundaryConditions)
    NULLIFY(boundaryConditionsVariable%dofConstraints)
    
    EXITS("BoundaryConditionsVariable_Initialise")
    RETURN
999 CALL BoundaryConditionsVariable_Finalise(boundaryConditionsVariable,dummyErr,dummyError,*998)
998 ERRORSEXITS("BoundaryConditionsVariable_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsVariable_Initialise

  !
  !================================================================================================================================
  !

  !>Finalise the Neumann condition information and deallocates all memory
  SUBROUTINE BoundaryConditionsVariable_NeumannFinalise(neumannBoundaryConditions,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsNeumannType), POINTER :: neumannBoundaryConditions !<A pointer to the Neumann boundary conditions to finalise 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(BoundaryConditionsNeumannType), POINTER :: boundaryConditionsNeumann

    ENTERS("BoundaryConditionsVariable_NeumannFinalise",err,error,*999)

    IF(ASSOCIATED(boundaryConditionsNeumann)) THEN
      IF(ALLOCATED(boundaryConditionsNeumann%setDOFs)) DEALLOCATE(boundaryConditionsNeumann%setDOFs)
      CALL BoundaryConditionsVariable_NeumannMatricesFinalise(boundaryConditionsNeumann,err,error,*999)
      DEALLOCATE(boundaryConditionsNeumann)
    ENDIF

    EXITS("BoundaryConditionsVariable_NeumannFinalise")
    RETURN
999 ERRORSEXITS("BoundaryConditionsVariable_NeumannFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsVariable_NeumannFinalise

  !
  !================================================================================================================================
  !

  !>Allocates and initialise the Neumann boundary conditions for a boundary conditions variable
  SUBROUTINE BoundaryConditionsVariable_NeumannInitialise(boundaryConditionsVariable,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to initialise Neumann conditions for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,numberOfValues,numberOfLocalDOFs
    TYPE(BoundaryConditionsNeumannType), POINTER :: boundaryConditionsNeumann
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("BoundaryConditionsVariable_NeumannInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(boundaryConditionsVariable)) &
      & CALL FlagError("Boundary conditions variable is not associated.",err,error,*998)
    IF(ASSOCIATED(boundaryConditionsVariable%neumannBoundaryConditions)) &
      & CALL FlagError("Neumann boundary conditions are already associated for the boundary conditions variable.",err,error,*998)
    
    numberOfValues=boundaryConditionsVariable%dofCounts(BOUNDARY_CONDITION_NEUMANN_POINT)+ &
      & boundaryConditionsVariable%dofCounts(BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)
    
    ALLOCATE(boundaryConditionsVariable%neumannBoundaryConditions,stat=err)
    IF(err/=0) CALL FlagError("Could not allocate Neumann boundary conditions.",err,error,*999)
    NULLIFY(boundaryConditionsNeumann)
    CALL BoundaryConditionsVariable_NeumannConditionsGet(boundaryConditionsVariable,boundaryConditionsNeumann,err,error,*999)
    NULLIFY(boundaryConditionsNeumann%integrationMatrix)
    NULLIFY(boundaryConditionsNeumann%pointValues)
    NULLIFY(boundaryConditionsNeumann%pointDofMapping)

    NULLIFY(fieldVariable)
    CALL BoundaryConditionsVariable_FieldVariableGet(boundaryConditionsVariable,fieldVariable,err,error,*999)
    CALL FieldVariable_NumberOfDOFsGet(fieldVariable,numberOfLocalDOFs,err,error,*999)
    ALLOCATE(boundaryConditionsNeumann%setDOFs(numberOfLocalDOFs),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate Neumann set DOFs.",err,error,*999)
    boundaryConditionsNeumann%setDOFs=0

    EXITS("BoundaryConditionsVariable_NeumannInitialise")
    RETURN
999 CALL BoundaryConditionsVariable_NeumannFinalise(boundaryConditionsVariable%neumannBoundaryConditions,dummyErr,dummyError,*998)
998 ERRORSEXITS("BoundaryConditionsVariable_NeumannInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsVariable_NeumannInitialise

  !
  !================================================================================================================================
  !

  !>Calculates integrated Neumann condition values from point values for a boundary conditions variable and
  !>updates the FIELD_INTEGRATED_NEUMANN_SET_TYPE parameter set for the field variable.
  SUBROUTINE BoundaryConditionsVariable_NeumannIntegrate(rhsBoundaryConditions,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsVariableType), POINTER, INTENT(IN) :: rhsBoundaryConditions !<The boundary conditions for the right hand side field variable
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: columnElementParameterIdx,componentNumber,derivIdx,derivativeNumber,dofType,dofTypeIdx,domainNumber, &
      & faceIdx,faceNumber,gaussIdx,globalDOF,interpolationType,lineIdx,localDOF,lineNumber,myGroupNodeNumber,neumannDOFIdx, &
      & neumannGlobalDOF,neumannLocalDerivNumber,neumannLocalNodeNumber,neumannLocalDOF,neumannNodeNumber,nodeIdx,nodeNumber, &
      & numberOfDimensions,numberOfGauss,numberOfNeumann,numberOfNodes,numberOfNodeDerivatives,numberOfNodeFaces, &
      & numberOfNodeLines,rhsScalingType,rowElementParameterIdx,versionNumber
    REAL(DP) :: columnPhi,gaussWeight,integratedValue,jacobian,rowPhi
    LOGICAL :: boundaryFace,boundaryLine,calculateFaces,calculateLines,dependentGeometry
    TYPE(BoundaryConditionsNeumannType), POINTER :: neumannConditions
    TYPE(BasisType), POINTER :: basis
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DistributedVectorType), POINTER :: integratedValues
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainFacesType), POINTER :: domainFaces
    TYPE(DomainLinesType), POINTER :: domainLines
    TYPE(DomainMappingType), POINTER :: rhsDomainMapping,rowsMapping
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(FieldType), POINTER :: geometricField,rhsField
    TYPE(FieldVariableType), POINTER :: geometricVariable,rhsVariable
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: interpolatedPointMetrics
    TYPE(FieldInterpolatedPointType), POINTER :: interpolatedPoint
    TYPE(FieldInterpolationParametersType), POINTER :: interpolationParameters,scalingParameters
    TYPE(QuadratureSchemeType), POINTER :: quadratureScheme
    TYPE(VARYING_STRING) :: localError
    TYPE(WorkGroupType), POINTER :: workGroup

    ENTERS("BoundaryConditionsVariable_NeumannIntegrate",err,error,*999)

    !Check that Neumann conditions are associated, otherwise do nothing
    NULLIFY(neumannConditions)
    CALL BoundaryConditionsVariable_NeumannConditionsExists(rhsBoundaryConditions,neumannConditions,err,error,*999)
    IF(ASSOCIATED(neumannConditions)) THEN
      NULLIFY(rhsVariable)
      CALL BoundaryConditionsVariable_FieldVariableGet(rhsBoundaryConditions,rhsVariable,err,error,*999)
!!TODO: We need to access the LHS mapping for the equations set. Just use variable for now. 
      !For rows we can re-use the RHS variable row mapping
      NULLIFY(rowsMapping)
      CALL FieldVariable_DomainMappingGet(rhsVariable,rowsMapping,err,error,*999)
      NULLIFY(rhsField)
      CALL FieldVariable_FieldGet(rhsVariable,rhsField,err,error,*999)
      CALL Field_ScalingTypeGet(rhsField,rhsScalingType,err,error,*999)
      NULLIFY(geometricField)
      CALL Field_GeometricGeneralFieldGet(rhsField,geometricField,dependentGeometry,err,error,*999)
      NULLIFY(geometricVariable)
      CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)

      CALL DistributedMatrix_AllValuesSet(neumannConditions%integrationMatrix,0.0_DP,err,error,*999)

      numberOfNeumann=rhsBoundaryConditions%dofCounts(BOUNDARY_CONDITION_NEUMANN_POINT) + &
        & rhsBoundaryConditions%dofCounts(BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)

      NULLIFY(decomposition)
      CALL Field_DecompositionGet(rhsField,decomposition,err,error,*999)
      NULLIFY(workGroup)
      CALL Decomposition_WorkGroupGet(decomposition,workGroup,err,error,*999)
      CALL WorkGroup_GroupNodeNumberGet(workGroup,myGroupNodeNumber,err,error,*999)

      ! Initialise field interpolation parameters for the geometric field, which are required for the
      ! face/line Jacobian and scale factors
      NULLIFY(interpolationParameters)
      CALL FieldVariable_InterpolationParameterInitialise(geometricVariable,interpolationParameters,err,error,*999)
      NULLIFY(scalingParameters)
      CALL FieldVariable_InterpolationParameterInitialise(rhsVariable,scalingParameters,err,error,*999)
      NULLIFY(interpolatedPoint)
      CALL Field_InterpolatedPointInitialise(interpolationParameters,interpolatedPoint,err,error,*999)
      NULLIFY(interpolatedPointMetrics)
      CALL Field_InterpolatedPointMetricsInitialise(interpolatedPoint,interpolatedPointMetrics,err,error,*999)

      ! Loop over all Neumann point DOFs, finding the boundary lines or faces they are on
      ! and integrating over them
      NULLIFY(rhsDomainMapping)
      CALL FieldVariable_DomainMappingGet(rhsVariable,rhsDomainMapping,err,error,*999)
      DO neumannDOFIdx=1,numberOfNeumann
        neumannGlobalDOF=neumannConditions%setDofs(neumannDOFIdx)
        CALL DomainMapping_DomainNumberFromGlobalGet(rhsDomainMapping,neumannGlobalDOF,1,domainNumber,err,error,*999)
        IF(domainNumber==myGroupNodeNumber) THEN
          CALL DomainMapping_LocalNumberFromGlobalGet(rhsDomainMapping,neumannGlobalDOF,1,neumannLocalDOF,err,error,*999)
          !Get Neumann DOF component and topology for that component
          CALL FieldVariable_DOFTypeGet(rhsVariable,neumannLocalDOF,dofType,dofTypeIdx,err,error,*999)
          SELECT CASE(dofType)
          CASE(FIELD_CONSTANT_INTERPOLATION)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(FIELD_NODE_BASED_INTERPOLATION)
            CALL FieldVariable_DOFParametersGetNode(rhsVariable,dofTypeIdx,versionNumber,derivativeNumber,neumannNodeNumber, &
              & componentNumber,err,error,*999)
          CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The DOF type of "//TRIM(NumberToVString(dofType,"*",err,error))// &
              & " for local DOF "//TRIM(NumberToVString(localDOF,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          NULLIFY(domain)
          CALL FieldVariable_ComponentDomainGet(rhsVariable,componentNumber,domain,err,error,*999)
          CALL Domain_NumberOfDimensionsGet(domain,numberOfDimensions,err,error,*999)
          NULLIFY(domainTopology)
          CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
          CALL FieldVariable_ComponentInterpolationGet(rhsVariable,componentNumber,interpolationType,err,error,*999)
          NULLIFY(decomposition)
          CALL Domain_DecompositionGet(domain,decomposition,err,error,*999)
          SELECT CASE(interpolationType)
          CASE(FIELD_CONSTANT_INTERPOLATION)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(FIELD_NODE_BASED_INTERPOLATION)
            NULLIFY(domainNodes)
            CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
            SELECT CASE(numberOfDimensions)
            CASE(1)
              CALL DistributedMatrix_ValuesSet(neumannConditions%integrationMatrix,neumannLocalDof,neumannDofIdx,1.0_DP, &
                & err,error,*999)
            CASE(2)
              CALL Decomposition_CalculateLinesGet(decomposition,calculateLines,err,error,*999)
              IF(.NOT.calculateLines) CALL FlagError("Decomposition does not have lines calculated.",err,error,*999)
              NULLIFY(domainLines)
              CALL DomainTopology_DomainLinesGet(domainTopology,domainLines,err,error,*999)
              CALL DomainNodes_NodeNumberOfLinesGet(domainNodes,neumannNodeNumber,numberOfNodeLines,err,error,*999)
              linesLoop: DO lineIdx=1,numberOfNodeLines
                CALL DomainNodes_NodeLineNumberGet(domainNodes,lineIdx,nodeNumber,lineNumber,err,error,*999)
                CALL DomainLines_LineBoundaryLineGet(domainLines,lineNumber,boundaryLine,err,error,*999)
               IF(.NOT.boundaryLine) CYCLE linesLoop
                CALL DomainLines_LineBasisGet(domainLines,lineNumber,basis,err,error,*999)
                CALL Basis_NumberOfLocalNodesGet(basis,numberOfNodes,err,error,*999)
                neumannLocalNodeNumber=0
                neumannLocalDerivNumber=0
                ! Check all nodes in line to find the local numbers for the Neumann DOF, and
                ! make sure we don't have an integrated_only condition set on the line
                DO nodeIdx=1,numberOfNodes
                  CALL DomainLines_LineNodeNumberGet(domainLines,nodeIdx,lineNumber,nodeNumber,err,error,*999)
                  CALL Basis_NodeNumberOfDerivativesGet(basis,nodeIdx,numberOfNodeDerivatives,err,error,*999)
                  DO derivIdx=1,numberOfNodeDerivatives
                    CALL DomainLines_DerivativeGlobalIndexGet(domainLines,derivIdx,nodeIdx,lineNumber,derivativeNumber, &
                      & err,error,*999)
                    CALL DomainLines_DerivativeVersionNumberGet(domainLines,derivIdx,nodeIdx,lineNumber,versionNumber, &
                      & err,error,*999)
                    CALL FieldVariable_LocalNodeDOFGet(rhsVariable,versionNumber,derivativeNumber,nodeNumber, &
                      & componentNumber,localDOF,err,error,*999)
                    CALL DomainMapping_LocalToGlobalGet(rhsDomainMapping,localDOF,globalDOF,err,error,*999)
                    IF(globalDOF==neumannGlobalDOF) THEN
                      neumannLocalNodeNumber=nodeIdx
                      neumannLocalDerivNumber=derivIdx
                    ELSE IF(rhsBoundaryConditions%conditionTypes(globalDof)==BOUNDARY_CONDITION_NEUMANN_INTEGRATED_ONLY) THEN
                      CYCLE linesLoop
                    ENDIF
                  ENDDO !derivIdx
                ENDDO !nodeIdx
                IF(neumannLocalNodeNumber==0) THEN
                  CALL FlagError("Could not find local Neumann node and derivative numbers in line.",err,error,*999)
                END IF

                ! Now perform actual integration
                NULLIFY(quadratureScheme)
                CALL Basis_QuadratureSchemeGet(basis,BASIS_DEFAULT_QUADRATURE_SCHEME,quadratureScheme,err,error,*999)
                CALL BasisQuadratureScheme_NumberOfGaussGet(quadratureScheme,numberOfGauss,err,error,*999)
                CALL Field_InterpolationParametersLineGet(FIELD_VALUES_SET_TYPE,lineNumber,interpolationParameters, &
                  & err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
                IF(rhsScalingType/=FIELD_NO_SCALING) &
                  & CALL Field_InterpolationParametersScaleFactorsLineGet(lineNumber,scalingParameters,err,error,*999)
 
                DO nodeIdx=1,numberOfNodes
                  CALL DomainLines_LineNodeNumberGet(domainLines,nodeIdx,lineNumber,nodeNumber,err,error,*999)
                  CALL Basis_NodeNumberOfDerivativesGet(basis,nodeIdx,numberOfNodeDerivatives,err,error,*999)
                  DO derivIdx=1,numberOfNodeDerivatives
                    CALL DomainLines_DerivativeGlobalIndexGet(domainLines,derivIdx,nodeIdx,lineNumber,derivativeNumber, &
                      & err,error,*999)
                    CALL DomainLines_DerivativeVersionNumberGet(domainLines,derivIdx,nodeIdx,lineNumber,versionNumber, &
                      & err,error,*999)
                    CALL FieldVariable_LocalNodeDOFGet(rhsVariable,versionNumber,derivativeNumber,nodeNumber, &
                      & componentNumber,localDOF,err,error,*999)
                    
                    CALL DomainMapping_LocalToGlobalGet(rowsMapping,localDOF,globalDOF,err,error,*999)
                    IF(rhsBoundaryConditions%conditionTypes(globalDof)==BOUNDARY_CONDITION_NEUMANN_POINT.OR. &
                      & rhsBoundaryConditions%conditionTypes(globalDof)==BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED) THEN
                      
                      CALL Basis_ElementParameterGet(basis,derivIdx,nodeIdx,rowElementParameterIdx,err,error,*999)
                      CALL Basis_ElementParameterGet(basis,neumannLocalDerivNumber,neumannLocalNodeNumber, &
                        & columnElementParameterIdx,err,error,*999)
                      
                      integratedValue=0.0_DP
                      
                      ! Loop over line gauss points, adding gauss weighted terms to the integral
                      DO gaussIdx=1,numberOfGauss
                        CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx,interpolatedPoint, &
                          & err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
                        CALL Field_InterpolatedPointMetricsCalculate(COORDINATE_JACOBIAN_LINE_TYPE,interpolatedPointMetrics, &
                          & err,error,*999)
                        
                        CALL BasisQuadratureScheme_GaussWeightGet(quadratureScheme,gaussIdx,gaussWeight,err,error,*999)
                        CALL FieldInterpolatedPointMetrics_JacobianGet(interpolatedPointMetrics,jacobian,err,error,*999)
                        
                        !Get basis function values at guass points
                        CALL BasisQuadratureScheme_GaussBasisFunctionGet(quadratureScheme,rowElementParameterIdx,NO_PART_DERIV, &
                          & gaussIdx,rowPhi,err,error,*999)
                        CALL BasisQuadratureScheme_GaussBasisFunctionGet(quadratureScheme,columnElementParameterIdx,NO_PART_DERIV, &
                          & gaussIdx,columnPhi,err,error,*999)
                        
                        !Add gauss point value to total line integral
                        integratedValue=integratedValue+rowPhi*columnPhi*gaussWeight*jacobian
                        
                      ENDDO !gaussIdx
                      
                      !Multiply by scale factors for dependent variable
                      IF(rhsScalingType/=FIELD_NO_SCALING) THEN
                        integratedValue=integratedValue* &
                          & scalingParameters%scaleFactors(rowElementParameterIdx,componentNumber)* &
                          & scalingParameters%scaleFactors(columnElementParameterIdx,componentNumber)
                      ENDIF
                      
                      !Add integral term to N matrix
                      CALL DistributedMatrix_ValuesAdd(neumannConditions%integrationMatrix,localDOF,neumannDOFIdx, &
                        & integratedValue,err,error,*999)
                      
                    ENDIF
                  ENDDO !derivIdx
                ENDDO !nodeIdx
              ENDDO linesLoop
            CASE(3)
              CALL Decomposition_CalculateFacesGet(decomposition,calculateFaces,err,error,*999)
              IF(.NOT.calculateFaces) CALL FlagError("Decomposition does not have faces calculated.",err,error,*999)
              NULLIFY(domainFaces)
              CALL DomainTopology_DomainFacesGet(domainTopology,domainFaces,err,error,*999)
              CALL DomainNodes_NodeNumberOfFacesGet(domainNodes,neumannNodeNumber,numberOfNodeFaces,err,error,*999)
              facesLoop: DO faceIdx=1,numberOfNodeFaces
                CALL DomainNodes_NodeFaceNumberGet(domainNodes,faceIdx,neumannNodeNumber,faceNumber,err,error,*999)
                CALL DomainFaces_FaceBoundaryFaceGet(domainFaces,faceNumber,boundaryface,err,error,*999)
                IF(.NOT.boundaryFace) CYCLE
                NULLIFY(basis)
                CALL DomainFaces_FaceBasisGet(domainFaces,faceNumber,basis,err,error,*999)
                CALL Basis_NumberOfLocalNodesGet(basis,numberOfNodes,err,error,*999)
                neumannLocalNodeNumber=0
                neumannLocalDerivNumber=0
                ! Check all nodes in the face to find the local numbers for the Neumann DOF, and
                ! make sure we don't have an integrated_only condition set on the face
                DO nodeIdx=1,numberOfNodes
                  CALL DomainFaces_FaceNodeNumberGet(domainFaces,nodeIdx,faceNumber,nodeNumber,err,error,*999)
                  CALL Basis_NodeNumberOfDerivativesGet(basis,nodeIdx,numberOfNodeDerivatives,err,error,*999)
                  DO derivIdx=1,numberOfNodeDerivatives
                    CALL DomainFaces_DerivativeGlobalIndexGet(domainFaces,derivIdx,nodeIdx,faceNumber,derivativeNumber, &
                      & err,error,*999)
                    CALL DomainFaces_DerivativeVersionNumberGet(domainFaces,derivIdx,nodeIdx,faceNumber,versionNumber, &
                      & err,error,*999)
                    CALL FieldVariable_LocalNodeDOFGet(rhsVariable,versionNumber,derivativeNumber,nodeNumber, &
                      & componentNumber,localDOF,err,error,*999)
                    CALL DomainMapping_LocalToGlobalGet(rhsDomainMapping,localDOF,globalDOF,err,error,*999)
                    IF(globalDOF==neumannGlobalDOF) THEN
                      neumannLocalNodeNumber=nodeIdx
                      neumannLocalDerivNumber=derivIdx
                    ELSE IF(rhsBoundaryConditions%conditionTypes(globalDOF)==BOUNDARY_CONDITION_NEUMANN_INTEGRATED_ONLY) THEN
                      CYCLE facesLoop
                    ENDIF
                  ENDDO !derivIdx
                ENDDO !nodeIdx
                IF(neumannLocalNodeNumber==0) &
                  & CALL FlagError("Could not find local Neumann node and derivative numbers in line.",err,error,*999)

                ! Now perform actual integration
                NULLIFY(quadratureScheme)
                CALL Basis_QuadratureSchemeGet(basis,BASIS_DEFAULT_QUADRATURE_SCHEME,quadratureScheme,err,error,*999)
                CALL BasisQuadratureScheme_NumberOfGaussGet(quadratureScheme,numberOfGauss,err,error,*999)
                CALL Field_InterpolationParametersFaceGet(FIELD_VALUES_SET_TYPE,faceNumber,interpolationParameters, &
                  & err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
                IF(rhsScalingType/=FIELD_NO_SCALING) &
                  & CALL Field_InterpolationParametersScaleFactorsFaceGet(faceNumber,scalingParameters,err,error,*999)

                DO nodeIdx=1,numberOfNodes
                  CALL DomainFaces_FaceNodeNumberGet(domainFaces,nodeIdx,faceNumber,nodeNumber,err,error,*999)
                  CALL Basis_NodeNumberOfDerivativesGet(basis,nodeIdx,numberOfNodeDerivatives,err,error,*999)
                  DO derivIdx=1,numberOfNodeDerivatives
                    CALL DomainFaces_DerivativeGlobalIndexGet(domainFaces,derivIdx,nodeIdx,faceNumber,derivativeNumber, &
                      & err,error,*999)
                    CALL DomainFaces_DerivativeVersionNumberGet(domainFaces,derivIdx,nodeIdx,faceNumber,versionNumber, &
                      & err,error,*999)
                    CALL FieldVariable_LocalNodeDOFGet(rhsVariable,versionNumber,derivativeNumber,nodeNumber, &
                      & componentNumber,localDOF,err,error,*999)
                    
                    CALL DomainMapping_LocalToGlobalGet(rowsMapping,localDOF,globalDOF,err,error,*999)
                    IF(rhsBoundaryConditions%conditionTypes(globalDof)==BOUNDARY_CONDITION_NEUMANN_POINT.OR. &
                      & rhsBoundaryConditions%conditionTypes(globalDof)==BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED) THEN
                      
                      CALL Basis_ElementParameterGet(basis,derivIdx,nodeIdx,rowElementParameterIdx,err,error,*999)
                      CALL Basis_ElementParameterGet(basis,neumannLocalDerivNumber,neumannLocalNodeNumber, &
                        & columnElementParameterIdx,err,error,*999)
                      
                      integratedValue=0.0_DP
                      ! Loop over line gauss points, adding gauss weighted terms to the integral
                      DO gaussIdx=1,numberOfGauss
                        CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx,interpolatedPoint, &
                          & err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
                        CALL Field_InterpolatedPointMetricsCalculate(COORDINATE_JACOBIAN_AREA_TYPE,interpolatedPointMetrics, &
                          & err,error,*999)
                        
                        CALL BasisQuadratureScheme_GaussWeightGet(quadratureScheme,gaussIdx,gaussWeight,err,error,*999)
                        CALL FieldInterpolatedPointMetrics_JacobianGet(interpolatedPointMetrics,jacobian,err,error,*999)
                        
                        !Get basis function values at guass points
                        CALL BasisQuadratureScheme_GaussBasisFunctionGet(quadratureScheme,rowElementParameterIdx,NO_PART_DERIV, &
                          & gaussIdx,rowPhi,err,error,*999)
                        CALL BasisQuadratureScheme_GaussBasisFunctionGet(quadratureScheme,columnElementParameterIdx,NO_PART_DERIV, &
                          & gaussIdx,columnPhi,err,error,*999)
                        
                        !Add gauss point value to total line integral
                        integratedValue=integratedValue+rowPhi*columnPhi*gaussWeight*jacobian
                        
                      ENDDO !gaussIdx
                      
                      ! Multiply by scale factors
                      IF(rhsScalingType/=FIELD_NO_SCALING) THEN
                        integratedValue=integratedValue* &
                          & scalingParameters%scaleFactors(rowElementParameterIdx,componentNumber)* &
                          & scalingParameters%scaleFactors(columnElementParameterIdx,componentNumber)
                      ENDIF
                      
                      ! Add integral term to N matrix
                      CALL DistributedMatrix_ValuesAdd(neumannConditions%integrationMatrix,localDOF,neumannDOFIdx, &
                        & integratedValue,err,error,*999)

                    ENDIF
                  ENDDO !derivIdx
                ENDDO !nodeIdx
              ENDDO facesLoop
            CASE DEFAULT
              CALL FlagError("The dimension is invalid for point Neumann conditions.",err,error,*999)
            END SELECT
          CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The interpolation type of "//TRIM(NumberToVString(interpolationType,"*",err,error))// &
              & " is invalid for component number "//TRIM(NumberToVString(componentNumber,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ENDIF
      ENDDO !neumannDOFIdx

      CALL DistributedMatrix_UpdateStart(neumannConditions%integrationMatrix,err,error,*999)
      CALL DistributedMatrix_UpdateFinish(neumannConditions%integrationMatrix,err,error,*999)
      
      CALL FieldVariable_ParameterSetVectorGet(rhsVariable,FIELD_INTEGRATED_NEUMANN_SET_TYPE,integratedValues,err,error,*999)
      CALL DistributedVector_AllValuesSet(integratedValues,0.0_DP,err,error,*999)
      ! Perform matrix multiplication, f = N q, to calculate force vector from integration matrix and point values
      CALL DistributedMatrix_MatrixByVectorAdd(neumannConditions%integrationMatrix,.FALSE.,neumannConditions%pointValues, &
        & DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,1.0_DP,integratedValues,err,error,*999)


      CALL FieldVariable_ParameterSetUpdateStart(rhsVariable,FIELD_INTEGRATED_NEUMANN_SET_TYPE,err,error,*999)
      
      IF(diagnostics1) THEN
        IF(dependentGeometry) THEN
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Using dependent field geometry",err,error,*999)
        ELSE
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Using undeformed geometry",err,error,*999)
        END IF
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfNeumann,6,6,neumannConditions%setDofs, &
          & '("  setDofs:",6(X,I8))', '(10X,6(X,I8))',err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Neumann point values",err,error,*999)
        CALL DistributedVector_Output(DIAGNOSTIC_OUTPUT_TYPE,neumannConditions%pointValues,err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Neumann integration matrix",err,error,*999)
        CALL DistributedMatrix_Output(DIAGNOSTIC_OUTPUT_TYPE,neumannConditions%integrationMatrix,err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Integrated values",err,error,*999)
        CALL DistributedVector_Output(DIAGNOSTIC_OUTPUT_TYPE,integratedValues,err,error,*999)
      END IF
      
      CALL FieldVariable_ParameterSetUpdateFinish(rhsVariable,FIELD_INTEGRATED_NEUMANN_SET_TYPE,err,error,*999)

    ENDIF !Neumann conditions associated

    EXITS("BoundaryConditionsVariable_NeumannIntegrate")
    RETURN
999 ERRORSEXITS("BoundaryConditionsVariable_NeumannIntegrate",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsVariable_NeumannIntegrate

  !
  !================================================================================================================================
  !

  !Finalise the Neumann condition matrices for Neumann boundary conditions
  SUBROUTINE BoundaryConditionsVariable_NeumannMatricesFinalise(boundaryConditionsNeumann,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsNeumannType), POINTER :: boundaryConditionsNeumann !<A pointer to the Neumann boundary conditions to finalise Neumann condition matrices for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("BoundaryConditionsVariable_NeumannMatricesFinalise",err,error,*999)

    IF(ASSOCIATED(boundaryConditionsNeumann)) THEN
      IF(ASSOCIATED(boundaryConditionsNeumann%integrationMatrix)) &
        & CALL DistributedMatrix_Destroy(boundaryConditionsNeumann%integrationMatrix,err,error,*999)
      IF(ASSOCIATED(boundaryConditionsNeumann%pointValues)) &
        & CALL DistributedVector_Destroy(boundaryConditionsNeumann%pointValues,err,error,*999)
      CALL DomainMapping_Finalise(boundaryConditionsNeumann%pointDofMapping,err,error,*999)
    ENDIF

    EXITS("BoundaryConditionsVariable_NeumannMatricesFinalise")
    RETURN
999 ERRORSEXITS("BoundaryConditionsVariable_NeumannMatricesFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsVariable_NeumannMatricesFinalise

  !
  !================================================================================================================================
  !

  !>Initialise the Neumann boundary conditions matrices and vectors. This must be done after we know which DOFs have Neumann point conditions so that we can work out the matrix sparsity pattern.
  SUBROUTINE BoundaryConditionsVariable_NeumannMatricesInitialise(boundaryConditionsVariable,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to initialise Neumann condition matrices for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnDOF,columnNodeNumber,componentNumber,derivIdx,derivativeNumber,dofType,dofTypeIdx,domainIdx, &
      & domainNumber,dummyErr,faceIdx,faceNumber,globalDOF,interpolationType,lineIdx,lineNumber,localDOF, &
      & localNeumannConditionIdx,localNodeIdx,localType,myGroupNodeNumber,neumannConditionNumber,neumannIdx,nodeIdx,nodeNumber, &
      & numberOfDimensions,numberOfDomains,numberOfGroupNodes,numberOfNodeDerivatives,numberOfNodeFaces,numberOfNodeLines, &
      & numberOfNodes,numberNonZeros,numberOfPointDOFs,numberRowEntries,numberOfVariables,totalNumberOfLocal,versionNumber
    INTEGER(INTG), ALLOCATABLE :: rowIndices(:), columnIndices(:), localDOFNumbers(:)
    REAL(DP) :: pointValue
    LOGICAL :: boundaryFace,boundaryLine,boundaryNode
    TYPE(BasisType), POINTER :: basis
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions
    TYPE(BoundaryConditionsNeumannType), POINTER :: boundaryConditionsNeumann
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainMappingType), POINTER :: rowsMapping, pointDOFMapping
    TYPE(DomainFacesType), POINTER :: domainFaces
    TYPE(DomainLinesType), POINTER :: domainLines
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping
    TYPE(FieldType), POINTER :: field
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(ListType), POINTER :: columnIndicesList, rowColumnIndicesList
    TYPE(VARYING_STRING) :: dummyError,localError
    TYPE(WorkGroupType), POINTER :: workGroup

    ENTERS("BoundaryConditionsVariable_NeumannMatricesInitialise",err,error,*999)

    NULLIFY(fieldVariable)
    CALL BoundaryConditionsVariable_FieldVariableGet(boundaryConditionsVariable,fieldVariable,err,error,*999)
    numberOfPointDofs=boundaryConditionsVariable%dofCounts(BOUNDARY_CONDITION_NEUMANN_POINT) + &
      & boundaryConditionsVariable%dofCounts(BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)
    NULLIFY(boundaryConditionsNeumann)
    CALL BoundaryConditionsVariable_NeumannConditionsGet(boundaryConditionsVariable,boundaryConditionsNeumann,err,error,*999)
!!TODO: We need to access the LHS mapping for the equations set. Just use variable for now. 
    !For rows we can re-use the RHS variable row mapping
    NULLIFY(rowsMapping)
    CALL FieldVariable_DomainMappingGet(fieldVariable,rowsMapping,err,error,*999)
    NULLIFY(workGroup)
    CALL DomainMapping_WorkGroupGet(rowsMapping,workGroup,err,error,*999)
    CALL WorkGroup_NumberOfGroupNodesGet(workGroup,numberOfGroupNodes,err,error,*999)

    !Create a domain mapping for the Neumann point DOFs, required for the distributed matrix columns
    NULLIFY(pointDOFMapping)
    CALL DomainMapping_Initialise(pointDofMapping,err,error,*999)
    CALL DomainMapping_WorkGroupSet(pointDofMapping,workGroup,err,error,*999)
    boundaryConditionsNeumann%pointDofMapping=>pointDofMapping
    !Calculate global to local mapping for Neumann DOFs
    pointDofMapping%numberOfGlobal=numberOfPointDofs
    ALLOCATE(pointDofMapping%globalToLocalMap(numberOfPointDofs),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate Neumann point DOF global to local mapping.",err,error,*999)
    ALLOCATE(localDofNumbers(0:numberOfGroupNodes-1),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate local Neumann DOF numbers.",err,error,*999)
    localDofNumbers=0
    
    IF(diagnostics2) CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Local numbering:",err,error,*999)
    
    DO neumannIdx=1,numberOfPointDofs
      globalDOF=boundaryConditionsNeumann%setDofs(neumannIdx)
      !Get domain information from the RHS variable domain mapping, but set new local numbers.
      CALL DomainMapping_NumberOfDomainsFromGlobalGet(rowsMapping,globalDOF,numberOfDomains,err,error,*999)
      pointDofMapping%globalToLocalMap(neumannIdx)%numberOfDomains=numberOfDomains
      ALLOCATE(pointDofMapping%globalToLocalMap(neumannIdx)%localNumber(numberOfDomains),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate Neumann DOF global to local map local number.",err,error,*999)
      ALLOCATE(pointDofMapping%globalToLocalMap(neumannIdx)%domainNumber(numberOfDomains),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate Neumann DOF global to local map domain number.",err,error,*999)
      ALLOCATE(pointDofMapping%globalToLocalMap(neumannIdx)%localType(numberOfDomains),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate Neumann DOF global to local map local type.",err,error,*999)
      IF(diagnostics2) CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Neumann point DOF index = ",neumannIdx,err,error,*999)
      DO domainIdx=1,numberOfDomains
        CALL DomainMapping_DomainNumberFromGlobalGet(rowsMapping,globalDOF,domainIdx,domainNumber,err,error,*999)
        CALL DomainMapping_LocalTypeFromGlobalGet(rowsMapping,globalDOF,domainIdx,localType,err,error,*999)
        pointDofMapping%globalToLocalMap(neumannIdx)%domainNumber(domainIdx)=domainNumber
        pointDofMapping%globalToLocalMap(neumannIdx)%localType(domainIdx)=localType
        IF(localType==DOMAIN_LOCAL_INTERNAL.OR.localType==DOMAIN_LOCAL_BOUNDARY) THEN
          localDofNumbers(domainNumber)=localDofNumbers(domainNumber)+1
          pointDofMapping%globalToLocalMap(neumannIdx)%localNumber(domainIdx)=localDofNumbers(domainNumber)
          IF(diagnostics2) THEN
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Global RHS variable DoF = ",globalDOF,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Domain number = ",domainNumber,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Local type = ",localType,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Local number = ",localDofNumbers(domainNumber),err,error,*999)
          ENDIF
        ENDIF
      ENDDO !domainIdx
    ENDDO !neumannIdx
    
    !Local DOFs must be numbered before ghost DOFs, so loop though again, this time numbering GHOST DOFs
    IF(diagnostics2) CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Ghost numbering:",err,error,*999)
    
    DO neumannIdx=1,numberOfPointDofs
      globalDOF=boundaryConditionsNeumann%setDofs(neumannIdx)
      CALL DomainMapping_NumberOfDomainsFromGlobalGet(rowsMapping,globalDOF,numberOfDomains,err,error,*999)
      IF(diagnostics2) CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Neumann point DOF index = ",neumannIdx,err,error,*999)
      DO domainIdx=1,numberOfDomains
        CALL DomainMapping_LocalTypeFromGlobalGet(rowsMapping,neumannIdx,domainIdx,localType,err,error,*999)
        IF(localDOF==DOMAIN_LOCAL_GHOST) THEN
          CALL DomainMapping_DomainNumberFromGlobalGet(rowsMapping,globalDOF,domainIdx,domainNumber,err,error,*999)
          localDofNumbers(domainNumber)=localDofNumbers(domainNumber)+1
          pointDofMapping%globalToLocalMap(neumannIdx)%localNumber(domainIdx)=localDofNumbers(domainNumber)
          IF(diagnostics2) THEN
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Global rhs var DOF = ",globalDOF,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Domain number = ",domainNumber,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Local number = ",localDofNumbers(domainNumber),err,error,*999)
          ENDIF
        ENDIF
      ENDDO !domainIdx
    ENDDO !neumannIdx
    
    CALL DomainMapping_LocalFromGlobalCalculate(pointDofMapping,err,error,*999)

    CALL DistributedMatrix_CreateStart(rowsMapping,pointDofMapping,boundaryConditionsNeumann%integrationMatrix,err,error,*999)
    NULLIFY(boundaryConditions)
    CALL BoundaryConditionsVariable_BoundaryConditionsGet(boundaryConditionsVariable,boundaryConditions,err,error,*999)
    SELECT CASE(boundaryConditions%neumannMatrixSparsity)
    CASE(BOUNDARY_CONDITION_SPARSE_MATRICES)
      ! Work out integration matrix sparsity structure
      ! For a single process, compressed column would be more memory efficient, but with
      ! multiple processes the number of Neumann point DOFs could be more than the number
      ! of local row DOFs, and multiplying a compressed row matrix by a vector is faster,
      ! so we will use compressed row storage
      CALL DomainMapping_TotalNumberOfLocalGet(rowsMapping,totalNumberOfLocal,err,error,*999)
      ALLOCATE(rowIndices(totalNumberOfLocal+1),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate Neumann integration matrix column indices.",err,error,*999)
      ! We don't know the number of non zeros before hand, so use a list to keep track of column indices
      NULLIFY(columnIndicesList)
      CALL List_CreateStart(columnIndicesList,err,error,*999)
      CALL List_DataTypeSet(columnIndicesList,LIST_INTG_TYPE,err,error,*999)
      CALL List_CreateFinish(columnIndicesList,err,error,*999)
      ! Stores the column indices for the current row
      NULLIFY(rowColumnIndicesList)
      CALL List_CreateStart(rowColumnIndicesList,err,error,*999)
      CALL List_DataTypeSet(rowColumnIndicesList,LIST_INTG_TYPE,err,error,*999)
      CALL List_MutableSet(rowColumnIndicesList,.TRUE.,err,error,*999)
      CALL List_CreateFinish(rowColumnIndicesList,err,error,*999)
      rowIndices(1)=1

      DO localDof=1,totalNumberOfLocal
        CALL FieldVariable_DOFTypeGet(fieldVariable,localDOF,dofType,dofTypeIdx,err,error,*999)
        SELECT CASE(dofType)
        CASE(FIELD_CONSTANT_INTERPOLATION)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(FIELD_NODE_BASED_INTERPOLATION)
          CALL FieldVariable_DOFParametersGetNode(fieldVariable,dofTypeIdx,versionNumber,derivativeNumber,nodeNumber, &
            & componentNumber,err,error,*999)
        CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The DOF type of "//TRIM(NumberToVString(dofType,"*",err,error))// &
            & " for local DOF "//TRIM(NumberToVString(localDOF,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT

        NULLIFY(domain)
        CALL FieldVariable_ComponentDomainGet(fieldVariable,componentNumber,domain,err,error,*999)
        NULLIFY(domainTopology)
        CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
        CALL FieldVariable_ComponentInterpolationGet(fieldVariable,componentNumber,interpolationType,err,error,*999)
        SELECT CASE(interpolationType)
        CASE(FIELD_CONSTANT_INTERPOLATION)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(FIELD_NODE_BASED_INTERPOLATION)
          NULLIFY(domainNodes)
          CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
          CALL DomainNodes_NodeBoundaryNodeGet(domainNodes,nodeNumber,boundaryNode,err,error,*999)
          IF(boundaryNode) THEN
            CALL Domain_NumberOfDimensionsGet(domain,numberOfDimensions,err,error,*999)
            SELECT CASE(numberOfDimensions)
            CASE(1)
              !Only one column used, as this is the same as setting an integrated value so no other DOFs are affected
              CALL DomainMapping_LocalToGlobalGet(rowsMapping,localDOF,globalDOF,err,error,*999)
              IF(boundaryConditionsVariable%conditionTypes(globalDof)==BOUNDARY_CONDITION_NEUMANN_POINT.OR. &
                & boundaryConditionsVariable%conditionTypes(globalDof)==BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED) THEN
                !Find the Neumann condition number
                neumannConditionNumber=0
                DO neumannIdx=1,numberOfPointDofs
                  IF(boundaryConditionsNeumann%setDofs(neumannIdx)==globalDOF) neumannConditionNumber=neumannIdx
                ENDDO !neumannIdx
                IF(neumannConditionNumber==0) THEN
                  localError="Could not find matching Neuamann condition number for global DOF "// &
                    & TRIM(NumberToVString(globalDof,"*",err,error))//" with a Neumann condition set."
                  CALL FlagError(localError,err,error,*999)
                ENDIF
                CALL List_ItemAdd(rowColumnIndicesList,neumannConditionNumber,err,error,*999)
              ENDIF
            CASE(2)
              !Loop over all lines for this node and find any DOFs that have a Neumann point condition set
              NULLIFY(domainLines)
              CALL DomainTopology_DomainLinesGet(domainTopology,domainLines,err,error,*999)
              CALL DomainNodes_NodeNumberOfLinesGet(domainNodes,nodeNumber,numberOfNodeLines,err,error,*999)
              DO lineIdx=1,numberOfNodeLines
                CALL DomainNodes_NodeLineNumberGet(domainNodes,lineIdx,nodeNumber,lineNumber,err,error,*999)
                CALL DomainLines_LineBoundaryLineGet(domainLines,lineNumber,boundaryLine,err,error,*999)
                IF(.NOT.boundaryLine) CYCLE
                NULLIFY(basis)
                CALL DomainLines_LineBasisGet(domainLines,lineNumber,basis,err,error,*999)
                CALL Basis_NumberOfLocalNodesGet(basis,numberOfNodes,err,error,*999)
                DO localNodeIdx=1,numberOfNodes
                  CALL DomainLines_LineNodeNumberGet(domainLines,nodeIdx,lineNumber,columnNodeNumber,err,error,*999)
                  CALL Basis_NodeNumberOfDerivativesGet(basis,nodeIdx,numberOfNodeDerivatives,err,error,*999)
                  DO derivIdx=1,numberOfNodeDerivatives
                    CALL DomainLines_DerivativeGlobalIndexGet(domainLines,derivIdx,nodeIdx,lineNumber,derivativeNumber, &
                      & err,error,*999)
                    CALL DomainLines_DerivativeVersionNumberGet(domainLines,derivIdx,nodeIdx,lineNumber,versionNumber, &
                      & err,error,*999)
                    CALL FieldVariable_LocalNodeDOFGet(fieldVariable,versionNumber,derivativeNumber,columnNodeNumber, &
                      & componentNumber,columnDOF,err,error,*999)
                    CALL DomainMapping_LocalToGlobalGet(rowsMapping,columnDOF,globalDOF,err,error,*999)
                    IF(boundaryConditionsVariable%conditionTypes(globalDof)==BOUNDARY_CONDITION_NEUMANN_POINT.OR. &
                      & boundaryConditionsVariable%conditionTypes(globalDof)==BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED) THEN
                      neumannConditionNumber=0
                      DO neumannIdx=1,numberOfPointDofs
                        IF(boundaryConditionsNeumann%setDofs(neumannIdx)==globalDOF) neumannConditionNumber=neumannIdx
                      ENDDO !neumannIdx
                      IF(neumannConditionNumber==0) THEN
                        localError="Could not find matching Neuamann condition number for global DOF "// &
                          & TRIM(NumberToVString(globalDOF,"*",err,error))//" with a Neumann condition set."
                        CALL FlagError(localError,err,error,*999)
                      ENDIF
                      CALL List_ItemAdd(rowColumnIndicesList,neumannConditionNumber,err,error,*999)
                    ENDIF
                  ENDDO !derivIdx
                ENDDO !localNodeIdx
              ENDDO !lineIdx
            CASE(3)
              NULLIFY(domainFaces)
              CALL DomainTopology_DomainFacesGet(domainTopology,domainFaces,err,error,*999)
              CALL DomainNodes_NodeNumberOfFacesGet(domainNodes,nodeNumber,numberOfNodeFaces,err,error,*999)
              ! Loop over all faces for this node and find any DOFs that have a Neumann point condition set
              DO faceIdx=1,numberOfNodeFaces
                CALL DomainNodes_NodeFaceNumberGet(domainNodes,faceIdx,nodeNumber,faceNumber,err,error,*999)
                CALL DomainFaces_FaceBoundaryFaceGet(domainFaces,faceNumber,boundaryface,err,error,*999)
                IF(.NOT.boundaryFace) CYCLE
                NULLIFY(basis)
                CALL DomainFaces_FaceBasisGet(domainFaces,faceNumber,basis,err,error,*999)
                CALL Basis_NumberOfLocalNodesGet(basis,numberOfNodes,err,error,*999)
                DO nodeIdx=1,numberOfNodes
                  CALL DomainFaces_FaceNodeNumberGet(domainFaces,nodeIdx,faceNumber,columnNodeNumber,err,error,*999)
                  CALL Basis_NodeNumberOfDerivativesGet(basis,nodeIdx,numberOfNodeDerivatives,err,error,*999)
                  DO derivIdx=1,numberOfNodeDerivatives
                    CALL DomainFaces_DerivativeGlobalIndexGet(domainFaces,derivIdx,nodeIdx,faceNumber,derivativeNumber, &
                      & err,error,*999)
                    CALL DomainFaces_DerivativeVersionNumberGet(domainFaces,derivIdx,nodeIdx,faceNumber,versionNumber, &
                      & err,error,*999)
                    CALL FieldVariable_LocalNodeDOFGet(fieldVariable,versionNumber,derivativeNumber,columnNodeNumber, &
                      & componentNumber,columnDOF,err,error,*999)
                    CALL DomainMapping_LocalToGlobalGet(rowsMapping,columnDOF,globalDOF,err,error,*999)
                    IF(boundaryConditionsVariable%conditionTypes(globalDOF)==BOUNDARY_CONDITION_NEUMANN_POINT.OR. &
                      & boundaryConditionsVariable%conditionTypes(globalDOF)==BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED) THEN
                      neumannConditionNumber=0
                      DO neumannIdx=1,numberOfPointDofs
                        IF(boundaryConditionsNeumann%setDofs(neumannIdx)==globalDof) neumannConditionNumber=neumannIdx
                      ENDDO !neumannIdx
                      IF(neumannConditionNumber==0) THEN
                        localError="Could not find matching Neuamann condition number for global DOF "// &
                          & TRIM(NumberToVString(globalDOF,"*",err,error))//" with a Neumann condition set."
                        CALL FlagError(localError,err,error,*999)
                      ENDIF
                      CALL List_ItemAdd(rowColumnIndicesList,neumannConditionNumber,err,error,*999)
                    ENDIF
                  ENDDO !derivIdx
                ENDDO !nodeIdx
              ENDDO !faceIdx
            CASE DEFAULT
              localError="The domain dimension of "//TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
                & " is invalid for point Neumann conditions."
              CALL FlagError(localError,err,error,*999)
            END SELECT !number of dimensions
          ENDIF !boundaryNode
        CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The interpolation type of "//TRIM(NumberToVString(interpolationType,"*",err,error))// &
            & " is invalid for component number "//TRIM(NumberToVString(componentNumber,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        END SELECT

        !Sort and remove duplicates
        CALL List_RemoveDuplicates(rowColumnIndicesList,err,error,*999)
        !Now add all column DOFs in this row that use Neumann conditions to the overall column indices
        CALL List_AppendList(columnIndicesList,rowColumnIndicesList,err,error,*999)
        CALL List_NumberOfItemsGet(rowColumnIndicesList,numberRowEntries,err,error,*999)
        rowIndices(localDof+1)=rowIndices(localDof)+numberRowEntries
        CALL List_ClearItems(rowColumnIndicesList,err,error,*999)
      ENDDO !local DOFs

      CALL List_Destroy(rowColumnIndicesList,err,error,*999)
      CALL List_DetachAndDestroy(columnIndicesList,numberNonZeros,columnIndices,err,error,*999)
      
      IF(diagnostics1) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Neumann integration matrix sparsity:",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Number non-zeros = ", numberNonZeros,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Number columns = ",numberOfPointDofs,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Number rows = ",totalNumberOfLocal,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfPointDofs+1,6,6,rowIndices, &
          & '("  Row indices: ",6(X,I6))', '(6X,6(X,I6))',err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberNonZeros,6,6,columnIndices, &
          & '("  Column indices: ",6(X,I6))', '(6X,6(X,I6))',err,error,*999)
      ENDIF

      CALL DistributedMatrix_StorageTypeSet(boundaryConditionsNeumann%integrationMatrix, &
        & DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,err,error,*999)
      CALL DistributedMatrix_NumberOfNonZerosSet(boundaryConditionsNeumann%integrationMatrix,numberNonZeros,err,error,*999)
      CALL DistributedMatrix_StorageLocationsSet(boundaryConditionsNeumann%integrationMatrix,rowIndices, &
        & columnIndices(1:numberNonZeros),err,error,*999)
      
      DEALLOCATE(localDOFNumbers)
      DEALLOCATE(rowIndices)
      DEALLOCATE(columnIndices)
    CASE(BOUNDARY_CONDITION_FULL_MATRICES)
      CALL DistributedMatrix_StorageTypeSet(boundaryConditionsNeumann%integrationMatrix,DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE, &
        & err,error,*999)
    CASE DEFAULT
      localError="The Neumann matrix sparsity type of "// &
        & TRIM(NumberToVString(boundaryConditions%neumannMatrixSparsity,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    CALL DistributedMatrix_TransposeTypeSet(boundaryConditionsNeumann%integrationMatrix, &
      & DISTRIBUTED_MATRIX_NO_TRANSPOSE_REQUIRED,err,error,*999)
    CALL DistributedMatrix_CreateFinish(boundaryConditionsNeumann%integrationMatrix,err,error,*999)
    
    !Set up vector of Neumann point values
    CALL DistributedVector_CreateStart(pointDofMapping,boundaryConditionsNeumann%pointValues,err,error,*999)
    CALL DistributedVector_CreateFinish(boundaryConditionsNeumann%pointValues,err,error,*999)
    NULLIFY(field)
    CALL FieldVariable_FieldGet(fieldVariable,field,err,error,*999)
    NULLIFY(decomposition)
    CALL Field_DecompositionGet(field,decomposition,err,error,*999)
    NULLIFY(workGroup)
    CALL Decomposition_WorkGroupGet(decomposition,workGroup,err,error,*999)
    CALL WorkGroup_GroupNodeNumberGet(workGroup,myGroupNodeNumber,err,error,*999)
    !Set point values vector from boundary conditions field parameter set
    DO neumannIdx=1,numberOfPointDofs
      globalDOF=boundaryConditionsNeumann%setDofs(neumannIdx)
      CALL DomainMapping_DomainNumberFromGlobalGet(rowsMapping,globalDOF,1,domainNumber,err,error,*999)
      IF(domainNumber==myGroupNodeNumber) THEN
        CALL DomainMapping_LocalNumberFromGlobalGet(rowsMapping,globalDOF,1,localDOF,err,error,*999)
        ! Set point DOF vector value
        CALL DomainMapping_LocalNumberFromGlobalGet(pointDOFMapping,neumannIdx,1,localNeumannConditionIdx,err,error,*999)
        CALL FieldVariable_ParameterSetGetLocalDOF(fieldVariable,FIELD_BOUNDARY_CONDITIONS_SET_TYPE,localDOF,pointValue, &
          & err,error,*999)
        CALL DistributedVector_ValuesSet(boundaryConditionsNeumann%pointValues,localNeumannConditionIdx,pointValue,err,error,*999)
      ENDIF
    ENDDO !neumannIdx
    CALL DistributedVector_UpdateStart(boundaryConditionsNeumann%pointValues,err,error,*999)
    CALL DistributedVector_UpdateFinish(boundaryConditionsNeumann%pointValues,err,error,*999)

    EXITS("BoundaryConditionsVariable_NeumannMatricesInitialise")
    RETURN
999 IF(ALLOCATED(rowIndices)) DEALLOCATE(rowIndices)
    IF(ALLOCATED(columnIndices)) DEALLOCATE(columnIndices)
    IF(ALLOCATED(localDofNumbers)) DEALLOCATE(localDofNumbers)
    CALL BoundaryConditionsVariable_NeumannMatricesFinalise(boundaryConditionsVariable%neumannBoundaryConditions, &
      & dummyErr,dummyError,*998)
998 ERRORS("BoundaryConditionsVariable_NeumannMatricesInitialise",err,error)
    EXITS("BoundaryConditionsVariable_NeumannMatricesInitialise")
    RETURN 1

  END SUBROUTINE BoundaryConditionsVariable_NeumannMatricesInitialise

  !
  !================================================================================================================================
  !

  !>Finalises the pressure incremented boundary condition and deallocates all memory
  SUBROUTINE BoundaryConditionsVariable_PressureIncrementedFinalise(boundaryConditionsPressureIncremented,err,error,*)
    !Argument variables
    TYPE(BoundaryConditionsPressureIncrementedType), POINTER :: boundaryConditionsPressureIncremented !<A pointer to the pressure incremented boundary conditions to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("BoundaryConditionsVariable_PressureIncrementedFinalise",err,error,*999)

    IF(ASSOCIATED(boundaryConditionsPressureIncremented)) THEN
      IF(ALLOCATED(boundaryConditionsPressureIncremented%pressureIncrementedDOFIndices)) &
        & DEALLOCATE(boundaryConditionsPressureIncremented%pressureIncrementedDOFIndices)
      DEALLOCATE(boundaryConditionsPressureIncremented)
    ENDIF    

    EXITS("BoundaryConditionsVariable_PressureIncrementedFinalise")
    RETURN
999 ERRORS("BoundaryConditionsVariable_PressureIncrementedFinalise",err,error)
    EXITS("BoundaryConditionsVariable_PressureIncrementedFinalise")
    RETURN 1

  END SUBROUTINE BoundaryConditionsVariable_PressureIncrementedFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the pressure incremented boundary condition.
  SUBROUTINE BoundaryConditionsVariable_PressureIncrementedInitialise(boundaryConditionsVariable,err,error,*)
    !Argument variables
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to initialise a boundary conditions dirichlet type for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,numberOfPressureIncrementedConditions
    TYPE(BoundaryConditionsPressureIncrementedType), POINTER :: boundaryConditionsPressureIncremented
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("BoundaryConditionsVariable_PressureIncrementedInitialise",err,error,*999)

    IF(.NOT.ASSOCIATED(boundaryConditionsVariable)) &
      & CALL FlagError("Boundary conditions variable is not associated.",err,error,*999)
    IF(ASSOCIATED(boundaryConditionsVariable%pressureIncrementedBoundaryConditions)) &
      & CALL FlagError("Pressure incremented boundary conditions are already associated for this boundary conditions variable." &
      & ,err,error,*999)
      
    ALLOCATE(boundaryConditionsVariable%pressureIncrementedBoundaryConditions,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate pressure incremented Boundary Conditions",err,error,*999)
    boundaryConditionsPressureIncremented=>boundaryConditionsVariable%pressureIncrementedBoundaryConditions
    numberOfPressureIncrementedConditions=boundaryConditionsVariable%dofCounts(BOUNDARY_CONDITION_PRESSURE_INCREMENTED)
    ALLOCATE(boundaryConditionsPressureIncremented%pressureIncrementedDOFIndices(numberOfPressureIncrementedConditions),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate Pressure incremented DOF indices array",err,error,*999)

    EXITS("BoundaryConditionsVariable_PressureIncrementedInitialise")
    RETURN
999 CALL BoundaryConditionsVariable_PressureIncrementedFinalise(boundaryConditionsVariable%pressureIncrementedBoundaryConditions, &
      & dummyErr,dummyError,*998)
998 ERRORS("BoundaryConditionsVariable_PressureIncrementedInitialise",err,error)
    EXITS("BoundaryConditionsVariable_PressureIncrementedInitialise")
    RETURN 1

  END SUBROUTINE BoundaryConditionsVariable_PressureIncrementedInitialise

  !
  !================================================================================================================================
  !

END MODULE BoundaryConditionsRoutines
