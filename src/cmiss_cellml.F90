!> \file
!> \author David Nickerson <nickerso@users.sourceforge.net>
!> \brief This module is a OpenCMISS(cm) buffer module to OpenCMISS(cellml).
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
!> Contributor(s): David Nickerson, Chris Bradley
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

!> This module is a OpenCMISS(cmfe) buffer module to OpenCMISS(cellml).
MODULE CmissCellML

  !Module imports
  USE ISO_C_BINDING

  USE BaseRoutines
  USE CellMLAccessRoutines
  
#ifdef WITH_CELLML
  USE CELLML_MODEL_DEFINITION
#endif
  
  ! Moved this usage declaration outside the preprocessor definition,
  ! as this file is included only if CELLML integration is selected.
  ! This fixes problems with the CMAKE FORTRAN parser. Its not detecting
  ! the file (=module) dependency correctly and hence breaks the build
  ! on some platforms.
  USE CMISSFortranC
  USE CmissMPI
  USE ComputationRoutines
  USE ComputationAccessRoutines
  USE Constants
  USE DecompositionAccessRoutines
  USE FieldRoutines
  USE FieldAccessRoutines
  USE ISO_VARYING_STRING
  USE InputOutput
  USE Kinds
  USE MeshAccessRoutines
#ifndef NOMPIMOD
  USE MPI
#endif
  USE RegionAccessRoutines
  USE Strings
  USE Types

#include "macros.h"

  IMPLICIT NONE

  PRIVATE

#ifdef NOMPIMOD
#include "mpif.h"
#endif

  !Module parameters

  !Module types

  !Module variables
 
  !Interfaces

  INTERFACE CellML_ModelImport
    MODULE PROCEDURE CellML_ModelImportC
    MODULE PROCEDURE CellML_ModelImportVS
  END INTERFACE CellML_ModelImport
  
  INTERFACE CellML_VariableSetAsKnown
    MODULE PROCEDURE CellML_VariableSetAsKnownC
    MODULE PROCEDURE CellML_VariableSetAsKnownVS
  END INTERFACE CellML_VariableSetAsKnown

  INTERFACE CellML_VariableSetAsWanted
    MODULE PROCEDURE CellML_VariableSetAsWantedC
    MODULE PROCEDURE CellML_VariableSetAsWantedVS
  END INTERFACE CellML_VariableSetAsWanted

  INTERFACE CellML_CreateCellMLToFieldMap
    MODULE PROCEDURE CellML_CreateCellMLToFieldMapC
    MODULE PROCEDURE CellML_CreateCellMLToFieldMapVS
  END INTERFACE CellML_CreateCellMLToFieldMap

  INTERFACE CellML_CreateFieldToCellMLMap
    MODULE PROCEDURE CellML_CreateFieldToCellMLMapC
    MODULE PROCEDURE CellML_CreateFieldToCellMLMapVS
  END INTERFACE CellML_CreateFieldToCellMLMap

  INTERFACE CellML_FieldComponentGet
    MODULE PROCEDURE CellML_FieldComponentGetC
    MODULE PROCEDURE CellML_FieldComponentGetVS
  END INTERFACE CellML_FieldComponentGet

  PUBLIC CellML_CellMLToFieldUpdate
  
  PUBLIC CellML_CreateStart,CellML_CreateFinish

  PUBLIC CellML_Destroy

  PUBLIC CellML_FieldToCellMLUpdate
  
  PUBLIC CellML_ModelImport

  PUBLIC CellML_VariableSetAsKnown,CellML_VariableSetAsWanted

  PUBLIC CellML_FieldMapsCreateStart,CellML_FieldMapsCreateFinish

  PUBLIC CellML_FieldModelDofSet
  
  PUBLIC CellML_CreateCellMLToFieldMap,CellML_CreateFieldToCellMLMap

  PUBLIC CellML_ModelsFieldCreateStart,CellML_ModelsFieldCreateFinish

  PUBLIC CellML_StateFieldCreateStart,CellML_StateFieldCreateFinish

  PUBLIC CellML_FieldComponentGet

  PUBLIC CellML_IntermediateFieldCreateFinish,CellML_IntermediateFieldCreateStart

  PUBLIC CellML_ParametersFieldCreateStart,CellML_ParametersFieldCreateFinish

  PUBLIC CellML_Generate

  PUBLIC CellMLEnvironments_Finalise,CellMLEnvironments_Initialise
  
CONTAINS

  !
  !=================================================================================================================================
  !


  !>Updates any mapped fields from the cellml fields
  SUBROUTINE CellML_CellMLToFieldUpdate(cellML,err,error,*)
    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML environment to update
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !< The error string
    !Local variables
    INTEGER(INTG) :: derivativeNumber,dofIdx,dofType,dof2ParamIdx,elementNumber,gaussNumber,gridNumber,mapIdx,modelIdx, &
      & nodeNumber,versionNumber
    INTEGER(INTG), POINTER :: modelsData(:)
    REAL(DP) :: dofValue
    TYPE(CellMLFieldMapsType), POINTER :: cellMLFieldMaps
    TYPE(CellMLIntermediateFieldType), POINTER :: cellMLIntermediateField
    TYPE(CellMLModelMapType), POINTER :: cellMLModelMap
    TYPE(CellMLModelMapsType), POINTER :: cellMLModelMaps
    TYPE(CellMLModelsFieldType), POINTER :: cellMLModelsField
    TYPE(CellMLParametersFieldType), POINTER :: cellMLParametersField
    TYPE(CellMLStateFieldType), POINTER :: cellMLStateField
    TYPE(FieldType), POINTER :: modelsField
    TYPE(FieldVariableType), POINTER :: modelsVariable
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("CellML_CellMLToFieldUpdate",err,error,*999)

    NULLIFY(cellMLModelsField)
    CALL CellML_CellMLModelsFieldGet(cellML,cellMLModelsField,err,error,*999)
    NULLIFY(cellMLFieldMaps)
    CALL CellML_CellMLFieldMapsGet(cellML,cellMLFieldMaps,err,error,*999)
    
    IF(cellMLModelsField%onlyOneModelIndex/=0) THEN
      !We have CellML models on this rank
      IF(cellMLModelsField%onlyOneModelIndex/=CELLML_MODELS_FIELD_NOT_CONSTANT) THEN
        !The CellML environement only uses one model and so we can optimise for this.
        NULLIFY(cellMLModelMaps)
        CALL CellMLFieldMaps_CellMLModelMapsGet(cellMLFieldMaps,cellMLModelsField%onlyOneModelIndex,cellMLModelMaps,err,error,*999)
        IF(diagnostics1) THEN
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"CellML to field update:",err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  CellML user number = ",cellML%userNumber,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  One model index = ",cellMLModelsField%onlyOneModelIndex,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of model maps = ",cellMLModelMaps%numberOfFieldsMappedTo, &
            & err,error,*999)
        ENDIF
        !Loop over the number of CellML to field maps
        DO mapIdx=1,cellMLModelMaps%numberOfFieldsMappedFrom
          NULLIFY(cellMLModelMap)
          CALL CellMLModelMaps_CellMLModelFromMapGet(cellMLModelMaps,mapIdx,cellMLModelMap,err,error,*999)
          SELECT CASE(cellMLModelMap%cellMLFieldType)
          CASE(CELLML_MODELS_FIELD)
            CALL FlagError("Cannot map models field.",err,error,*999)
          CASE(CELLML_STATE_FIELD)
            NULLIFY(cellMLStateField)
            CALL CellML_CellMLStateFieldGet(cellML,cellMLStateField,err,error,*999)
            CALL Field_ParametersToFieldParametersCopy(cellMLStateField%stateField, &
              & FIELD_U_VARIABLE_TYPE,cellMLModelMap%cellMLParameterSet,cellMLModelMap%cellMLVariableNumber, &
              & cellMLModelMap%field,cellMLModelMap%variableType,cellMLModelMap%fieldParameterSet, &
              & cellMLModelMap%componentNumber,err,error,*999)
          CASE(CELLML_INTERMEDIATE_FIELD)
            NULLIFY(cellMLIntermediateField)
            CALL CellML_CellMLIntermediateFieldGet(cellML,cellMLIntermediateField,err,error,*999)
            CALL Field_ParametersToFieldParametersCopy(cellMLIntermediateField%intermediateField, &
              & FIELD_U_VARIABLE_TYPE,cellMLModelMap%cellMLParameterSet,cellMLModelMap%cellMLVariableNumber, &
              & cellMLModelMap%field,cellMLModelMap%variableType,cellMLModelMap%fieldParameterSet, &
              & cellMLModelMap%componentNumber,err,error,*999)
          CASE(CELLML_PARAMETERS_FIELD)
            NULLIFY(cellMLParametersField)
            CALL CellML_CellMLParametersFieldGet(cellML,cellMLParametersField,err,error,*999)
            CALL Field_ParametersToFieldParametersCopy(cellMLParametersField%parametersField, &
              & FIELD_U_VARIABLE_TYPE,cellMLModelMap%cellMLParameterSet,cellMLModelMap%cellMLVariableNumber, &
              & cellMLModelMap%field,cellMLModelMap%variableType,cellMLModelMap%fieldParameterSet, &
              & cellMLModelMap%componentNumber,err,error,*999)
          CASE DEFAULT
            localError="The CellML to field model map CellML field type of "// &
              & TRIM(NumberToVString(cellMLModelMap%cellMLFieldType,"*",err,error))// &
              & " is invalid for map index "//TRIM(NumberToVString(mapIdx,"*",err,error))//" of model index "// &
              & TRIM(NumberToVString(cellMLModelsField%onlyOneModelIndex,"*",err,error))// &
              & " of CellML user number "//TRIM(NumberToVString(cellML%userNumber,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          IF(diagnostics1) THEN
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Map index : ",mapIdx,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    CellML field type      = ",cellMLModelMap%cellMLFieldType, &
              & err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    CellML parameter set   = ",cellMLModelMap%cellMLParameterSet, &
              & err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    CellML variable number = ",cellMLModelMap%cellMLVariableNumber, &
              & err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Field user number      = ",cellMLModelMap%field%userNumber, &
              & err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Field variable type    = ",cellMLModelMap%variableType, &
              & err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Field parameter set    = ",cellMLModelMap%fieldParameterSet, &
              & err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Field component number = ",cellMLModelMap%componentNumber, &
              & err,error,*999)                   
          ENDIF
        ENDDO !mapIdx
      ELSE
        !More than one model used in the models field.
        NULLIFY(modelsField)
        CALL CellMLModelsField_ModelsFieldGet(cellMLModelsField,modelsField,err,error,*999)
        NULLIFY(modelsVariable)
        CALL Field_VariableGet(modelsField,FIELD_U_VARIABLE_TYPE,modelsVariable,err,error,*999)
        NULLIFY(modelsData)
        CALL FieldVariable_ParameterSetDataGet(modelsVariable,FIELD_VALUES_SET_TYPE,modelsData,err,error,*999)
        IF(diagnostics1) THEN
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"CellML to field update:",err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  CellML user number = ",cellML%userNumber,err,error,*999)
        ENDIF
        !Loop over the dofs in the models field
        DO dofIdx=1,modelsVariable%totalNumberOfDofs
          modelIdx=modelsData(dofIdx)
          IF(modelIdx>0) THEN
            NULLIFY(cellMLModelMaps)
            CALL CellMLFieldMaps_CellMLModelMapsGet(cellMLFieldMaps,modelIdx,cellMLModelMaps,err,error,*999)
            dofType=modelsVariable%dofToParamMap%DOFType(1,dofIdx)
            dof2ParamIdx=modelsVariable%dofToParamMap%DOFType(2,dofIdx)
            SELECT CASE(dofType)
            CASE(FIELD_CONSTANT_INTERPOLATION)
              !Do nothing
            CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
              elementNumber=modelsVariable%dofToParamMap%elementDOF2ParamMap(1,dof2ParamIdx)
            CASE(FIELD_NODE_BASED_INTERPOLATION)
              versionNumber=modelsVariable%dofToParamMap%nodeDOF2ParamMap(1,dof2ParamIdx)
              derivativeNumber=modelsVariable%dofToParamMap%nodeDOF2ParamMap(2,dof2ParamIdx)
              nodeNumber=modelsVariable%dofToParamMap%nodeDOF2ParamMap(3,dof2ParamIdx)
            CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
              gridNumber=modelsVariable%dofToParamMap%gridPointDOF2ParamMap(1,dof2ParamIdx)
            CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
              gaussNumber=modelsVariable%dofToParamMap%gaussPointDOF2ParamMap(1,dof2ParamIdx)
              elementNumber=modelsVariable%dofToParamMap%gaussPointDOF2ParamMap(2,dof2ParamIdx)
            CASE DEFAULT
              localError="The DOF type of "//TRIM(NumberToVString(dofType,"*",err,error))// &
                & " for local DOF number "//TRIM(NumberToVString(dofIdx,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
            IF(diagnostics1) THEN
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  DOF index : ",dofIdx,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Model index : ",modelIdx,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      DOF type = ",dofType,err,error,*999)
              SELECT CASE(dofType)
              CASE(FIELD_CONSTANT_INTERPOLATION)
                !Do nothing
              CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Element number = ",elementNumber,err,error,*999)
              CASE(FIELD_NODE_BASED_INTERPOLATION)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Version number = ",versionNumber,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Derivative number = ",derivativeNumber, &
                  & err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Node number = ",nodeNumber,err,error,*999)
              CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Grid number = ",gridNumber,err,error,*999)
              CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Gauss number = ",gaussNumber,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Element number = ",elementNumber,err,error,*999)
              CASE DEFAULT
                localError="The DOF type of "//TRIM(NumberToVString(dofType,"*",err,error))// &
                  & " for local DOF number "//TRIM(NumberToVString(dofIdx,"*",err,error))//" is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of model maps = ",cellMLModelMaps%numberOfFieldsMappedFrom, &
                & err,error,*999)
            ENDIF
            !Loop over the number of CellML to field maps
            DO mapIdx=1,cellMLModelMaps%numberOfFieldsMappedFrom
              NULLIFY(cellMLModelMap)
              CALL CellMLModelMaps_CellMLModelFromMapGet(cellMLModelMaps,mapIdx,cellMLModelMap,err,error,*999)
              !Get the CellML field DOF value
              SELECT CASE(cellMLModelMap%cellMLFieldType)
              CASE(CELLML_MODELS_FIELD)
                CALL FlagError("Cannot map models field.",err,error,*999)
              CASE(CELLML_STATE_FIELD)
                NULLIFY(cellMLStateField)
                CALL CellML_CellMLStateFieldGet(cellML,cellMLStateField,err,error,*999)
                SELECT CASE(dofType)
                CASE(FIELD_CONSTANT_INTERPOLATION)
                  CALL Field_ParameterSetGetConstant(cellMLStateField%stateField,FIELD_U_VARIABLE_TYPE, &
                    & cellMLModelMap%cellMLParameterSet,cellMLModelMap%cellMLVariableNumber,dofValue,err,error,*999)
                CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                  CALL Field_ParameterSetGetLocalElement(cellMLStateField%stateField,FIELD_U_VARIABLE_TYPE, &
                    & cellMLModelMap%cellMLParameterSet,elementNumber,cellMLModelMap%cellMLVariableNumber,dofValue, &
                    & err,error,*999)
                CASE(FIELD_NODE_BASED_INTERPOLATION)
                  CALL Field_ParameterSetGetLocalNode(cellMLStateField%stateField,FIELD_U_VARIABLE_TYPE, &
                    & cellMLModelMap%cellMLParameterSet,versionNumber,derivativeNumber,nodeNumber, &
                    & cellMLModelMap%cellMLVariableNumber,dofValue,err,error,*999)
                CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                  CALL Field_ParameterSetGetLocalGaussPoint(cellMLStateField%stateField,FIELD_U_VARIABLE_TYPE, &
                    & cellMLModelMap%cellMLParameterSet,gaussNumber,elementNumber,cellMLModelMap%cellMLVariableNumber, &
                    & dofValue,err,error,*999)
                CASE DEFAULT
                  localError="The DOF type of "//TRIM(NumberToVString(dofType,"*",err,error))// &
                    & " for local DOF number "//TRIM(NumberToVString(dofIdx,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE(CELLML_INTERMEDIATE_FIELD)
                NULLIFY(cellMLIntermediateField)
                CALL CellML_CellMLIntermediateFieldGet(cellML,cellMLIntermediateField,err,error,*999)
                SELECT CASE(dofType)
                CASE(FIELD_CONSTANT_INTERPOLATION)
                  CALL Field_ParameterSetGetConstant(cellMLIntermediateField%intermediateField, &
                    & FIELD_U_VARIABLE_TYPE,cellMLModelMap%cellMLParameterSet,cellMLModelMap%cellMLVariableNumber, &
                    & dofValue,err,error,*999)
                CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                  CALL Field_ParameterSetGetLocalElement(cellMLIntermediateField%intermediateField, &
                    & FIELD_U_VARIABLE_TYPE,cellMLModelMap%cellMLParameterSet,elementNumber, &
                    & cellMLModelMap%cellMLVariableNumber,dofValue,err,error,*999)
                CASE(FIELD_NODE_BASED_INTERPOLATION)
                  CALL Field_ParameterSetGetLocalNode(cellMLIntermediateField%intermediateField, &
                    & FIELD_U_VARIABLE_TYPE,cellMLModelMap%cellMLParameterSet,versionNumber,derivativeNumber, &
                    & nodeNumber,cellMLModelMap%cellMLVariableNumber,dofValue,err,error,*999)
                CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                  CALL Field_ParameterSetGetLocalGaussPoint(cellmlIntermediateField%intermediateField, &
                    & FIELD_U_VARIABLE_TYPE,cellMLModelMap%cellMLParameterSet,gaussNumber,elementNumber, &
                    & cellMLModelMap%cellMLVariableNumber,dofValue,err,error,*999)
                CASE DEFAULT
                  localError="The DOF type of "//TRIM(NumberToVString(dofType,"*",err,error))// &
                    & " for local DOF number "//TRIM(NumberToVString(dofIdx,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE(CELLML_PARAMETERS_FIELD)
                NULLIFY(cellMLParametersField)
                CALL CellML_CellMLParametersFieldGet(cellML,cellMLParametersField,err,error,*999)
                SELECT CASE(dofType)
                CASE(FIELD_CONSTANT_INTERPOLATION)
                  CALL Field_ParameterSetGetConstant(cellMLParametersField%parametersField, &
                  & FIELD_U_VARIABLE_TYPE,cellMLModelMap%cellMLParameterSet,cellMLModelMap%cellMLVariableNumber, &
                  & dofValue,err,error,*999)
                CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                  CALL Field_ParameterSetGetLocalElement(cellMLParametersField%parametersField, &
                    & FIELD_U_VARIABLE_TYPE,cellMLModelMap%cellMLParameterSet,elementNumber, &
                    & cellMLModelMap%cellMLVariableNumber,dofValue,err,error,*999)
                CASE(FIELD_NODE_BASED_INTERPOLATION)
                  CALL Field_ParameterSetGetLocalNode(cellMLParametersField%parametersField, &
                    & FIELD_U_VARIABLE_TYPE,cellMLModelMap%cellMLParameterSet,versionNumber,derivativeNumber, &
                    & nodeNumber,cellMLModelMap%cellMLVariableNumber,dofValue,err,error,*999)
                CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                  CALL Field_ParameterSetGetLocalGaussPoint(cellMLParametersField%parametersField, &
                  & FIELD_U_VARIABLE_TYPE,cellMLModelMap%cellMLParameterSet,gaussNumber,elementNumber, &
                  & cellMLModelMap%cellMLVariableNumber,dofValue,err,error,*999)
                CASE DEFAULT
                  localError="The DOF type of "//TRIM(NumberToVString(dofType,"*",err,error))// &
                    & " for local DOF number "//TRIM(NumberToVString(dofIdx,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE DEFAULT
                localError="The CellML to field model map CellML field type of "// &
                  & TRIM(NumberToVString(cellMLModelMap%cellMLFieldType,"*",err,error))// &
                  & " is invalid for map index "//TRIM(NumberToVString(mapIdx,"*",err,error))// &
                  & " of model index "//TRIM(NumberToVString(cellMLModelsField%onlyOneModelIndex,"*",err,error))// &
                  & " of CellML user number "//TRIM(NumberToVString(cellML%userNumber,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              END SELECT
              !Update the OpenCMISS mapped field DOF value
              SELECT CASE(dofType)
              CASE(FIELD_CONSTANT_INTERPOLATION)
                CALL Field_ParameterSetUpdateConstant(cellMLModelMap%field,cellMLModelMap%variableType, &
                  & cellMLModelMap%fieldParameterSet,cellMLModelMap%componentNumber,dofValue,err,error,*999)
              CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                CALL Field_ParameterSetUpdateLocalElement(cellMLModelMap%field,cellMLModelMap%variableType, &
                  & cellMLModelMap%fieldParameterSet,elementNumber,cellMLModelMap%componentNumber,dofValue, &
                  & err,error,*999)
              CASE(FIELD_NODE_BASED_INTERPOLATION)
                CALL Field_ParameterSetUpdateLocalNode(cellMLModelMap%field,cellMLModelMap%variableType, &
                  & cellMLModelMap%fieldParameterSet,versionNumber,derivativeNumber,nodeNumber, &
                  & cellMLModelMap%componentNumber,dofValue,err,error,*999)
              CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                CALL FlagError("Not implemented.",err,error,*999)
              CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                CALL Field_ParameterSetUpdateLocalGaussPoint(cellMLModelMap%field,cellMLModelMap%variableType, &
                  & cellMLModelMap%fieldParameterSet,gaussNumber,elementNumber,cellMLModelMap%componentNumber, &
                  & dofValue,err,error,*999)
              CASE DEFAULT
                localError="The DOF type of "//TRIM(NumberToVString(dofType,"*",err,error))// &
                  & " for local DOF number "//TRIM(NumberToVString(dofIdx,"*",err,error))//" is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
              IF(diagnostics1) THEN
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Map index : ",mapIdx,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        CellML field type      = ", &
                  & cellMLModelMap%cellMLFieldType,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        CellML parameter set   = ", &
                  & cellMLModelMap%cellMLParameterSet,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        CellML variable number = ", &
                  & cellMLModelMap%cellMLVariableNumber,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Field user number      = ", &
                  & cellMLModelMap%field%userNumber,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Field variable type    = ", &
                  & cellMLModelMap%variableType,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Field parameter set    = ", &
                  & cellMLModelMap%fieldParameterSet,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Field component number = ", &
                  & cellMLModelMap%componentNumber,err,error,*999)                   
              ENDIF
            ENDDO !mapIdx
          ENDIF !modelIdx>0
        ENDDO !dofIdx              
        CALL FieldVariable_ParameterSetDataRestore(modelsVariable,FIELD_VALUES_SET_TYPE,modelsData,err,error,*999)
      ENDIF
    ENDIF

    EXITS("CellML_CellMLToFieldUpdate")
    RETURN
999 ERRORSEXITS("CellML_CellMLToFieldUpdate",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_CellMLToFieldUpdate

  !
  !=================================================================================================================================
  !

  !>Set up the CellML environment in the given region.
  !! For a given region, create a CellML environment that will be used to define CellML models in OpenCMISS. This will simply create and initialise an empty CellML environment object in the specified region with the specified unique identifier number. Also set some flag to indicate the CellML environment object is in the process of being created and should not yet be used.
  !! \todo Should be passing in a region? or leave that for the field mapping?
  SUBROUTINE CellML_CreateStart(cellMLUserNumber,region,cellML,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: cellMLUserNumber !<The user specified ID for this CellML environment.
    TYPE(RegionType), POINTER, INTENT(IN) :: region !<The region this CellML environment belongs to.
    TYPE(CellMLType), POINTER :: cellML !<The newly created CellML environment object.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: cellMLIdx
    TYPE(CellMLType), POINTER :: newCellML
    TYPE(CellMLEnvironmentsType), POINTER :: cellMLEnvironments
    TYPE(CellMLPtrType), ALLOCATABLE :: newEnvironments(:)
    TYPE(VARYING_STRING) :: localError

    NULLIFY(newCellML)

    ENTERS("CellML_CreateStart",err,error,*998)

    IF(ASSOCIATED(cellML)) CALL FlagError("CellML is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated.",err,error,*999)
    
    NULLIFY(cellML)
    CALL CellML_UserNumberFind(cellMLUserNumber,region,cellML,err,error,*999)
    IF(ASSOCIATED(cellML)) THEN
      localError="CellML environment number "//TRIM(NumberToVString(cellMLUserNumber,"*",err,error))// &
        & " has already been created on region number "//TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    NULLIFY(cellMLEnvironments)
    CALL Region_CellMLEnvironmentsGet(region,cellMLEnvironments,err,error,*999)
    !Allocate new cellml
    NULLIFY(newCellML)
    CALL CellML_Initialise(newCellML,err,error,*999)
    !Set cellml defaults
    newCellML%userNumber=cellMLUserNumber
    newCellML%globalNumber=cellMLEnvironments%numberOfEnvironments+1
    newCellML%environments=>cellMLEnvironments
    newCellML%region=>region
    !Add cellml to the list of cellml environments
    ALLOCATE(newEnvironments(cellMLEnvironments%numberOfEnvironments+1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new CellML environments.",err,error,*999)
    DO cellMLIdx=1,cellMLEnvironments%numberOfEnvironments
      newEnvironments(cellMLIdx)%ptr=>cellMLEnvironments%environments(cellMLIdx)%ptr
    ENDDO !cellMLIdx
    newEnvironments(cellMLEnvironments%numberOfEnvironments+1)%ptr=>newCellML
    CALL MOVE_ALLOC(newEnvironments,cellMLEnvironments%environments)
    cellMLEnvironments%numberOfEnvironments=cellMLEnvironments%numberOfEnvironments+1
    cellML=>newCellML

    EXITS("CellML_CreateStart")
    RETURN
999 NULLIFY(cellML)
    IF(ALLOCATED(newEnvironments)) DEALLOCATE(newEnvironments)
998 ERRORSEXITS("CellML_CreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_CreateStart

  !
  !=================================================================================================================================
  !

  !>Finish creating the CellML environment.
  !>At this point we know all the variables that are known and wanted so can generate the final code.
  !! Check the provided CellML environment object and if it all looks good clear the "in progress" flag to indicate the object is now ready for use.
  SUBROUTINE CellML_CreateFinish(cellML,err,error,*)
    
    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<The CellML environment object to check and finialise creation of.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(CellMLModelType), POINTER :: cellMLModel
    INTEGER(C_INT) :: errorCode
    INTEGER(INTG) :: modelIdx
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("CellML_CreateFinish",err,error,*999)

    CALL CellML_AssertNotFinished(cellML,err,error,*999)
    !Check that we have set up the models
    IF(cellML%numberOfModels<=0) &
      & CALL FlagError("Invalid setup. No models have been imported into the CellML environment.",err,error,*999)
    
    DO modelIdx=1,cellML%numberOfModels
      NULLIFY(cellMLModel)
      CALL CellML_CellMLModelGet(cellML,modelIdx,cellMLModel,err,error,*999)
      
#ifdef WITH_CELLML
        
      !CALL CELLML_MODEL_DEFINITION_SET_SAVE_TEMP_FILES(cellMLModel%ptr,1)
      errorCode=CELLML_MODEL_DEFINITION_INSTANTIATE(cellMLModel%ptr)
      IF(errorCode/=0) THEN
        localError="Error "//TRIM(NumberToVString(errorCode,"*",err,error))// &
          & " while instantiating CellML model index "//TRIM(NumberToVString(modelIdx,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      cellMLModel%numberOfState=CELLML_MODEL_DEFINITION_GET_N_RATES(cellMLModel%ptr)
      IF(cellMLModel%numberOfState>cellML%maximumNumberOfState) &
        & cellML%maximumNumberOfState=cellMLModel%numberOfState
      cellMLModel%numberOfParameters=CELLML_MODEL_DEFINITION_GET_N_CONSTANTS(cellMLModel%ptr)
      IF(cellMLModel%numberOfParameters>cellML%maximumNumberOfParameters) &
        & cellML%maximumNumberOfParameters=cellMLModel%numberOfParameters
      cellMLModel%numberOfIntermediate=CELLML_MODEL_DEFINITION_GET_N_ALGEBRAIC(cellMLModel%ptr)
      IF(cellMLModel%numberOfIntermediate>cellML%maximumNumberOfIntermediate)  &
        & cellML%maximumNumberOfIntermediate=cellMLModel%numberOfIntermediate
      
#else
      
      CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*999)
      
#endif
        
    ENDDO !modelIdx
    cellML%cellMLFinished=.TRUE.
    
    EXITS("CellML_CreateFinish")
    RETURN
999 ERRORSEXITS("CellML_CreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_CreateFinish

  !
  !=================================================================================================================================
  !
  
  !>Destroys the given CellML environment.
  SUBROUTINE CellML_Destroy(cellML,err,error,*)
    
    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML environment to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !< The error string
    !Local variables
    INTEGER(INTG) :: cellMLIdx,cellMLPosition
    TYPE(CellMLEnvironmentsType), POINTER :: cellMLEnvironments
    TYPE(CellMLPtrType), ALLOCATABLE :: newEnvironments(:)
    
    ENTERS("CellML_Destroy",err,error,*999)

    IF(.NOT.ASSOCIATED(cellML)) CALL FlagError("CellML is not associated.",err,error,*999)

    NULLIFY(cellMLEnvironments)
    CALL CellML_CellMLEnvironmentsGet(cellML,cellMLEnvironments,err,error,*999)    
    cellMLPosition=cellML%globalNumber
    CALL CellML_Finalise(cellML,err,error,*999)
    IF(cellMLEnvironments%numberOfEnvironments>1) THEN
      ALLOCATE(newEnvironments(cellMLEnvironments%numberOfEnvironments-1),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocated new CellML environments.",err,error,*999)
      DO cellMLIdx=1,cellMLEnvironments%numberOfEnvironments
        IF(cellMLIdx<cellMLPosition) THEN
          newEnvironments(cellMLIdx)%ptr=>cellMLEnvironments%environments(cellMLIdx)%ptr
        ELSE IF(cellMLIdx>cellMLPosition) THEN
          cellMLEnvironments%environments(cellMLIdx)%ptr%globalNumber=cellMLEnvironments%environments(cellMLIdx)%ptr%globalNumber-1
          newEnvironments(cellMLIdx-1)%ptr=>cellMLEnvironments%environments(cellMLIdx)%ptr              
        ENDIF
      ENDDO !cellMLIdx
      CALL MOVE_ALLOC(newEnvironments,cellMLEnvironments%environments)
      cellMLEnvironments%numberOfEnvironments=cellMLEnvironments%numberOfEnvironments-1
    ELSE
      IF(ALLOCATED(cellMLEnvironments%environments)) DEALLOCATE(cellMLEnvironments%environments)
      cellMLEnvironments%numberOfEnvironments=0
    ENDIF

    EXITS("CellML_Destroy")
    RETURN
999 IF(ALLOCATED(newEnvironments)) DEALLOCATE(newEnvironments)
    ERRORSEXITS("CellML_Destroy",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_Destroy

  !
  !=================================================================================================================================
  !

  !>Updates any cellml fields from the mapped fields
  SUBROUTINE CellML_FieldToCellMLUpdate(cellML,err,error,*)
    
    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML environment to update
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !< The error string
    !Local variables
    INTEGER(INTG) :: derivativeNumber,dofIdx,dofType,dof2ParamIdx,elementNumber,gaussNumber,gridNumber,mapIdx,modelIdx, &
      & nodeNumber,versionNumber
    INTEGER(INTG), POINTER :: modelsData(:)
    REAL(DP) :: dofValue
    TYPE(CellMLFieldMapsType), POINTER :: cellMLFieldMaps
    TYPE(CellMLIntermediateFieldType), POINTER :: cellMLIntermediateField
    TYPE(CellMLModelMapType), POINTER :: cellMLModelMap
    TYPE(CellMLModelsFieldType), POINTER :: cellMLModelsField
    TYPE(CellMLParametersFieldType), POINTER :: cellMLParametersField
    TYPE(CellMLStateFieldType), POINTER :: cellMLStateField
    TYPE(CellMLModelMapsType), POINTER :: cellMLModelMaps
    TYPE(FieldType), POINTER :: modelsField
    TYPE(FieldVariableType), POINTER :: modelsVariable
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("CellML_FieldToCellMLUpdate",err,error,*999)

    NULLIFY(cellMLModelsField)
    CALL CellML_CellMLModelsFieldGet(cellML,cellMLModelsField,err,error,*999)
    NULLIFY(cellMLFieldMaps)
    CALL CellML_CellMLFieldMapsGet(cellML,cellMLFieldMaps,err,error,*999)
    
    IF(cellML%modelsField%onlyOneModelIndex/=0) THEN
      !We have CellML models on this rank
      IF(cellML%modelsField%onlyOneModelIndex/=CELLML_MODELS_FIELD_NOT_CONSTANT) THEN
        !The CellML environement only uses one model and so we can optimise for this.
        NULLIFY(cellMLModelMaps)
        CALL CellMLFieldMaps_CellMLModelMapsGet(cellMLFieldMaps,cellMLModelsField%onlyOneModelIndex,cellMLModelMaps,err,error,*999)
        IF(diagnostics1) THEN
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Field to CellML update:",err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  CellML user number = ",cellML%userNumber,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  One model index = ",cellMLModelsField%onlyOneModelIndex,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of model maps = ",cellMLModelMaps%numberOfFieldsMappedTo, &
            & err,error,*999)
        ENDIF
        !Loop over the number of field to CellML maps
        DO mapIdx=1,cellMLModelMaps%numberOfFieldsMappedTo
          NULLIFY(cellMLModelMap)
          CALL CellMLModelMaps_CellMLModelToMapGet(cellMLModelMaps,mapIdx,cellMLModelMap,err,error,*999)
          SELECT CASE(cellMLModelMap%cellMLFieldType)
          CASE(CELLML_MODELS_FIELD)
            CALL FlagError("Cannot map models field.",err,error,*999)
          CASE(CELLML_STATE_FIELD)
            NULLIFY(cellMLStateField)
            CALL CellML_CellMLStateFieldGet(cellML,cellMLStateField,err,error,*999)
            CALL Field_ParametersToFieldParametersCopy(cellMLModelMap%field,cellMLModelMap%variableType, &
              & cellMLModelMap%fieldParameterSet,cellMLModelMap%componentNumber,cellMLStateField%stateField, &
              & FIELD_U_VARIABLE_TYPE,cellMLModelMap%cellMLParameterSet,cellMLModelMap%cellMLVariableNumber,err,error,*999)
          CASE(CELLML_INTERMEDIATE_FIELD)
            NULLIFY(cellMLIntermediateField)
            CALL CellML_CellMLIntermediateFieldGet(cellML,cellMLIntermediateField,err,error,*999)
            CALL Field_ParametersToFieldParametersCopy(cellMLModelMap%field,cellMLModelMap%variableType, &
              & cellMLModelMap%fieldParameterSet,cellMLModelMap%componentNumber,cellMLIntermediateField%intermediateField, &
              & FIELD_U_VARIABLE_TYPE,cellMLModelMap%cellMLParameterSet,cellMLModelMap%cellMLVariableNumber,err,error,*999)
          CASE(CELLML_PARAMETERS_FIELD)
            NULLIFY(cellMLParametersField)
            CALL CellML_CellMLParametersFieldGet(cellML,cellMLParametersField,err,error,*999)
            CALL Field_ParametersToFieldParametersCopy(cellMLModelMap%field,cellMLModelMap%variableType, &
              & cellMLModelMap%fieldParameterSet,cellMLModelMap%componentNumber,cellMLParametersField%parametersField, &
              & FIELD_U_VARIABLE_TYPE,cellMLModelMap%cellMLParameterSet,cellMLModelMap%cellMLVariableNumber,err,error,*999)
          CASE DEFAULT
            localError="The CellML to field model map CellML field type of "// &
              & TRIM(NumberToVString(cellMLModelMap%cellMLFieldType,"*",err,error))// &
              & " is invalid for map index "//TRIM(NumberToVString(mapIdx,"*",err,error))//" of model index "// &
              & TRIM(NumberToVString(cellMLModelsField%onlyOneModelIndex,"*",err,error))// &
              & " of CellML environment number "//TRIM(NumberToVString(cellML%userNumber,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          IF(diagnostics1) THEN
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Map index : ",mapIdx,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Field user number      = ",cellMLModelMap%field%userNumber, &
              & err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Field variable type    = ",cellMLModelMap%variableType, &
              & err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Field parameter set    = ",cellMLModelMap%fieldParameterSet, &
              & err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Field component number = ",cellMLModelMap%componentNumber, &
              & err,error,*999)                   
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    CellML field type      = ",cellMLModelMap%cellMLFieldType, &
              & err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    CellML parameter set   = ",cellMLModelMap%cellMLParameterSet, &
              & err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    CellML variable number = ",cellMLModelMap%cellMLVariableNumber, &
              & err,error,*999)
          ENDIF
        ENDDO !mapIdx
      ELSE
        !More than one model is used in the models field
        NULLIFY(modelsField)
        CALL CellMLModelsField_ModelsFieldGet(cellMLModelsField,modelsField,err,error,*999)
        NULLIFY(modelsVariable)
        CALL Field_VariableGet(modelsField,FIELD_U_VARIABLE_TYPE,modelsVariable,err,error,*999)
        NULLIFY(modelsData)
        CALL FieldVariable_ParameterSetDataGet(modelsVariable,FIELD_VALUES_SET_TYPE,modelsData,err,error,*999)
        IF(diagnostics1) THEN
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Field to CellML update:",err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  CellML user number = ",cellML%userNumber,err,error,*999)
        ENDIF
        !Loop over the dofs in the models field
        DO dofIdx=1,modelsVariable%totalNumberOfDofs
          modelIdx=modelsData(dofIdx)
          IF(modelIdx>0) THEN
            NULLIFY(cellMLModelMaps)
            CALL CellMLFieldMaps_CellMLModelMapsGet(cellMLFieldMaps,modelIdx,cellMLModelMaps,err,error,*999)
            dofType=modelsVariable%dofToParamMap%DOFType(1,dofIdx)
            dof2ParamIdx=modelsVariable%dofToParamMap%DOFType(2,dofIdx)
            SELECT CASE(dofType)
            CASE(FIELD_CONSTANT_INTERPOLATION)
              !Do nothing
            CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
              elementNumber=modelsVariable%dofToParamMap%elementDOF2ParamMap(1,dof2ParamIdx)
            CASE(FIELD_NODE_BASED_INTERPOLATION)
              versionNumber=modelsVariable%dofToParamMap%nodeDOF2ParamMap(1,dof2ParamIdx)
              derivativeNumber=modelsVariable%dofToParamMap%nodeDOF2ParamMap(2,dof2ParamIdx)
              nodeNumber=modelsVariable%dofToParamMap%nodeDOF2ParamMap(3,dof2ParamIdx)
            CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
              gridNumber=modelsVariable%dofToParamMap%gridPointDOF2ParamMap(1,dof2ParamIdx)
            CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
              gaussNumber=modelsVariable%dofToParamMap%gaussPointDOF2ParamMap(1,dof2ParamIdx)
              elementNumber=modelsVariable%dofToParamMap%gaussPointDOF2ParamMap(2,dof2ParamIdx)
            CASE DEFAULT
              localError="The DOF type of "//TRIM(NumberToVString(dofType,"*",err,error))// &
                & " for local DOF number "//TRIM(NumberToVString(dofIdx,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
            IF(diagnostics1) THEN
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  DOF index : ",dofIdx,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Model index : ",modelIdx,err,error,*999)
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      DOF type = ",dofType,err,error,*999)
              SELECT CASE(dofType)
              CASE(FIELD_CONSTANT_INTERPOLATION)
                !Do nothing
              CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Element number = ",elementNumber,err,error,*999)
              CASE(FIELD_NODE_BASED_INTERPOLATION)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Version number = ",versionNumber,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Derivative number = ",derivativeNumber, &
                  & err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Node number = ",nodeNumber,err,error,*999)
              CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Grid number = ",gridNumber,err,error,*999)
              CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Gauss number = ",gaussNumber,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Element number = ",elementNumber,err,error,*999)
              CASE DEFAULT
                localError="The DOF type of "//TRIM(NumberToVString(dofType,"*",err,error))// &
                  & " for local DOF number "//TRIM(NumberToVString(dofIdx,"*",err,error))//" is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of model maps = ",cellMLModelMaps%numberOfFieldsMappedTo, &
                & err,error,*999)
            ENDIF
            !Loop over the number of field to CellML maps
            DO mapIdx=1,cellMLModelMaps%numberOfFieldsMappedTo
              NULLIFY(cellMLModelMap)
              CALL CellMLModelMaps_CellMLModelToMapGet(cellMLModelMaps,mapIdx,cellMLModelMap,err,error,*999)
              !Get the OpenCMISS mapped field DOF value
              SELECT CASE(dofType)
              CASE(FIELD_CONSTANT_INTERPOLATION)
                CALL Field_ParameterSetGetConstant(cellMLModelMap%field,cellMLModelMap%variableType, &
                  & cellMLModelMap%fieldParameterSet,cellMLModelMap%componentNumber,dofValue,err,error,*999)
              CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                CALL Field_ParameterSetGetLocalElement(cellMLModelMap%field,cellMLModelMap%variableType, &
                  & cellMLModelMap%fieldParameterSet,elementNumber,cellMLModelMap%componentNumber,dofValue, &
                  & err,error,*999)
              CASE(FIELD_NODE_BASED_INTERPOLATION)
                CALL Field_ParameterSetGetLocalNode(cellMLModelMap%field,cellMLModelMap%variableType, &
                  & cellMLModelMap%fieldParameterSet,versionNumber,derivativeNumber,nodeNumber, &
                  & cellMLModelMap%componentNumber,dofValue,err,error,*999)
              CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                CALL FlagError("Not implemented.",err,error,*999)
              CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                CALL Field_ParameterSetGetLocalGaussPoint(cellMLModelMap%field,cellMLModelMap%variableType, &
                  & cellMLModelMap%fieldParameterSet,gaussNumber,elementNumber, &
                  & cellMLModelMap%componentNumber,dofValue,err,error,*999)
              CASE DEFAULT
                localError="The DOF type of "//TRIM(NumberToVString(dofType,"*",err,error))// &
                  & " for local DOF number "//TRIM(NumberToVString(dofIdx,"*",err,error))//" is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
              !Update the CellML field DOF value
              SELECT CASE(cellMLModelMap%cellMLFieldType)
              CASE(CELLML_MODELS_FIELD)
                CALL FlagError("Cannot map models field.",err,error,*999)
              CASE(CELLML_STATE_FIELD)
                NULLIFY(cellMLStateField)
                CALL CellML_CellMLStateFieldGet(cellML,cellMLStateField,err,error,*999)
                SELECT CASE(dofType)
                CASE(FIELD_CONSTANT_INTERPOLATION)
                  CALL Field_ParameterSetUpdateConstant(cellMLStateField%stateField,FIELD_U_VARIABLE_TYPE, &
                    & cellMLModelMap%cellMLParameterSet,cellMLModelMap%cellMLVariableNumber,dofValue,err,error,*999)
                CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                  CALL Field_ParameterSetUpdateLocalElement(cellMLStateField%stateField,FIELD_U_VARIABLE_TYPE, &
                    & cellMLModelMap%cellMLParameterSet,elementNumber,cellMLModelMap%cellMLVariableNumber,dofValue, &
                    & err,error,*999)
                CASE(FIELD_NODE_BASED_INTERPOLATION)
                  CALL Field_ParameterSetUpdateLocalNode(cellMLStateField%stateField,FIELD_U_VARIABLE_TYPE, &
                    & cellMLModelMap%cellMLParameterSet,versionNumber,derivativeNumber,nodeNumber, &
                    & cellMLModelMap%cellMLVariableNumber,dofValue,err,error,*999)
                CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                  CALL Field_ParameterSetUpdateLocalGaussPoint(cellMLStateField%stateField,FIELD_U_VARIABLE_TYPE, &
                    & cellMLModelMap%cellMLParameterSet,gaussNumber,elementNumber,cellMLModelMap%cellMLVariableNumber, &
                    & dofValue,err,error,*999)
                CASE DEFAULT
                  localError="The DOF type of "//TRIM(NumberToVString(dofType,"*",err,error))// &
                    & " for local DOF number "//TRIM(NumberToVString(dofIdx,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE(CELLML_INTERMEDIATE_FIELD)
                NULLIFY(cellMLIntermediateField)
                CALL CellML_CellMLIntermediateFieldGet(cellML,cellMLIntermediateField,err,error,*999)
                SELECT CASE(dofType)
                CASE(FIELD_CONSTANT_INTERPOLATION)
                  CALL Field_ParameterSetUpdateConstant(cellMLIntermediateField%intermediateField, &
                    & FIELD_U_VARIABLE_TYPE,cellMLModelMap%cellMLParameterSet,cellMLModelMap%cellMLVariableNumber, &
                    & dofValue,err,error,*999)
                CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                  CALL Field_ParameterSetUpdateLocalElement(cellMLIntermediateField%intermediateField, &
                    & FIELD_U_VARIABLE_TYPE,cellMLModelMap%cellMLParameterSet,elementNumber, &
                    & cellMLModelMap%cellMLVariableNumber,dofValue,err,error,*999)
                CASE(FIELD_NODE_BASED_INTERPOLATION)
                  CALL Field_ParameterSetUpdateLocalNode(cellMLIntermediateField%intermediateField, &
                    & FIELD_U_VARIABLE_TYPE,cellMLModelMap%cellMLParameterSet,versionNumber,derivativeNumber, &
                    & nodeNumber,cellMLModelMap%cellMLVariableNumber,dofValue,err,error,*999)
                CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                  CALL Field_ParameterSetUpdateLocalGaussPoint(cellMLIntermediateField%intermediateField, &
                    & FIELD_U_VARIABLE_TYPE,cellMLModelMap%cellMLParameterSet,gaussNumber,elementNumber, &
                    & cellMLModelMap%cellMLVariableNumber,dofValue,err,error,*999)
                CASE DEFAULT
                  localError="The DOF type of "//TRIM(NumberToVString(dofType,"*",err,error))// &
                    & " for local DOF number "//TRIM(NumberToVString(dofIdx,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE(CELLML_PARAMETERS_FIELD)
                NULLIFY(cellMLParametersField)
                CALL CellML_CellMLParametersFieldGet(cellML,cellMLParametersField,err,error,*999)
                SELECT CASE(dofType)
                CASE(FIELD_CONSTANT_INTERPOLATION)
                  CALL Field_ParameterSetUpdateConstant(cellMLParametersField%parametersField, &
                    & FIELD_U_VARIABLE_TYPE,cellMLModelMap%cellMLParameterSet,cellMLModelMap%cellMLVariableNumber, &
                    & dofValue,err,error,*999)
                CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                  CALL Field_ParameterSetUpdateLocalElement(cellMLParametersField%parametersField, &
                    & FIELD_U_VARIABLE_TYPE,cellMLModelMap%cellMLParameterSet,elementNumber, &
                    & cellMLModelMap%cellMLVariableNumber,dofValue,err,error,*999)
                CASE(FIELD_NODE_BASED_INTERPOLATION)
                  CALL Field_ParameterSetUpdateLocalNode(cellMLParametersField%parametersField, &
                    & FIELD_U_VARIABLE_TYPE,cellMLModelMap%cellMLParameterSet,versionNumber,derivativeNumber, &
                    & nodeNumber,cellMLModelMap%cellMLVariableNumber,dofValue,err,error,*999)
                CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                  CALL Field_ParameterSetUpdateLocalGaussPoint(cellML%ParametersField%parametersField, &
                    & FIELD_U_VARIABLE_TYPE,cellMLModelMap%cellMLParameterSet,gaussNumber,elementNumber, &
                    & cellMLModelMap%cellMLVariableNumber,dofValue,err,error,*999)
                CASE DEFAULT
                  localError="The DOF type of "//TRIM(NumberToVString(dofType,"*",err,error))// &
                    & " for local DOF number "//TRIM(NumberToVString(dofIdx,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE DEFAULT
                localError="The CellML to field model map CellML field type of "// &
                  & TRIM(NumberToVString(cellMLModelMap%cellMLFieldType,"*",err,error))// &
                  & " is invalid for map index "//TRIM(NumberToVString(mapIdx,"*",err,error))// &
                  & " of model index "//TRIM(NumberToVString(cellMLModelsField%onlyOneModelIndex,"*",err,error))// &
                  & " of CellML environment number "//TRIM(NumberToVString(cellML%userNumber,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              END SELECT
              IF(diagnostics1) THEN
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Map index : ",mapIdx,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        CellML field type      = ", &
                  & cellMLModelMap%cellMLFieldType,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        CellML parameter set   = ", &
                  & cellMLModelMap%cellMLParameterSet,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        CellML variable number = ", &
                  & cellMLModelMap%cellMLVariableNumber,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Field user number      = ", &
                  & cellMLModelMap%field%userNumber,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Field variable type    = ", &
                  & cellMLModelMap%variableType,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Field parameter set    = ", &
                  & cellMLModelMap%fieldParameterSet,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Field component number = ", &
                  & cellMLModelMap%componentNumber,err,error,*999)                   
              ENDIF
            ENDDO !mapIdx
          ENDIF !modelIdx>0
        ENDDO !dofIdx              
        CALL Field_ParameterSetDataRestore(modelsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & modelsData,err,error,*999)
      ENDIF
    ENDIF

    EXITS("CellML_FieldToCellMLUpdate")
    RETURN
999 ERRORSEXITS("CellML_FieldToCellMLUpdate",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_FieldToCellMLUpdate

  !
  !=================================================================================================================================
  !

  !>Finalise a CellML environment and deallocate all memory.
  SUBROUTINE CellML_Finalise(cellML,err,error,*)
    
    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML environment to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !< The error string
    !Local variables
    INTEGER(INTG) :: modelIdx
    
    ENTERS("CellML_Finalise",err,error,*999)

    IF(ASSOCIATED(cellML)) THEN
      IF(ALLOCATED(cellML%models)) THEN
        DO modelIdx=1,SIZE(cellML%models,1)
          CALL CellML_ModelFinalise(cellML%models(modelIdx)%ptr,err,error,*999)
        ENDDO !modelIdx
        DEALLOCATE(cellML%models)
      ENDIF
      CALL CellML_FieldMapsFinalise(cellML%fieldMaps,err,error,*999)
      CALL CellML_ModelsFieldFinalise(cellML%modelsField,err,error,*999)
      CALL CellML_StateFieldFinalise(cellML%stateField,err,error,*999)
      CALL CellML_IntermediateFieldFinalise(cellML%intermediateField,err,error,*999)
      CALL CellML_ParametersFieldFinalise(cellML%parametersField,err,error,*999)
      DEALLOCATE(cellML)
    ENDIF

    EXITS("CellML_Finalise")
    RETURN
999 ERRORSEXITS("CellML_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_Finalise

  !
  !=================================================================================================================================
  !

  !>Initialise a CellML environment and deallocate all memory.
  SUBROUTINE CellML_Initialise(cellML,err,error,*)
    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML environment to initialise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !< The error string
    !Local variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("CellML_Initialise",err,error,*998)

#ifdef WITH_CELLML

    IF(ASSOCIATED(cellML)) CALL FlagError("CellML environment is already associated.",err,error,*998)
    
    ALLOCATE(cellML,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate CellML environment.",err,error,*999)
    cellML%globalNumber=0
    cellML%userNumber=0
    NULLIFY(cellML%region)
    NULLIFY(cellML%environments)
    cellML%cellMLFinished=.FALSE.
    cellML%numberOfModels=0
    cellML%maximumNumberOfState=0
    cellML%maximumNumberOfParameters=0
    cellML%maximumNumberOfIntermediate=0
    NULLIFY(cellML%fieldMaps)
    NULLIFY(cellML%modelsField)
    NULLIFY(cellML%stateField)
    NULLIFY(cellML%intermediateField)
    NULLIFY(cellML%parametersField)
    cellML%cellMLGenerated=.FALSE.

#else

    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*998)

#endif

    EXITS("CellML_Initialise")
    RETURN
999 CALL CellML_Finalise(cellML,dummyErr,dummyError,*998)
998 ERRORSEXITS("CellML_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_Initialise

  !
  !=================================================================================================================================
  !

  !>Finish creating the field maps for a CellML environment.
  SUBROUTINE CellML_FieldMapsCreateFinish(cellML,err,error,*)
    
    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<The CellML environment to finish creating the maps for. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: modelIdx
    TYPE(CellMLFieldMapsType), POINTER :: cellMLFieldMaps
    TYPE(CellMLModelMapsType), POINTER :: cellMLModelMaps
    TYPE(VARYING_STRING) :: localError
     
    ENTERS("CellML_FieldMapsCreateFinish",err,error,*999)

#ifdef WITH_CELLML

    CALL CellML_AssertIsFinished(cellML,err,error,*999)

    NULLIFY(cellMLFieldMaps)
    CALL CellML_CellMLFieldMapsGet(cellML,cellMLFieldMaps,err,error,*999)
    CALL CellMLFieldMaps_AssertNotFinished(cellMLFieldMaps,err,error,*999)
    !Check that each model has field mappings
    DO modelIdx=1,cellML%numberOfModels
      NULLIFY(cellMLModelMaps)
      CALL CellMLFieldMaps_CellMLModelMapsGet(cellMLFieldMaps,modelIdx,cellMLModelMaps,err,error,*999)
      IF(cellMLModelMaps%numberOfFieldsMappedTo==0.AND.cellMLModelMaps%numberOfFieldsMappedFrom==0) THEN
        localError="Invalid setup. CellML model index "//TRIM(NumberToVString(modelIdx,"*",err,error))// &
          & " does not have any mappings to or from a field."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDDO !modelIdx
    cellML%fieldMaps%cellMLFieldMapsFinished=.TRUE.
    
#else

    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*999)

#endif

    EXITS("CellML_FieldMapsCreateFinish")
    RETURN
999 ERRORSEXITS("CellML_FieldMapsCreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_FieldMapsCreateFinish

  !
  !=================================================================================================================================
  !

  !>Start the creation of field maps for a CellML environment.
  SUBROUTINE CellML_FieldMapsCreateStart(cellML,err,error,*)
    
    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<The CellML environment to create the maps for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
  
    ENTERS("CellML_FieldMapsCreateStart",err,error,*999)

#ifdef WITH_CELLML

    CALL CellML_AssertIsFinished(cellML,err,error,*999)
    
    !Initialise the field maps
    CALL CellML_FieldMapsInitialise(cellML,err,error,*999)

#else

    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*999)

#endif

    EXITS("CellML_FieldMapsCreateStart")
    RETURN
999 ERRORSEXITS("CellML_FieldMapsCreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_FieldMapsCreateStart

  !
  !=================================================================================================================================
  !

  !>Finalise a CellML maps and deallocate all memory.
  SUBROUTINE CellML_FieldMapsFinalise(cellMLFieldMaps,err,error,*)
    
    !Argument variables
    TYPE(CellMLFieldMapsType), POINTER :: cellMLFieldMaps !<A pointer to the CellML field maps to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !< The error string
    !Local variables
    INTEGER(INTG) :: modelIdx
    
    ENTERS("CellML_FieldMapsFinalise",err,error,*999)

#ifdef WITH_CELLML

    IF(ASSOCIATED(cellMLFieldMaps)) THEN
      IF(ALLOCATED(cellMLFieldMaps%modelMaps)) THEN
        DO modelIdx=1,SIZE(cellMLFieldMaps%modelMaps,1)
          CALL CellML_ModelMapsFinalise(cellMLFieldMaps%modelMaps(modelIdx)%ptr,err,error,*999)
        ENDDO !modelIdx
        DEALLOCATE(cellMLFieldMaps%modelMaps)
      ENDIF
      DEALLOCATE(cellMLFieldMaps)
    ENDIF

#else

    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*999)

#endif

    EXITS("CellML_FieldMapsFinalise")
    RETURN
999 ERRORSEXITS("CellML_FieldMapsFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_FieldMapsFinalise

  !
  !=================================================================================================================================
  !

  !>Initialise a CellML field mpas.
  SUBROUTINE CellML_FieldMapsInitialise(cellML,err,error,*)
    
    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML environment to initialise the field maps for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !< The error string
    !Local variables
    INTEGER(INTG) :: dummyErr,modelIdx
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("CellML_FieldMapsInitialise",err,error,*998)

#ifdef WITH_CELLML

    CALL CellML_AssertIsFinished(cellML,err,error,*999)
    
    IF(ASSOCIATED(cellML%fieldMaps)) CALL FlagError("CellML environment field maps is already associated.",err,error,*999)
    
    ALLOCATE(cellML%fieldMaps,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate CellML field maps.",err,error,*999)
    cellML%fieldMaps%cellML=>cellML
    cellML%fieldMaps%cellMLFieldMapsFinished=.FALSE.
    NULLIFY(cellML%fieldMaps%sourceGeometricField)
    NULLIFY(cellML%fieldMaps%sourceFieldVariable)
    NULLIFY(cellML%fieldMaps%sourceFieldDomain)
    cellML%fieldMaps%sourceFieldInterpolationType=0
    ALLOCATE(cellML%fieldMaps%modelMaps(cellML%numberOfModels),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate CellML environment field maps model maps.",err,error,*999)
    DO modelIdx=1,cellML%numberOfModels
      NULLIFY(cellML%fieldMaps%modelMaps(modelIdx)%ptr)
      CALL CellML_ModelMapsInitialise(cellML%fieldMaps,cellML%fieldMaps%modelMaps(modelIdx)%ptr,err,error,*999)
    ENDDO !modelIdx
    
#else

    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*998)

#endif

    EXITS("CellML_FieldMapsInitialise")
    RETURN
999 CALL CellML_FieldMapsFinalise(cellML%fieldMaps,dummyErr,dummyError,*998)
998 ERRORSEXITS("CellML_FieldMapsInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_FieldMapsInitialise

  !
  !=================================================================================================================================
  !

  !>Finalise a CellML model and deallocate all memory.
  SUBROUTINE CellML_ModelFinalise(cellMLModel,err,error,*)
    
    !Argument variables
    TYPE(CellMLModelType), POINTER :: cellMLModel !<A pointer to the CellML model to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !< The error string
    !Local variables
     
    ENTERS("CellML_ModelFinalise",err,error,*999)

#ifdef WITH_CELLML

    IF(ASSOCIATED(cellMLModel)) THEN
      IF(C_ASSOCIATED(cellMLModel%ptr)) CALL DESTROY_CELLML_MODEL_DEFINITION(cellMLModel%ptr)
      cellMLModel%modelId=""
      DEALLOCATE(cellMLModel)
    ENDIF

#else

    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*999)

#endif

    EXITS("CellML_ModelFinalise")
    RETURN
999 ERRORSEXITS("CellML_ModelFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_ModelFinalise

  !
  !=================================================================================================================================
  !

  !>Initialise a CellML model.
  SUBROUTINE CellML_ModelInitialise(cellMLModel,err,error,*)
    
    !Argument variables
    TYPE(CellMLModelType), POINTER :: cellMLModel !<On return, a pointer to the CellML initialised model. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !< The error string
    !Local variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("CellML_ModelInitialise",err,error,*998)

#ifdef WITH_CELLML

    IF(ASSOCIATED(cellMLModel)) CALL FlagError("CellML model is already associated.",err,error,*998)
   
    ALLOCATE(cellMLModel,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate CellML model.",err,error,*999)
    NULLIFY(cellMLModel%cellML)
    cellMLModel%globalNumber=0
    cellMLModel%modelId=""
    cellMLModel%ptr=C_NULL_PTR
    cellMLModel%numberOfState=0
    cellMLModel%numberOfIntermediate=0
    cellMLModel%numberOfParameters=0
 
#else

    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*998)

#endif

    EXITS("CellML_ModelInitialise")
    RETURN
999 CALL CellML_ModelFinalise(cellMLModel,dummyErr,dummyError,*998)
998 ERRORSEXITS("CellML_ModelInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_ModelInitialise

  !
  !=================================================================================================================================
  !

  !>Finalise a CellML model map and deallocate all memory.
  SUBROUTINE CellML_ModelMapFinalise(cellMLModelMap,err,error,*)
    !Argument variables
    TYPE(CellMLModelMapType), POINTER :: cellMLModelMap !<A pointer to the CellML model map to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !< The error string
    !Local variables
    
    ENTERS("CellML_ModelMapFinalise",err,error,*999)

#ifdef WITH_CELLML

    IF(ASSOCIATED(cellMLModelMap)) THEN
      CALL ERASE(cellMLModelMap%variableId)
      DEALLOCATE(cellMLModelMap)
    ENDIF

#else

    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*999)

#endif

    EXITS("CellML_ModelMapFinalise")
    RETURN
999 ERRORSEXITS("CellML_ModelMapFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_ModelMapFinalise

  !
  !=================================================================================================================================
  !

  !>Initialise a CellML model map.
  SUBROUTINE CellML_ModelMapInitialise(cellMLModelMap,err,error,*)
    
    !Argument variables
    TYPE(CellMLModelMapType), POINTER :: cellMLModelMap !<On return, a pointer to the CellML initialised model map field. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !< The error string
    !Local variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("CellML_ModelMapInitialise",err,error,*998)

#ifdef WITH_CELLML

    IF(ASSOCIATED(cellMLModelMap)) CALL FlagError("CellML model map field is already associated.",err,error,*998)
    
    ALLOCATE(cellMLModelMap,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate CellML model map.",err,error,*999)
    cellMLModelMap%cellMLMapType = 0
    NULLIFY(cellMLModelMap%field)
    cellMLModelMap%variableType=0
    cellMLModelMap%componentNumber=0
    cellMLModelMap%fieldParameterSet=0
    cellMLModelMap%variableId=""
    cellMLModelMap%cellMLFieldType=0
    cellMLModelMap%cellMLVariableNumber=0
    cellMLModelMap%cellMLParameterSet=0
  
#else

    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*998)

#endif

    EXITS("CellML_ModelMapInitialise")
    RETURN
999 CALL CellML_ModelMapFinalise(cellMLModelMap,dummyErr,dummyError,*998)
998 ERRORSEXITS("CellML_ModelMapInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_ModelMapInitialise

  !
  !=================================================================================================================================
  !

  !>Finalise a CellML model maps and deallocate all memory.
  SUBROUTINE CellML_ModelMapsFinalise(cellMLModelMaps,err,error,*)
    
    !Argument variables
    TYPE(CellMLModelMapsType), POINTER :: cellMLModelMaps !<A pointer to the CellML model maps to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !< The error string
    !Local variables
    INTEGER(INTG) :: mapIdx
    
    ENTERS("CellML_ModelMapsFinalise",err,error,*999)

#ifdef WITH_CELLML

    IF(ASSOCIATED(cellMLModelMaps)) THEN
      IF(ALLOCATED(cellMLModelMaps%fieldsMappedTo)) THEN
        DO mapIdx=1,SIZE(cellMLModelMaps%fieldsMappedTo,1)
          CALL CellML_ModelMapFinalise(cellMLModelMaps%fieldsMappedTo(mapIdx)%ptr,err,error,*999)
        ENDDO !mapIdx
        DEALLOCATE(cellMLModelMaps%fieldsMappedTo)
      ENDIF
      IF(ALLOCATED(cellMLModelMaps%fieldsMappedFrom)) THEN
        DO mapIdx=1,SIZE(cellMLModelMaps%fieldsMappedFrom,1)
          CALL CellML_ModelMapFinalise(cellMLModelMaps%fieldsMappedFrom(mapIdx)%ptr,err,error,*999)
        ENDDO !mapIdx
        DEALLOCATE(cellMLModelMaps%fieldsMappedFrom)
      ENDIF
      DEALLOCATE(cellMLModelMaps)
    ENDIF

#else

    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*999)

#endif

    EXITS("CellML_ModelMapsFinalise")
    RETURN
999 ERRORSEXITS("CellML_ModelMapsFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_ModelMapsFinalise

  !
  !=================================================================================================================================
  !

  !>Initialise a CellML model maps.
  SUBROUTINE CellML_ModelMapsInitialise(cellMLFieldMaps,cellMLModelMaps,err,error,*)
    
    !Argument variables
    TYPE(CellMLFieldMapsType), POINTER :: cellMLFieldMaps !<A pointer to the field maps to initialise the model maps for
    TYPE(CellMLModelMapsType), POINTER :: cellMLModelMaps !<A pointer to the model maps to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !< The error string
    !Local variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("CellML_ModelMapsInitialise",err,error,*998)

#ifdef WITH_CELLML

    IF(ASSOCIATED(cellMLModelMaps)) CALL FlagError("CellML model maps is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cellMLFieldMaps)) CALL FlagError("CellML field maps is not associated.",err,error,*999)
    
    ALLOCATE(cellMLModelMaps,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate CellML model maps.",err,error,*999)
    cellMLModelMaps%cellMLFieldMaps=>cellMLFieldMaps
    cellMLModelMaps%numberOfFieldsMappedTo=0
    cellMLModelMaps%numberOfFieldsMappedFrom=0    

#else

    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*998)

#endif

    EXITS("CellML_ModelMapsInitialise")
    RETURN
999 CALL CellML_ModelMapsFinalise(cellMlModelMaps,dummyErr,dummyError,*998)
998 ERRORSEXITS("CellML_ModelMapsInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_ModelMapsInitialise

  !
  !=================================================================================================================================
  !

  !>Import the specified CellML model into the given CellML environment object.
  SUBROUTINE CellML_ModelImportC(cellML,URI,modelIndex,err,error,*)
    
    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<The CellML environment object into which we want to import the specified model.
    CHARACTER(LEN=*) :: URI !<The (absolute? relative?) URI of the model to import. As per tracker item 2013 comment 8 the uri should now simply point to a CellML document. Can use a relative URL which will be interpreted relative to the CWD of the executed application.
    INTEGER(INTG), INTENT(OUT) :: modelIndex !<On return, the index for this model within the given CellML environment object.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local variables
    INTEGER(INTG) :: cellMLIdx
    TYPE(CellMLModelType), POINTER :: newModel
    TYPE(CellMLModelPtrType), ALLOCATABLE :: newModels(:)
    CHARACTER(256) :: cURI
    INTEGER(INTG) :: cURIL

    ENTERS("CellML_ModelImportC",err,error,*999)

#ifdef WITH_CELLML

    NULLIFY(newModel)
    modelIndex=0

    IF(.NOT.ASSOCIATED(cellML)) CALL FlagError("CellML environment is not associated.",err,error,*999)
    
    !set up new model object
    !Allocate new CellML model
    NULLIFY(newModel)
    CALL CellML_ModelInitialise(newModel,err,error,*999)
    cURIL = LEN_TRIM(URI)
    WRITE(cURI,'(A,A)') URI(1:cURIL),C_NULL_CHAR
    newModel%ptr = CREATE_CELLML_MODEL_DEFINITION(cURI)
    IF(.NOT.C_ASSOCIATED(newModel%ptr)) CALL FlagError("Error instantiating CellML model.",err,error,*999)
    newModel%globalNumber=cellML%numberOfModels+1
    newModel%modelId=URI(1:cURIL)
    !Add model to this CellML environment's list of models
    ALLOCATE(newModels(cellML%numberOfModels+1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new CellML models array.",err,error,*999)
    DO cellMLIdx=1,cellML%numberOfModels
      newModels(cellMLIdx)%ptr=>cellML%models(cellMLIdx)%ptr
    ENDDO !cellMLIdx
    newModels(cellML%numberOfModels+1)%ptr=>newModel
    CALL MOVE_ALLOC(newModels,cellML%models)
    cellML%numberOfModels=cellML%numberOfModels+1
    modelIndex=cellML%numberOfModels

#else

    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*999)

#endif

    EXITS("CellML_ModelImportC")
    RETURN
999 ERRORSEXITS("CellML_ModelImportC",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_ModelImportC

  !
  !=================================================================================================================================
  !

  !>Import the specified CellML model into the given CellML environment object.
  SUBROUTINE CellML_ModelImportVS(cellML,URI,modelIndex,err,error,*)
    
    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<The CellML environment object into which we want to import the specified model.
    TYPE(VARYING_STRING), INTENT(IN) :: URI !<The (absolute? relative?) URI of the model to import. As per tracker item 2013 comment 8 the URI should now simply point to a CellML document. Can use a relative URL which will be interpreted relative to the CWD of the executed application.
    INTEGER(INTG), INTENT(OUT) :: modelIndex !<On return, the index for this model within the given CellML environment object.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local variables

    ENTERS("CellML_ModelImportVS",err,error,*999)

#ifdef WITH_CELLML

    CALL CellML_ModelImport(cellML,CHAR(URI),modelIndex,err,error,*999)

#else

    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*999)

#endif

    EXITS("CellML_ModelImportVS")
    RETURN
999 ERRORSEXITS("CellML_ModelImportVS",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_ModelImportVS

  !
  !=================================================================================================================================
  !

  !>Sets a CellML model variable to be known - i.e., the variable's value will be set by an OpenCMISS field
  SUBROUTINE CellML_VariableSetAsKnownC(cellML,modelIndex,variableID,err,error,*)

    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<The CellML environment object in which to create the map.
    INTEGER(INTG), INTENT(IN) :: modelIndex !<The index of the CellML model in which to find the given variable.
    CHARACTER(LEN=*), INTENT(IN) :: variableID !<The CellML variable to set as known (in the format 'component_name/variable_name').
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local variables
    CHARACTER(LEN=MAXSTRLEN) :: cName
    INTEGER(INTG) :: cNameL
    INTEGER(C_INT) :: errorCode
    TYPE(CellMLModelType), POINTER :: cellMLModel
    TYPE(VARYING_STRING) :: localError

    ENTERS("CellML_VariableSetAsKnownC",err,error,*999)

#ifdef WITH_CELLML

    CALL CellML_AssertNotFinished(cellML,err,error,*999)

    NULLIFY(cellMLModel)
    CALL CellML_CellMLModelGet(cellML,modelIndex,cellMLModel,err,error,*999)
    
    !All input arguments are ok.
    cNameL = LEN_TRIM(variableID)
    WRITE(cName,'(A,A)') variableID(1:cNameL),C_NULL_CHAR
    errorCode = CELLML_MODEL_DEFINITION_SET_VARIABLE_AS_KNOWN(cellMLModel%ptr,cName)
    IF(errorCode/=0) THEN
      localError="Error "//TRIM(NumberToVString(errorCode,"*",err,error))// &
        & ". The specified variable can not be set as known: "//variableID
      CALL FlagError(localError,err,error,*999)
    ENDIF

#else

    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*999)

#endif

    EXITS("CellML_VariableSetAsKnownC")
    RETURN
999 ERRORSEXITS("CellML_VariableSetAsKnownC",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_VariableSetAsKnownC

  !
  !=================================================================================================================================
  !

  !>Sets a CellML model variable to be known - i.e., the variable's value will be set by an OpenCMISS field
  SUBROUTINE CellML_VariableSetAsKnownVS(cellML,modelUserNumber,variableID,err,error,*)

    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<The CellML environment object in which to create the map.
    INTEGER(INTG), INTENT(IN) :: modelUserNumber !<The index of the CellML model in which to find the given variable.
    TYPE(VARYING_STRING), INTENT(IN) :: variableID !<The CellML variable to set as known (in the format 'component_name/variable_name').
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local variables

    ENTERS("CellML_VariableSetAsKnownVS",err,error,*999)

#ifdef WITH_CELLML

    CALL CellML_VariableSetAsKnown(cellML,modelUserNumber,CHAR(variableID),err,error,*999)

#else

    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*999)

#endif

    EXITS("CellML_VariableSetAsKnownVS")
    RETURN
999 ERRORSEXITS("CellML_VariableSetAsKnownVS",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_VariableSetAsKnownVS

  !
  !=================================================================================================================================
  !

  !>Sets a CellML model variable to be wanted - i.e., the variable's value will used by an OpenCMISS field
  SUBROUTINE CellML_VariableSetAsWantedC(cellML,modelIndex,variableID,err,error,*)

    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<The CellML environment object in which to create the map.
    INTEGER(INTG), INTENT(IN) :: modelIndex !<The index of the CellML model in which to find the given variable.
    CHARACTER(LEN=*), INTENT(IN) :: variableID !<The CellML variable to set as wanted (in the format 'component_name/variable_name').
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local variables
    INTEGER(INTG) :: cNameL
    INTEGER(C_INT) :: errorCode
    CHARACTER(LEN=MAXSTRLEN) :: cName
    TYPE(CellMLModelType), POINTER :: cellMLModel
    TYPE(VARYING_STRING) :: localError

    ENTERS("CellML_VariableSetAsWantedC",err,error,*999)

#ifdef WITH_CELLML

    CALL CellML_AssertNotFinished(cellML,err,error,*999)

    NULLIFY(cellMLModel)
    CALL CellML_CellMLModelGet(cellML,modelIndex,cellMLModel,err,error,*999)
    !All input arguments are ok.
    cNameL = LEN_TRIM(variableID)
    WRITE(cName,'(A,A)') variableID(1:cNameL),C_NULL_CHAR
    errorCode = CELLML_MODEL_DEFINITION_SET_VARIABLE_AS_WANTED(cellMLModel%ptr,cName)
    IF(errorCode/=0) THEN
      localError="Error "//TRIM(NumberToVString(errorCode,"*",err,error))// &
        & ". The specified variable can not be set as wanted: "//variableID
      CALL FlagError(localError,err,error,*999)
    ENDIF
 
#else

    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*999)

#endif

    EXITS("CellML_VariableSetAsWantedC")
    RETURN
999 ERRORSEXITS("CellML_VariableSetAsWantedC",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_VariableSetAsWantedC

  !
  !=================================================================================================================================
  !

  !>Sets a CellML model variable to be wanted - i.e., the variable's value will be used by an OpenCMISS field
  SUBROUTINE CellML_VariableSetAsWantedVS(cellML,modelUserNumber,variableID,err,error,*)

    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<The CellML environment object in which to create the map.
    INTEGER(INTG), INTENT(IN) :: modelUserNumber !<The index of the CellML model in which to find the given variable.
    TYPE(VARYING_STRING), INTENT(IN) :: variableID !<The CellML variable to set as wanted (in the format 'component_name/variable_name').
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local variables

    ENTERS("CellML_VariableSetAsWantedVS",err,error,*999)

#ifdef WITH_CELLML

    CALL CellML_VariableSetAsWanted(cellML,modelUserNumber,CHAR(variableID),err,error,*999)

#else

    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*999)

#endif

    EXITS("CellML_VariableSetAsWantedVS")
    RETURN
999 ERRORSEXITS("CellML_VariableSetAsWantedVS",err,error)
    RETURN 1
  END SUBROUTINE CellML_VariableSetAsWantedVS

  !
  !=================================================================================================================================
  !

  !>Create a CellML model variable to field variable component map.
  SUBROUTINE CellML_CreateCellMLToFieldMapC(cellML,modelIndex,variableID,cellMLParameterSetType,field,variableType, &
    & componentNumber,fieldParamSetType,err,error,*)
    
    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<The CellML environment object in which to create the map.
    INTEGER(INTG), INTENT(IN) :: modelIndex !<The index of the CellML model to map from.
    CHARACTER(LEN=*), INTENT(IN) :: variableID !<The ID of the CellML variable in the given model to map from.
    INTEGER(INTG), INTENT(IN) :: cellMLParameterSetType !<The CellML parameter set to map from.
    TYPE(FieldType), POINTER, INTENT(IN) :: field !<The field to map to.
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable to map to.
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The field variable component number to map to.
    INTEGER(INTG), INTENT(IN) :: fieldParamSetType !<The field variable parameter set to map to.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local variables
    !INTEGER(INTG) :: cNameL
    INTEGER(C_INT) :: errorC
    INTEGER(INTG) :: cellMLVariableType,cellMLFieldType,cellMLVariableNumber,fieldInterpolationType,mapIdx
    CHARACTER(LEN=1,KIND=C_CHAR) :: cName(MAXSTRLEN)
    TYPE(CellMLFieldMapsType), POINTER :: cellMLFieldMaps
    TYPE(CellMLModelType), POINTER :: cellMLModel
    TYPE(CellMLModelMapType), POINTER :: newCellMLModelMap
    TYPE(CellMLModelMapPtrType), ALLOCATABLE :: newFieldsMappedFrom(:)
    TYPE(CellMLModelMapsType), POINTER :: cellMLModelMaps
    TYPE(DomainType), POINTER :: fieldDomain,sourceFieldDomain
    TYPE(FieldType), POINTER :: fieldGeometricField,sourceGeometricField
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(FieldParameterSetType), POINTER :: fieldParameterSet
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("CellML_CreateCellMLToFieldMapC",err,error,*999)

#ifdef WITH_CELLML

    CALL CellML_AssertIsFinished(cellML,err,error,*999)
    NULLIFY(fieldGeometricField)
    CALL Field_GeometricFieldGet(field,fieldGeometricField,err,error,*999)
    NULLIFY(fieldVariable)
    CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
    NULLIFY(fieldParameterSet)
    CALL FieldVariable_ParameterSetGet(fieldVariable,fieldParamSetType,fieldParameterSet,err,error,*999)
    NULLIFY(fieldDomain)
    CALL FieldVariable_ComponentDomainGet(fieldVariable,componentNumber,fieldDomain,err,error,*999)
    CALL FieldVariable_ComponentInterpolationGet(fieldVariable,componentNumber,fieldInterpolationType,err,error,*999)
    NULLIFY(cellMLModel)
    CALL CellML_CellMLModelGet(cellML,modelIndex,cellMLModel,err,error,*999)
    NULLIFY(cellMLFieldMaps)
    CALL CellML_CellMLFieldMapsGet(cellML,cellMLFieldMaps,err,error,*999)
    CALL CellMLFieldMaps_AssertNotFinished(cellMLfieldMaps,err,error,*999)
    NULLIFY(cellMLModelMaps)
    CALL CellMLFieldMaps_CellMLModelMapsGet(cellMLFieldMaps,modelIndex,cellMLModelMaps,err,error,*999)
    CALL CMISSF2CString(variableID,cName)
    !All input arguments are ok.Get the type of the variable being mapped
    errorC = CELLML_MODEL_DEFINITION_GET_VARIABLE_TYPE(cellMLModel%ptr,cName,cellMLVariableType)
    IF(errorC /= 0) THEN
      localError="Error "//TRIM(NumberToVString(errorC,"*",err,error))// &
        & ". Failed to get the type of CellML variable: "//variableID
      CALL FlagError(localError,err,error,*999)
    ENDIF
    cellMLFieldType=CellML_MapCellMLVariableToFieldType(cellMLVariableType,err,error)
    IF(err/=0) GOTO 999
    CALL CellML_FieldComponentGet(cellML,modelIndex,cellMLFieldType,variableID,cellMLVariableNumber,err,error,*999)
    !cNameL = LEN_TRIM(variableID)
    !WRITE(cName,'(A,A)') cName(1:cNameL),C_NULL_CHAR
    !cellMLVariableNumber=CELLML_MODEL_DEFINITION_ADD_MAPPING_TO_FIELD(cellMLModel%ptr,cName)
    !Now check that the mapped field is consistent with the other mapped fields for the model.
    sourceGeometricField=>cellMLFieldMaps%sourceGeometricField
    IF(ASSOCIATED(sourceGeometricField)) THEN
      IF(.NOT.ASSOCIATED(sourceGeometricField,fieldGeometricField)) THEN
        localError="The geometric field for field user number "// &
          & TRIM(NumberToVString(field%userNumber,"*",err,error))// &
          & " does not match the geometric field for other field variable components mapped" // &
          & " in the CellML environment."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      NULLIFY(sourceFieldDomain)
      CALL CellMLFieldMaps_SourceFieldDomainGet(cellMLFieldMaps,sourceFieldDomain,err,error,*999)
      IF(.NOT.ASSOCIATED(sourceFieldDomain,fieldDomain)) THEN
        localError="The domain for component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
          & " of variable type "//TRIM(NumberToVString(variableType,"*",err,error))// &
          & " of field user number "//TRIM(NumberToVString(field%userNumber,"*",err,error))// &
          & " does not match the domain for other field variable components mapped in the CellML environment."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      IF(cellMLFieldMaps%sourceFieldInterpolationType/=fieldInterpolationType) THEN
        localError="The interpolation type of "//TRIM(NumberToVString(fieldInterpolationType,"*",err,error))// &
          & " for component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
          & " of variable type "//TRIM(NumberToVString(variableType,"*",err,error))// &
          & " of field user number "//TRIM(NumberToVString(field%userNumber,"*",err,error))// &
          & " does not match the interpolation type of "// &
          & TRIM(NumberToVString(cellMLFieldMaps%sourceFieldInterpolationType,"*",err,error))// &
          & " used in other field variable components mapped in the CellML environment."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      cellMLFieldMaps%sourceGeometricField=>fieldGeometricField
      cellMLFieldMaps%sourceFieldVariable=>fieldVariable
      cellMLFieldMaps%sourceFieldDomain=>fieldDomain
      cellMLFieldMaps%sourceFieldInterpolationType=fieldInterpolationType
    ENDIF
    !Everything is OK so create the model map field. Get the type of the variable being mapped
    errorC = CELLML_MODEL_DEFINITION_GET_VARIABLE_TYPE(cellMLModel%ptr,cName,cellMLVariableType)
    IF(errorC /= 0) THEN
      localError="Error "//TRIM(NumberToVString(errorC,"*",err,error))//". Failed to get the TYPE of CellML variable: "// &
        & variableID
      CALL FlagError(localError,err,error,*999)
    ENDIF
    cellMLFieldType=CellML_MapCellMLVariableToFieldType(cellMLVariableType,err,error)
    IF(err/=0) GOTO 999
    CALL CellML_FieldComponentGet(cellML,modelIndex,cellMLFieldType,variableID,cellMLVariableNumber,err,error,*999)
    NULLIFY(newCellMLModelMap)
    CALL CellML_ModelMapInitialise(newCellMLModelMap,err,error,*999)
    newCellMLModelMap%cellMLMapType=CELLML_MAP_FROM_FIELD_TYPE
    newCellMLModelMap%field=>field
    newCellMLModelMap%variableType=variableType
    newCellMLModelMap%componentNumber=componentNumber
    newCellMLModelMap%fieldParameterSet=fieldParamSetType
    newCellMLModelMap%variableId=variableID
    newCellMLModelMap%cellMLFieldType=cellMLFieldType
    newCellMLModelMap%cellMLVariableNumber=cellMLVariableNumber
    newCellMLModelMap%cellMLParameterSet=cellMLParameterSetType
    !Put this model map field into the list of to field maps
    ALLOCATE(newFieldsMappedFrom(cellMLModelMaps%numberOfFieldsMappedFrom+1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new fields mapped from.",err,error,*999)
    DO mapIdx=1,cellMLModelMaps%numberOfFieldsMappedFrom
      newFieldsMappedFrom(mapIdx)%ptr=>cellMLModelMaps%fieldsMappedFrom(mapIdx)%ptr
    ENDDO !mapIdx
    newFieldsMappedFrom(cellMLModelMaps%numberOfFieldsMappedFrom+1)%ptr=>newCellMLModelMap
    CALL MOVE_ALLOC(newFieldsMappedFrom,cellMLModelMaps%fieldsMappedFrom)
    cellMLModelMaps%numberOfFieldsMappedFrom=cellMLModelMaps%numberOfFieldsMappedFrom+1          
    
    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"CellML model variable -> field map:",err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE," CellML model :",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"   CellML User number      = ",cellML%userNumber,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"   Model index             = ",modelIndex,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"   Variable ID             = ",variableID,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"   CellML field type       = ",cellMLFieldType,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"   CellML variable number  = ",cellMLVariableNumber,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"   CellML parameter set    = ",cellMLParameterSetType,err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE," Field :",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"   User number             = ",field%userNumber,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"   Variable type           = ",variableType,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"   Component number        = ",componentNumber,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"   Parameter set           = ",fieldParamSetType,err,error,*999)
    ENDIF

#else

    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*999)

#endif

    EXITS("CellML_CreateCellMLToFieldMapC")
    RETURN
999 ERRORSEXITS("CellML_CreateCellMLToFieldMapC",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_CreateCellMLToFieldMapC

  !
  !=================================================================================================================================
  !

  !>Create a CellML model variable to field variable component map.
  SUBROUTINE CellML_CreateCellMLToFieldMapVS(cellML,modelUserNumber,variableID,cellMLParameterSetType, &
    & field,variableType,componentNumber,fieldParamSetType,err,error,*)
    
    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<The CellML environment object in which to create the map.
    INTEGER(INTG), INTENT(IN) :: modelUserNumber !<The user number of the CellML model to map from.
    TYPE(VARYING_STRING), INTENT(IN) :: variableID !<The ID of the CellML variable in the given model to map from.
    INTEGER(INTG), INTENT(IN) :: cellMLParameterSetType !<The CellML parameter set to map from.
    TYPE(FieldType), POINTER, INTENT(IN) :: field !<The field to map to.
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable to map to.
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The field variable component number to map to.
    INTEGER(INTG), INTENT(IN) :: fieldParamSetType !<The field variable parameter set to map to.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local variables

    ENTERS("CellML_CreateCellMLToFieldMapVS",err,error,*999)

#ifdef WITH_CELLML

    CALL CellML_CreateCellMLToFieldMap(cellML,modelUserNumber,CHAR(variableID),cellMLParameterSetType, &
      & field,variableType,componentNumber,fieldParamSetType,err,error,*999)

#else

    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*999)

#endif

    EXITS("CellML_CreateCellMLToFieldMapVS")
    RETURN
999 ERRORSEXITS("CellML_CreateCellMLToFieldMapVS",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_CreateCellMLToFieldMapVS

  !
  !=================================================================================================================================
  !

  !>Create a field variable component to CellML model variable map.
  SUBROUTINE CellML_CreateFieldToCellMLMapC(cellML,field,variableType,componentNumber,fieldParamSetType, &
    & modelIndex,variableID,cellMLParameterSetType,err,error,*)
    
    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<The CellML environment object in which to create the map.
    TYPE(FieldType), POINTER, INTENT(IN) :: field !<The field to map from.
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to map from.
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The field variable component number to map from.
    INTEGER(INTG), INTENT(IN) :: fieldParamSetType !<The field variable parameter set to map from.
    INTEGER(INTG), INTENT(IN) :: modelIndex !<The index of the CellML model to map to.
    CHARACTER(LEN=*), INTENT(IN) :: variableID !<The ID of the CellML variable in the given model to map to.
    INTEGER(INTG), INTENT(IN) :: cellMLParameterSetType !<The CellML parameter set to map to.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local variables
    INTEGER(INTG) :: cellMLFieldType,cellMLVariableNumber,cellMLVariableType,fieldInterpolationType,mapIdx
    !INTEGER(INTG) :: cNameL
    INTEGER(C_INT) :: errorC
    CHARACTER(LEN=1,KIND=C_CHAR) :: cName(MAXSTRLEN)
    TYPE(CellMLFieldMapsType), POINTER :: cellMLFieldMaps
    TYPE(CellMLModelType), POINTER :: cellMLModel
    TYPE(CellMLModelMapType), POINTER :: newCellMLModelMap
    TYPE(CellMLModelMapPtrType), ALLOCATABLE :: newFieldsMappedTo(:)
    TYPE(CellMLModelMapsType), POINTER :: cellMLModelMaps
    TYPE(DomainType), POINTER :: fieldDomain,sourceFieldDomain
    TYPE(FieldType), POINTER :: fieldGeometricField,sourceGeometricField
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(FieldParameterSetType), POINTER :: fieldParameterSet
    TYPE(VARYING_STRING) :: localError
   
    ENTERS("CellML_CreateFieldToCellMLMapC",err,error,*999)

#ifdef WITH_CELLML

    CALL CellML_AssertIsFinished(cellML,err,error,*999)
    NULLIFY(fieldGeometricField)
    CALL Field_GeometricFieldGet(field,fieldGeometricField,err,error,*999)
    NULLIFY(fieldVariable)
    CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
    NULLIFY(fieldParameterSet)
    CALL FieldVariable_ParameterSetGet(fieldVariable,fieldParamSetType,fieldParameterSet,err,error,*999)
    NULLIFY(fieldDomain)
    CALL FieldVariable_ComponentDomainGet(fieldVariable,componentNumber,fieldDomain,err,error,*999)
    CALL FieldVariable_ComponentInterpolationGet(fieldVariable,componentNumber,fieldInterpolationType,err,error,*999)
    NULLIFY(cellMLModel)
    CALL CellML_CellMLModelGet(cellML,modelIndex,cellMLModel,err,error,*999)
    NULLIFY(cellMLFieldMaps)
    CALL CellML_CellMLFieldMapsGet(cellML,cellMLFieldMaps,err,error,*999)
    CALL CellMLFieldMaps_AssertNotFinished(cellMLfieldMaps,err,error,*999)
    NULLIFY(cellMLModelMaps)
    CALL CellMLFieldMaps_CellMLModelMapsGet(cellMLFieldMaps,modelIndex,cellMLModelMaps,err,error,*999)
    CALL CMISSF2CString(variableID,cName)    
    !All input arguments are ok. Get the type of the variable being mapped
    errorC = CELLML_MODEL_DEFINITION_GET_VARIABLE_TYPE(cellMLModel%ptr,cName,cellMLVariableType)
    IF(errorC /= 0) THEN
      localError="Error "//TRIM(NumberToVString(errorC,"*",err,error))// &
        & ". Failed to get the type of CellML variable: "//variableID
      CALL FlagError(localError,err,error,*999)
    ENDIF
    cellMLFieldType=CellML_MapCellMLVariableToFieldType(cellMLVariableType,err,error)
    IF(err/=0) GOTO 999
    CALL CellML_FieldComponentGet(cellML,modelIndex,cellMLFieldType,variableID,cellMLVariableNumber,err,error,*999)
    !cNameL = LEN_TRIM(variableID)
    !WRITE(cName,'(A,A)') cName(1:cNameL),C_NULL_CHAR
    !cellMLVariableNumber=CELLML_MODEL_DEFINITION_ADD_MAPPING_TO_FIELD(cellMLModel%ptr,cName)
    !Now check that the mapped field is consistent with the other mapped fields for the model.
    sourceGeometricField=>cellMLFieldMaps%sourceGeometricField
    IF(ASSOCIATED(sourceGeometricField)) THEN
      IF(.NOT.ASSOCIATED(sourceGeometricField,fieldGeometricField)) THEN
        localError="The geometric field for field user number "//TRIM(NumberToVString(field%userNumber,"*",err,error))// &
          & " does not match the geometric field for other field variable components mapped in the" // &
          & " CellML environment."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      NULLIFY(sourceFieldDomain)
      CALL CellMLFieldMaps_SourceFieldDomainGet(cellMLFieldMaps,sourceFieldDomain,err,error,*999)
      IF(.NOT.ASSOCIATED(sourceFieldDomain,fieldDomain)) THEN
        localError="The domain for component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
          & " of variable type "//TRIM(NumberToVString(variableType,"*",err,error))// &
          & " of field user number "//TRIM(NumberToVString(field%userNumber,"*",err,error))// &
          & " does not match the domain for other field variable components mapped in the CellML environment."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      IF(cellMLFieldMaps%sourceFieldInterpolationType/=fieldInterpolationType) THEN
        localError="The interpolation type of "//TRIM(NumberToVString(fieldInterpolationType,"*",err,error))// &
          & " for component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
          & " of variable type "//TRIM(NumberToVString(variableType,"*",err,error))// &
          & " of field user number "//TRIM(NumberToVString(field%userNumber,"*",err,error))// &
          & " does not match the interpolation type of "// &
          & TRIM(NumberToVString(cellMLFieldMaps%sourceFieldInterpolationType,"*",err,error))// &
          & " used in other field variable components mapped in the CellML environment."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      cellMLFieldMaps%sourceGeometricField=>fieldGeometricField
      cellMLFieldMaps%sourceFieldVariable=>fieldVariable
      cellMLFieldMaps%sourceFieldDomain=>fieldDomain
      cellMLFieldMaps%sourceFieldInterpolationType=fieldInterpolationType
    ENDIF
    !Everything is OK so create the model map field.
    NULLIFY(newCellMLModelMap)
    CALL CellML_ModelMapInitialise(newCellMLModelMap,err,error,*999)
    newCellMLModelMap%cellMLMapType=CELLML_MAP_TO_FIELD_TYPE
    newCellMLModelMap%field=>field
    newCellMLModelMap%variableType=variableType
    newCellMLModelMap%componentNumber=componentNumber
    newCellMLModelMap%fieldParameterSet=fieldParamSetType
    newCellMLModelMap%variableId=variableID
    newCellMLModelMap%cellMLFieldType=cellMLFieldType            
    newCellMLModelMap%cellMLVariableNumber=cellMLVariableNumber
    newCellMLModelMap%cellMLParameterSet=cellMLParameterSetType
    !Put this model map field into the list of to field maps
    ALLOCATE(newFieldsMappedTo(cellMLModelMaps%numberOfFieldsMappedTo+1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new fields mapped to.",err,error,*999)
    DO mapIdx=1,cellMLModelMaps%numberOfFieldsMappedTo
      newFieldsMappedTo(mapIdx)%ptr=>cellMLModelMaps%fieldsMappedTo(mapIdx)%ptr
    ENDDO !mapIdx
    newFieldsMappedTo(cellMLModelMaps%numberOfFieldsMappedTo+1)%ptr=>newCellMLModelMap
    CALL MOVE_ALLOC(newFieldsMappedTo,cellMLModelMaps%fieldsMappedTo)
    cellMLModelMaps%numberOfFieldsMappedTo=cellMLModelMaps%numberOfFieldsMappedTo+1
    
    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"CellML model field -> CellML map:",err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE," Field :",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"   User number             = ",field%userNumber,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"   Variable type           = ",variableType,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"   Component number        = ",componentNumber,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"   Parameter set           = ",fieldParamSetType,err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE," CellML model :",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"   CellML User number      = ",cellML%userNumber,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"   Model index             = ",modelIndex,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"   Variable ID             = ",variableID,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"   CellML field type       = ",cellMLFieldType,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"   CellML variable number  = ",cellMLVariableNumber,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"   CellML parameter set    = ",cellMLParameterSetType,err,error,*999)
    ENDIF

#else

    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*999)

#endif

    EXITS("CellML_CreateFieldToCellMLMapC")
    RETURN
999 ERRORSEXITS("CellML_CreateFieldToCellMLMapC",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_CreateFieldToCellMLMapC

  !
  !=================================================================================================================================
  !

  !>Create a field variable component to CellML model variable map.
  SUBROUTINE CellML_CreateFieldToCellMLMapVS(cellML,field,variableType,componentNumber,fieldParamSetType, &
    & modelUserNumber,variableID,cellMLParameterSetType,err,error,*)
    
    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<The CellML environment object in which to create the map.
    TYPE(FieldType), POINTER, INTENT(IN) :: field !<The field to map from.
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to map from.
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The field variable component number to map from.
    INTEGER(INTG), INTENT(IN) :: fieldParamSetType !<The field variable parameter set to map from.
    INTEGER(INTG), INTENT(IN) :: modelUserNumber !<The user number of the CellML model to map to.
    TYPE(VARYING_STRING), INTENT(IN) :: variableID !<The ID of the CellML variable in the given model to map to.
    INTEGER(INTG), INTENT(IN) :: cellMLParameterSetType !<The CellML parameter set to map to.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local variables

    ENTERS("CellML_CreateFieldToCellMLMapVS",err,error,*999)

#ifdef WITH_CELLML

    CALL CellML_CreateFieldToCellMLMap(cellML,field,variableType,componentNumber,fieldParamSetType, &
      & modelUserNumber,CHAR(variableID),cellMLParameterSetType,err,error,*999)
    
#else

    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*999)

#endif

    EXITS("CellML_CreateFieldToCellMLMapVS")
    RETURN
999 ERRORSEXITS("CellML_CreateFieldToCellMLMapVS",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_CreateFieldToCellMLMapVS

  !
  !=================================================================================================================================
  !

  !>Set the dof in a field specified by a model DOF and component to a value.
  SUBROUTINE CellML_FieldModelDofSet(modelVariable,modelDofIdx,field,variableType,parameterSetIdx,componentIdx, &
     value,err,error,*)
    
    !Argument variables
    TYPE(FieldVariableType), POINTER :: modelVariable !<A pointer to the model field variable
    INTEGER(INTG), INTENT(IN) :: modelDofIdx !<The model dof to set.
    TYPE(FieldType), POINTER :: field !<A pointer to the field to set the value for
    INTEGER(INTG), INTENT(IN) :: variableType !<The variable type to set the value for
    INTEGER(INTG), INTENT(IN) :: parameterSetIdx !<The parameter set index of the field variable to set.
    INTEGER(INTG), INTENT(IN) :: componentIdx !<The component index of the field variable to set.
    REAL(DP), INTENT(IN) :: value !<The value to set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local variables
    INTEGER(INTG) :: derivativeNumber,dofParamIdx,dofType,elementNumber,gaussNumber,gridNumber,nodeNumber,versionNumber
    TYPE(VARYING_STRING) :: localError

    ENTERS("CellML_FieldModelDofSet",err,error,*999)

#ifdef WITH_CELLML

    IF(.NOT.ASSOCIATED(modelVariable)) CALL FlagError("Model variable is not asssociated.",err,error,*999)
    
    IF(modelDofIdx<1.OR.modelDofIdx>modelVariable%numberOfDofs) THEN
      localError="The model DOF index of "//TRIM(NumberToVString(modelDofIdx,"*",err,error))// &
        & " is invalid. The DOF index needs to be > 0 and <= "// &
        & TRIM(NumberToVString(modelVariable%numberOfDofs,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    dofType=modelVariable%dofToParamMap%DOFType(1,modelDofIdx)
    dofParamIdx=modelVariable%dofToParamMap%DOFType(2,modelDofIdx)
    SELECT CASE(dofType)
    CASE(FIELD_CONSTANT_DOF_TYPE)
      CALL Field_ParameterSetUpdateConstant(field,variableType,parameterSetIdx,componentIdx,VALUE,err,error,*999)
    CASE(FIELD_ELEMENT_DOF_TYPE)
      elementNumber=modelVariable%dofToParamMap%elementDOF2ParamMap(1,dofParamIdx)
      CALL Field_ParameterSetUpdateLocalElement(field,variableType,parameterSetIdx,elementNumber,componentIdx,VALUE, &
        & err,error,*999)
    CASE(FIELD_NODE_DOF_TYPE)
      versionNumber=modelVariable%dofToParamMap%nodeDOF2ParamMap(1,dofParamIdx)
      derivativeNumber=modelVariable%dofToParamMap%nodeDOF2ParamMap(2,dofParamIdx)
      nodeNumber=modelVariable%dofToParamMap%nodeDOF2ParamMap(3,dofParamIdx)
      CALL Field_ParameterSetUpdateLocalNode(field,variableType,parameterSetIdx,versionNumber,derivativeNumber,NodeNumber, &
        & componentIdx,VALUE,err,error,*999)
    CASE(FIELD_GRID_POINT_DOF_TYPE)
      gridNumber=modelVariable%dofToParamMap%gridPointDOF2ParamMap(1,dofParamIdx)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(FIELD_GAUSS_POINT_DOF_TYPE)
      gaussNumber=modelVariable%dofToParamMap%gaussPointDOF2ParamMap(1,dofParamIdx)
      elementNumber=modelVariable%dofToParamMap%gaussPointDOF2ParamMap(2,dofParamIdx)
      CALL Field_ParameterSetUpdateLocalGaussPoint(field,variableType,parameterSetIdx,gaussNumber,elementNumber,componentIdx, &
        & VALUE,err,error,*999)
    CASE DEFAULT
      localError="The DOF type of "//TRIM(NumberToVString(dofType,"*",err,error))// &
        & " for DOF number "//TRIM(NumberToVString(modelDofIdx,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
#else

    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*999)

#endif

    EXITS("CellML_FieldModelDofSet")
    RETURN
999 ERRORSEXITS("CellML_FieldModelDofSet",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_FieldModelDofSet

  !
  !=================================================================================================================================
  !

  !>Checks a CellML environment models field for correctness.
  SUBROUTINE CellMLModelsField_Check(cellMLModelsField,err,error,*)
    !Argument variables
    TYPE(CellMLModelsFieldType), POINTER :: cellMLModelsField !<A pointer to the CellML environment models field to check.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !< The error string
    !Local variables
    INTEGER(INTG) :: modelIdx,sourceDOFIdx,firstDOFIdx,mpiIError,onlyOneModelIndex,groupCommunicator
    INTEGER(INTG), POINTER :: modelsData(:)
    TYPE(CellMLType), POINTER :: cellML
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(FieldType), POINTER :: modelsField
    TYPE(FieldVariableType), POINTER :: modelsVariable
    TYPE(VARYING_STRING) :: localError
    TYPE(WorkGroupType), POINTER :: workGroup
    
    ENTERS("CellMLModelsField_Check",err,error,*999)

#ifdef WITH_CELLML

    CALL CellMLModelsField_AssertIsFinished(cellMLModelsField,err,error,*999)
    IF(cellMLModelsField%onlyOneModelIndex==CELLML_MODELS_FIELD_NOT_CHECKED) THEN
      !Models field has not been checked before.
      NULLIFY(cellML)
      CALL CellMLModelsField_CellMLGet(cellMLModelsField,cellML,err,error,*999)
      NULLIFY(modelsField)
      CALL CellMLModelsField_ModelsFieldGet(cellMLModelsField,modelsField,err,error,*999)
      NULLIFY(decomposition)
      CALL Field_DecompositionGet(modelsField,decomposition,err,error,*999)
      NULLIFY(workGroup)
      CALL Decomposition_WorkGroupGet(decomposition,workGroup,err,error,*999)
      CALL WorkGroup_GroupCommunicatorGet(workGroup,groupCommunicator,err,error,*999)
      NULLIFY(modelsVariable)
      CALL Field_VariableGet(modelsField,FIELD_U_VARIABLE_TYPE,modelsVariable,err,error,*999)
      IF(modelsVariable%numberOfDofs<=0) CALL FlagError("CellML models field variable does not have any DOFs.",err,error,*999)
      NULLIFY(modelsData)
      CALL FieldVariable_ParameterSetDataGet(modelsVariable,FIELD_VALUES_SET_TYPE,modelsData,err,error,*999)
      !check for the first non-zero model index
      sourceDOFIdx=0
      firstDOFIdx=1
      DO sourceDOFIdx=1,modelsVariable%totalNumberOfDofs
        modelIdx=modelsData(sourceDOFIdx)
        IF(modelIdx>=0) THEN
          cellMLModelsField%onlyOneModelIndex=modelIdx
          firstDOFIdx=sourceDOFIdx
          IF(modelIdx>0) EXIT
        ELSE
          localError="The model index of "//TRIM(NumberToVString(modelIdx,"*",err,error))// &
            & " is invalid for source DOF 1. The model index must be >= 0 and <= "// &
            & TRIM(NumberToVString(cellML%numberOfModels,"*",err,error))//"."
          CALL FlagError("The models field has not been set for DOF 1.",err,error,*999)
        ENDIF
      ENDDO !sourceDOFIdx
      IF(modelIdx<0.OR.modelIdx>cellML%numberOfModels) THEN
        localError="The model index of "//TRIM(NumberToVString(modelIdx,"*",err,error))// &
          & " is invalid for source DOF 1. The model index must be >= 0 and <= "// &
          & TRIM(NumberToVString(cellML%numberOfModels,"*",err,error))//"."
        CALL FlagError("The models field has not been set for DOF 1.",err,error,*999)
      ENDIF
      DO sourceDOFIdx=(firstDOFIdx+1),modelsVariable%totalNumberOfDofs
        modelIdx=modelsData(sourceDOFIdx)
        IF(modelIdx<0.OR.modelIdx>cellML%numberOfModels) THEN
          localError="The model index of "//TRIM(NumberToVString(modelIdx,"*",err,error))// &
            & " is invalid for source DOF "//TRIM(NumberToVString(sourceDOFIdx,"*",err,error))// &
            & ". The model index must be >= 0 and <= "// &
            & TRIM(NumberToVString(cellML%numberOfModels,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        IF(modelIdx/=cellMLModelsField%onlyOneModelIndex.AND.modelIdx/=0) THEN
          cellMLModelsField%onlyOneModelIndex=CELLML_MODELS_FIELD_NOT_CONSTANT
          EXIT
        ENDIF
      ENDDO !sourceDOFIdx
      onlyOneModelIndex=0
      CALL MPI_ALLREDUCE(cellMLModelsField%onlyOneModelIndex,onlyOneModelIndex,1,MPI_INTEGER,MPI_MAX,groupCommunicator,mpiIerror)
      CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",mpiIerror,err,error,*999)
      IF(onlyOneModelIndex==0) CALL FlagError("Models field does not have any models set.",err,error,*999)
!!TODO: Do we need to make sure it is the same model number on different ranks? The only one model optimisation is to ensure
!!that we don't have to reference the models field inside dof loops on the rank??? 
      CALL FieldVariable_ParameterSetDataRestore(modelsVariable,FIELD_VALUES_SET_TYPE,modelsData,err,error,*999)
    ENDIF

#else

    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*999)

#endif

    EXITS("CellMLModelsField_Check")
    RETURN
999 ERRORSEXITS("CellMLModelsField_Check",err,error)
    RETURN 1
    
  END SUBROUTINE CellMLModelsField_Check

  !
  !=================================================================================================================================
  !

  !>Start the creation of the models field for the given CellML environment.
  SUBROUTINE CellML_ModelsFieldCreateStart(modelsFieldUserNumber,cellML,modelsField,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: modelsFieldUserNumber !<The unique identifier for the models field to be created for the given CellML environment object.
    TYPE(CellMLType), POINTER :: cellML !<The CellML environment object for which we need to create the models field.
    TYPE(FieldType), POINTER :: modelsField !<If associated on entry, a pointer to the user created models field which has the same user number as the specified models field user number. If not associated on entry, on exit, a pointer to the created models field for the CellML environment.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local variables
    INTEGER(INTG) :: sourceMeshComponentNumber
    TYPE(CellMLFieldMapsType), POINTER :: cellMLFieldMaps
    TYPE(CellMLModelsFieldType), POINTER :: cellMLModelsField
    TYPE(DecompositionType), POINTER :: modelsFieldDecomposition,sourceFieldDecomposition
    TYPE(DomainType), POINTER :: sourceFieldDomain
    TYPE(FieldType), POINTER :: field,modelsGeometricField,sourceGeometricField
    TYPE(RegionType), POINTER :: region,modelsFieldRegion
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("CellML_ModelsFieldCreateStart",err,error,*999)

#ifdef WITH_CELLML

    CALL CellML_AssertIsFinished(cellML,err,error,*999)
    NULLIFY(cellMLFieldMaps)
    CALL CellML_CellMLFieldMapsGet(cellML,cellMLFieldMaps,err,error,*999)
    CALL CellMLFieldMaps_AssertIsFinished(cellMLFieldMaps,err,error,*999)
    IF(ASSOCIATED(cellML%modelsField)) CALL FlagError("The CellML environment models field is already associated.",err,error,*999)

    NULLIFY(region)
    CALL CellML_RegionGet(cellML,region,err,error,*999)
    NULLIFY(sourceGeometricField)
    CALL CellMLFieldMaps_SourceGeometricFieldGet(cellMLFieldMaps,sourceGeometricField,err,error,*999)
    NULLIFY(sourceFieldDomain)
    CALL CellMLFieldMaps_SourceFieldDomainGet(cellMLFieldMaps,sourceFieldDomain,err,error,*999)
    CALL Domain_MeshComponentNumberGet(sourceFieldDomain,sourceMeshComponentNumber,err,error,*999)
    
    IF(ASSOCIATED(modelsField)) THEN
      !Check the field has been finished
      CALL Field_AssertIsFinished(modelsField,err,error,*999)
      !Check the user numbers match
      IF(modelsFieldUserNumber/=modelsField%userNumber) THEN
        localError="The specified models field user number of "//TRIM(NumberToVString(modelsFieldUserNumber,"*",err,error))// &
          & " does not match the user number of the specified models field of "// &
          & TRIM(NumberToVString(modelsField%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      NULLIFY(modelsFieldRegion)
      CALL Field_RegionGet(modelsField,modelsFieldRegion,err,error,*999)
      !Check the field is defined on the same region as the CellML region
      IF(modelsFieldRegion%userNumber/=region%userNumber) THEN
        localError="Invalid region setup. The specified models field has been created on region number "// &
          & TRIM(NumberToVString(modelsFieldRegion%userNumber,"*",err,error))// &
          & " and the specified CellML environment has been created on region number "// &
          & TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      NULLIFY(modelsGeometricField)
      CALL Field_GeometricFieldGet(modelsField,modelsGeometricField,err,error,*999)
      !Check the specified models field has the same geometric field as the source field
      IF(.NOT.ASSOCIATED(sourceGeometricField,modelsGeometricField)) THEN
        CALL FlagError("The specified models field does not have the same geometric field as the "// &
          & "geometric field for the specified CellML environment.",err,error,*999)
      ENDIF
      !Check the specified models field has the same decomposition as the source field
      NULLIFY(sourceFieldDecomposition)
      CALL Domain_DecompositionGet(sourceFieldDomain,sourceFieldDecomposition,err,error,*999)
      NULLIFY(modelsFieldDecomposition)
      CALL Field_DecompositionGet(modelsField,modelsFieldDecomposition,err,error,*999)
      IF(.NOT.ASSOCIATED(sourceFieldDecomposition,modelsFieldDecomposition)) THEN
        CALL FlagError("The specified models field does not have the same decomposition as the source "// &
          & "domain decomposition for the specified CellML environment.",err,error,*999)
      ENDIF
    ELSE
      !Check the user number has not already been used for a field in this region.
      NULLIFY(field)
      CALL Field_UserNumberFind(modelsFieldUserNumber,region,field,err,error,*999)
      IF(ASSOCIATED(field)) THEN
        localError="The specified models field user number of "//TRIM(NumberToVString(modelsFieldUserNumber,"*",err,error))// &
          & "has already been used to create a field on region number "//TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDIF
    CALL CellML_ModelsFieldInitialise(cellML,err,error,*999)
    NULLIFY(cellMLModelsField)
    CALL CellML_CellMLModelsFieldGet(cellML,cellMLModelsField,err,error,*999)
    IF(ASSOCIATED(modelsField)) THEN
      !Now check the supplied field.
      CALL Field_DataTypeCheck(modelsField,FIELD_U_VARIABLE_TYPE,FIELD_INTG_TYPE,err,error,*999)
      CALL Field_TypeCheck(modelsField,FIELD_GENERAL_TYPE,err,error,*999)
      CALL Field_NumberOfVariablesCheck(modelsField,1,err,error,*999)
      CALL Field_VariableTypesCheck(modelsField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
      CALL Field_NumberOfComponentsCheck(modelsField,FIELD_U_VARIABLE_TYPE,1,err,error,*999)
      CALL Field_ComponentMeshComponentCheck(modelsField,FIELD_U_VARIABLE_TYPE,1,sourceMeshComponentNumber,err,error,*999)
      CALL Field_ComponentInterpolationCheck(modelsField,FIELD_U_VARIABLE_TYPE,1,cellMLFieldMaps%sourceFieldInterpolationType, &
        & err,error,*999)
    ELSE
      !Create the CellML environment models field
      cellMLModelsField%modelsFieldAutoCreated=.TRUE.     
      CALL Field_CreateStart(modelsFieldUserNumber,region,cellMLModelsField%modelsField,err,error,*999)
      CALL Field_DataTypeSetAndLock(cellMLModelsField%modelsField,FIELD_U_VARIABLE_TYPE,FIELD_INTG_TYPE,err,error,*999)
      CALL Field_LabelSet(cellMLModelsField%modelsField,"CellMLModelsField",err,error,*999)
      CALL Field_TypeSetAndLock(cellMLModelsField%modelsField,FIELD_GENERAL_TYPE,err,error,*999)
      CALL Field_DecompositionSetAndLock(cellMLModelsField%modelsField,sourceFieldDecomposition,err,error,*999)
      CALL Field_GeometricFieldSetAndLock(cellMLModelsField%modelsField,sourceGeometricField,err,error,*999)
      CALL Field_NumberOfVariablesSetAndLock(cellMLModelsField%modelsField,1,err,error,*999)
      CALL Field_VariableTypesSetAndLock(cellMLModelsField%modelsField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
      CALL Field_VariableLabelSet(cellMLModelsField%modelsField,FIELD_U_VARIABLE_TYPE,"ModelMap",err,error,*999)
      CALL Field_DOFOrderTypeSet(cellMLModelsField%modelsField,FIELD_U_VARIABLE_TYPE,FIELD_CONTIGUOUS_COMPONENT_DOF_ORDER, &
        & err,error,*999)
      CALL Field_NumberOfComponentsSetAndLock(cellMLModelsField%modelsField,FIELD_U_VARIABLE_TYPE,1, err,error,*999)
      CALL Field_ComponentLabelSet(cellMLModelsField%modelsField,FIELD_U_VARIABLE_TYPE,1,"ModelUserNumber",err,error,*999)
      CALL Field_ComponentMeshComponentSetAndLock(cellMLModelsField%modelsField,FIELD_U_VARIABLE_TYPE,1,sourceMeshComponentNumber, &
        & err,error,*999)
      CALL Field_ComponentInterpolationSetAndLock(cellMLModelsField%modelsField,FIELD_U_VARIABLE_TYPE,1, &
        & cellMLFieldMaps%sourceFieldInterpolationType,err,error,*999)
    ENDIF
    !Set pointers
    IF(cellMLModelsField%modelsFieldAutoCreated) THEN            
      modelsField=>cellMLModelsField%modelsField
    ELSE
      cellMLModelsField%modelsField=>modelsField
    ENDIF
    
#else
    
    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*999)
    
#endif
    
    EXITS("CellML_ModelsFieldCreateStart")
    RETURN
999 ERRORSEXITS("CellML_ModelsFieldCreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_ModelsFieldCreateStart

  !
  !=================================================================================================================================
  !

  !>Finish the creation of the models field for the given CellML environment.
  SUBROUTINE CellML_ModelsFieldCreateFinish(cellML,err,error,*)
    
    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<The CellML environment object for which we need to finish creation of the models field.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local variables
    TYPE(CellMLModelsFieldType), POINTER :: cellMLModelsField

    ENTERS("CellML_ModelsFieldCreateFinish",err,error,*999)

#ifdef WITH_CELLML

    NULLIFY(cellMLModelsField)
    CALL CellML_CellMLModelsFieldGet(cellML,cellMLModelsField,err,error,*999)
    CALL CellMLModelsField_AssertNotFinished(cellMLModelsField,err,error,*999)
    
    !Finish the models field creation
    IF(cellMLModelsField%modelsFieldAutoCreated) CALL Field_CreateFinish(cellMLModelsField%modelsField,err,error,*999)
    cellMLModelsField%modelsFieldFinished=.TRUE.
    !Default the models field to the first model
    CALL Field_ComponentValuesInitialise(cellMLModelsField%modelsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1, &
      & err,error,*999)

#else

    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*999)

#endif

    EXITS("CellML_ModelsFieldCreateFinish")
    RETURN
999 ERRORSEXITS("CellML_ModelsFieldCreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_ModelsFieldCreateFinish

  !
  !=================================================================================================================================
  !

  !>Finalise a CellML environment models field and deallocate all memory.
  SUBROUTINE CellML_ModelsFieldFinalise(cellMLModelsField,err,error,*)
    
    !Argument variables
    TYPE(CellMLModelsFieldType), POINTER :: cellMLModelsField !<A pointer to the CellML environment models field to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !< The error string
    !Local variables
    
    ENTERS("CellML_ModelsFieldFinalise",err,error,*999)

#ifdef WITH_CELLML

    IF(ASSOCIATED(cellMLModelsField)) THEN
      DEALLOCATE(cellMLModelsField)
    ENDIF

#else

    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*999)

#endif

    EXITS("CellML_ModelsFieldFinalise")
    RETURN
999 ERRORSEXITS("CellML_ModelsFieldFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_ModelsFieldFinalise

  !
  !=================================================================================================================================
  !

  !>Initialise a CellML environment models field.
  SUBROUTINE CellML_ModelsFieldInitialise(cellML,err,error,*)
    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML environment to initialise the models field for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !< The error string
    !Local variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("CellML_ModelsFieldInitialise",err,error,*998)

#ifdef WITH_CELLML

    IF(.NOT.ASSOCIATED(cellML))  CALL FlagError("CellML environment is not associated.",err,error,*998)    
    IF(ASSOCIATED(cellML%modelsField)) CALL FlagError("CellML environment models field is already associated.",err,error,*998)
     
    ALLOCATE(cellML%modelsField,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate CellML environment models field.",err,error,*999)
    cellML%modelsField%cellML=>cellML
    cellML%modelsField%modelsFieldFinished=.FALSE.
    cellML%modelsField%modelsFieldAutoCreated=.FALSE.
    NULLIFY(cellML%modelsField%modelsField)
    cellML%modelsField%onlyOneModelIndex=CELLML_MODELS_FIELD_NOT_CHECKED

#else

    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*998)

#endif

    EXITS("CellML_ModelsFieldInitialise")
    RETURN
999 CALL CellML_ModelsFieldFinalise(cellML%modelsField,dummyErr,dummyError,*998)
998 ERRORSEXITS("CellML_ModelsFieldInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_ModelsFieldInitialise

  !
  !=================================================================================================================================
  !

  !>Start the creation of the state field for the given CellML environment.
  SUBROUTINE CellML_StateFieldCreateStart(stateFieldUserNumber,cellML,stateField,err,error,*)
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: stateFieldUserNumber !<The unique identifier for the state field to be created for the given CellML environment object.
    TYPE(CellMLType), POINTER :: cellML !<The CellML environment object for which to create the state field.
    TYPE(FieldType), POINTER :: stateField  !<If associated on entry, a pointer to the user created state field which has the same user number as the specified state field user number. If not associated on entry, on exit, a pointer to the created state field for the CellML environment.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local variables
    INTEGER(INTG) :: componentIdx,sourceMeshComponentNumber
    TYPE(CellMLFieldMapsType), POINTER :: cellMLFieldMaps
    TYPE(CellMLStateFieldType), POINTER :: cellMLStateField
    TYPE(DecompositionType), POINTER :: sourceFieldDecomposition,stateFieldDecomposition
    TYPE(DomainType), POINTER :: sourceFieldDomain
    TYPE(FieldType), POINTER :: field,sourceGeometricField,stateGeometricFIeld
    TYPE(RegionType), POINTER :: region,stateFieldRegion
    TYPE(VARYING_STRING) :: localError

    ENTERS("CellML_StateFieldCreateStart",err,error,*999)

#ifdef WITH_CELLML

    CALL CellML_AssertIsFinished(cellML,err,error,*999)
    NULLIFY(cellMLFieldMaps)
    CALL CellML_CellMLFieldMapsGet(cellML,cellMLFieldMaps,err,error,*999)
    CALL CellMLFieldMaps_AssertIsFinished(cellMLFieldMaps,err,error,*999)
    IF(ASSOCIATED(cellML%stateField)) CALL FlagError("The CellML environment state field is already associated.",err,error,*999)

    NULLIFY(region)
    CALL CellML_RegionGet(cellML,region,err,error,*999)
    NULLIFY(sourceGeometricField)
    CALL CellMLFieldMaps_SourceGeometricFieldGet(cellMLFieldMaps,sourceGeometricField,err,error,*999)
    NULLIFY(sourceFieldDomain)
    CALL CellMLFieldMaps_SourceFieldDomainGet(cellMLFieldMaps,sourceFieldDomain,err,error,*999)
    NULLIFY(sourceFieldDecomposition)
    CALL Domain_DecompositionGet(sourceFieldDomain,sourceFieldDecomposition,err,error,*999)
    CALL Domain_MeshComponentNumberGet(sourceFieldDomain,sourceMeshComponentNumber,err,error,*999)
   
    IF(ASSOCIATED(stateField)) THEN
      !Check the field has been finished
      CALL Field_AssertIsFinished(stateField,err,error,*999)
      !Check the user numbers match
      IF(stateFieldUserNumber/=stateField%userNumber) THEN
        localError="The specified state field user number of "// &
          & TRIM(NumberToVString(stateFieldUserNumber,"*",err,error))// &
          & " does not match the user number of the specified state field of "// &
          & TRIM(NumberToVString(stateField%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      !Check that the state field has the same region as the CellML environment
      NULLIFY(stateFieldRegion)
      CALL Field_RegionGet(stateField,stateFieldRegion,err,error,*999)
      IF(stateFieldRegion%userNumber/=region%userNumber) THEN
        localError="Invalid region setup. The specified state field has been created on region number "// &
          & TRIM(NumberToVString(stateFieldRegion%userNumber,"*",err,error))// &
          & " and the CellML environment has been created on region number "// &
          & TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      !Check the specified models field has the same geometric field as the source field
      NULLIFY(stateGeometricField)
      CALL Field_GeometricFieldGet(stateField,stateGeometricField,err,error,*999)
      IF(.NOT.ASSOCIATED(sourceGeometricField,stateGeometricField)) THEN
        CALL FlagError("The specified state field does not have the same geometric field as the "// &
          & "geometric field for the specified CellML environment.",err,error,*999)
      ENDIF
      !Check the specified models field has the same decomposition as the source field
      NULLIFY(stateFieldDecomposition)
      CALL Field_DecompositionGet(stateField,stateFieldDecomposition,err,error,*999)
      IF(.NOT.ASSOCIATED(sourceFieldDecomposition,stateFieldDecomposition)) THEN
        CALL FlagError("The specified state field does not have the same decomposition as the source "// &
          & "domain decomposition for the specified CellML environment.",err,error,*999)
      ENDIF
    ELSE
      !Check the user number has not already been used for a field in this region.
      NULLIFY(field)
      CALL Field_UserNumberFind(stateFieldUserNumber,region,field,err,error,*999)
      IF(ASSOCIATED(field)) THEN
        localError="The specified state field user number of "//TRIM(NumberToVString(stateFieldUserNumber,"*",err,error))// &
          & "has already been used to create a field on region number "//TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDIF
    CALL CellML_StateFieldInitialise(cellML,err,error,*999)
    NULLIFY(cellMLStateField)
    CALL CellML_CellMLStateFieldGet(cellML,cellMLStateField,err,error,*999)
    IF(ASSOCIATED(stateField)) THEN
      !Now check the supplied field.
      CALL Field_DataTypeCheck(stateField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
      CALL Field_TypeCheck(stateField,FIELD_GENERAL_TYPE,err,error,*999)
      CALL Field_NumberOfVariablesCheck(stateField,1,err,error,*999)
      CALL Field_VariableTypesCheck(stateField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
      CALL Field_NumberOfComponentsCheck(stateField,FIELD_U_VARIABLE_TYPE,cellML%maximumNumberOfState,err,error,*999)
      DO componentIdx=1,cellML%maximumNumberOfState
        CALL Field_ComponentMeshComponentCheck(stateField,FIELD_U_VARIABLE_TYPE,componentIdx,sourceMeshComponentNumber, &
          & err,error,*999)
        CALL Field_ComponentInterpolationCheck(stateField,FIELD_U_VARIABLE_TYPE,componentIdx, &
          & cellMLFieldMaps%sourceFieldInterpolationType,err,error,*999)
      ENDDO !componentIdx
    ELSE
      cellMLStateField%stateFieldAutoCreated=.TRUE.
      !Create the CellML environment models field
      CALL Field_CreateStart(stateFieldUserNumber,region,cellMLStateField%stateField,err,error,*999)
      CALL Field_DataTypeSetAndLock(cellMLStateField%stateField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
      CALL Field_LabelSet(cellMLStateField%stateField,"CellMLStateField",err,error,*999)
      CALL Field_TypeSetAndLock(cellMLStateField%stateField,FIELD_GENERAL_TYPE,err,error,*999)
      CALL Field_DecompositionSetAndLock(cellMLStateField%stateField,sourceFieldDecomposition,err,error,*999)
      CALL Field_GeometricFieldSetAndLock(cellMLStateField%stateField,sourceGeometricField,err,error,*999)
      CALL Field_NumberOfVariablesSetAndLock(cellMLStateField%stateField,1,err,error,*999)
      CALL Field_VariableTypesSetAndLock(cellMLStateField%stateField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
      CALL Field_VariableLabelSet(cellMLStateField%stateField,FIELD_U_VARIABLE_TYPE,"StateVariable",err,error,*999)
      CALL Field_DOFOrderTypeSet(cellMLStateField%stateField,FIELD_U_VARIABLE_TYPE,FIELD_CONTIGUOUS_COMPONENT_DOF_ORDER, &
        & err,error,*999)
      CALL Field_NumberOfComponentsSetAndLock(cellML%stateField%stateField,FIELD_U_VARIABLE_TYPE,cellML%maximumNumberOfState, &
        & err,error,*999)
      DO componentIdx=1,cellML%maximumNumberOfState
        CALL Field_ComponentMeshComponentSetAndLock(cellMLStateField%stateField,FIELD_U_VARIABLE_TYPE,componentIdx, &
          & sourceMeshComponentNumber,err,error,*999)
        CALL Field_ComponentInterpolationSetAndLock(cellML%stateField%stateField,FIELD_U_VARIABLE_TYPE,componentIdx, &
          & cellMLFieldMaps%sourceFieldInterpolationType,err,error,*999)
      ENDDO !componentIdx
    ENDIF
    !Set pointers
    IF(cellMLStateField%stateFieldAutoCreated) THEN            
      stateField=>cellMLStateField%stateField
    ELSE
      cellMLStateField%stateField=>stateField
    ENDIF
 
#else

    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*999)

#endif

    EXITS("CellML_StateFieldCreateStart")
    RETURN
999 ERRORSEXITS("CellML_StateFieldCreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_StateFieldCreateStart

  !
  !=================================================================================================================================
  !

  !>Finialse the creation of the state field for the given CellML environment.
  !! Finish creating the state variable field for the provided CellML environment.
  !! - default field values are set for all components in the field at all DOFs based on the initial_value's from the CellML model
  SUBROUTINE CellML_StateFieldCreateFinish(cellML,err,error,*)
    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<The CellML environment object for which to finalise the creation of the state field.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local variables
    INTEGER(INTG) :: cellMLVariableType,errorCode,modelIdx,modelsDOFIdx,stateComponentIdx
    INTEGER(INTG), POINTER :: modelsData(:)
    REAL(DP) :: initialValue
    TYPE(CellMLFieldMapsType), POINTER :: cellMLFieldMaps
    TYPE(CellMLModelType), POINTER :: cellMLModel
    TYPE(CellMLModelsFieldType), POINTER :: cellMLModelsField
    TYPE(CellMLStateFieldType), POINTER :: cellMLStateField
    TYPE(FieldVariableType), POINTER :: modelsVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("CellML_StateFieldCreateFinish",err,error,*999)

#ifdef WITH_CELLML

    NULLIFY(cellMLStateField)
    CALL CellML_CellMLStateFieldGet(cellML,cellMLStateField,err,error,*999)
    CALL CellMLStateField_AssertNotFinished(cellMLStateField,err,error,*999)
    NULLIFY(cellMLModelsField)
    CALL CellML_CellMLModelsFieldGet(cellML,cellMLModelsField,err,error,*999)
    CALL CellMLModelsField_AssertIsFinished(cellMLModelsField,err,error,*999)
    CALL CellMLModelsField_Check(cellML%modelsField,err,error,*999)
    !Finish the state field creation
    IF(cellMLStateField%stateFieldAutoCreated) CALL Field_CreateFinish(cellMLStateField%stateField,err,error,*999)              
    !Set the default field values to the initial CellML values.
    IF(cellML%modelsField%onlyOneModelIndex/=0) THEN
      !If there are CellML models on this rank
      IF(cellML%modelsField%onlyOneModelIndex/=CELLML_MODELS_FIELD_NOT_CONSTANT) THEN
        !Only one model so optimise
        NULLIFY(cellMLModel)
        CALL CellML_CellMLModelGet(cellML,cellMLModelsField%onlyOneModelIndex,cellMLModel,err,error,*999)
        DO stateComponentIdx=1,cellMLModel%numberOfState
          cellMLVariableType=CellML_MapCellMLFieldTypeToVariableType(CELLML_STATE_FIELD,err,error)
          IF(err/=0) GOTO 999
          errorCode = CELLML_MODEL_DEFINITION_GET_INITIAL_VALUE_BY_INDEX(cellMLModel%ptr,cellMLVariableType,&
            & stateComponentIdx,initialValue)
          IF(errorCode/=0) THEN
            !problem getting the initial value
            localError="Error "//TRIM(NumberToVString(errorCode,"*",err,error))// &
              & ". Failed to get an initial value for the state variable with component number "//&
              & TRIM(NumberToVString(stateComponentIdx,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          CALL Field_ComponentValuesInitialise(cellMLStateField%stateField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & stateComponentIdx,initialValue,err,error,*999)
        ENDDO !stateComponentIdx
      ELSE
        !Multiple models so go through each dof.
        NULLIFY(cellMLFieldMaps)
        CALL CellML_CellMLFieldMapsGet(cellML,cellMLFieldMaps,err,error,*999)
        NULLIFY(modelsVariable)
        CALL Field_VariableGet(cellMLModelsField%modelsField,FIELD_U_VARIABLE_TYPE,modelsVariable,err,error,*999)
        NULLIFY(modelsData)
        CALL FieldVariable_ParameterSetDataGet(modelsVariable,FIELD_VALUES_SET_TYPE,modelsData,err,error,*999)
        DO modelsDOFIdx=1,modelsVariable%numberOfDofs
          modelIdx=modelsData(modelsDOFIdx)
          IF(modelIdx>0) THEN
            NULLIFY(cellMLModel)
            CALL CellML_CellMLModelGet(cellML,modelIdx,cellMLModel,err,error,*999)
            DO stateComponentIdx=1,cellMLModel%numberOfState
              cellMLVariableType=CellML_MapCellMLFieldTypeToVariableType(CELLML_STATE_FIELD,err,error)              
              IF(err/=0) GOTO 999
              errorCode = CELLML_MODEL_DEFINITION_GET_INITIAL_VALUE_BY_INDEX(cellMLModel%ptr,cellMLVariableType,&
                & stateComponentIdx,initialValue)
              IF(errorCode/=0) THEN
                !problem getting the initial value
                localError="Error "//TRIM(NumberToVString(errorCode,"*",err,error))// &
                  & ". Failed to get an initial value for the state variable with component number "//&
                  & TRIM(NumberToVString(stateComponentIdx,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              ENDIF
              CALL CellML_FieldModelDofSet(modelsVariable,modelsDOFIdx,cellMLStateField%stateField, &
                & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,stateComponentIdx,initialValue,err,error,*999)
            ENDDO !stateComponentIdx
          ENDIF
        ENDDO !modelsDOFIdx
        CALL FieldVariable_ParameterSetDataRestore(modelsVariable,FIELD_VALUES_SET_TYPE,modelsData,err,error,*999)
      ENDIF
    ENDIF
    cellML%stateField%stateFieldFinished=.TRUE.

#else
    
    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*999)
    
#endif

    EXITS("CellML_StateFieldCreateFinish")
    RETURN
999 ERRORSEXITS("CellML_StateFieldCreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_StateFieldCreateFinish

  !
  !=================================================================================================================================
  !

  !>Finalise a CellML environment state field and deallocate all memory.
  SUBROUTINE CellML_StateFieldFinalise(cellMLStateField,err,error,*)
    !Argument variables
    TYPE(CellMLStateFieldType), POINTER :: cellMLStateField !<A pointer to the CellML environment state field to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !< The error string
    !Local variables
    
    ENTERS("CellML_StateFieldFinalise",err,error,*999)

#ifdef WITH_CELLML

    IF(ASSOCIATED(cellMLStateField)) THEN
      DEALLOCATE(cellMLStateField)
    ENDIF

#else

    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*999)

#endif

    EXITS("CellML_StateFieldFinalise")
    RETURN
999 ERRORSEXITS("CellML_StateFieldFinalise",err,error)
    RETURN 1
  END SUBROUTINE CellML_StateFieldFinalise

  !
  !=================================================================================================================================
  !

  !>Initialise a CellML environment models field.
  SUBROUTINE CellML_StateFieldInitialise(cellML,err,error,*)
    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML environment to initialise the models field for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !< The error string
    !Local variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("CellML_StateFieldInitialise",err,error,*998)

#ifdef WITH_CELLML

    IF(ASSOCIATED(cellML)) CALL FlagError("CellML environment is not associated.",err,error,*998)
    IF(ASSOCIATED(cellML%stateField)) CALL FlagError("CellML environment state field is already associated.",err,error,*998)
      
    ALLOCATE(cellML%stateField,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate CellML environment state field.",err,error,*999)
    cellML%stateField%cellML=>cellML
    cellML%stateField%stateFieldFinished=.FALSE.
    cellML%stateField%stateFieldAutoCreated=.FALSE.
    NULLIFY(cellML%stateField%stateField)

#else

    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*998)

#endif

    EXITS("CellML_StateFieldInitialise")
    RETURN
999 CALL CellML_StateFieldFinalise(cellML%stateField,dummyErr,dummyError,*998)
998 ERRORSEXITS("CellML_StateFieldInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_StateFieldInitialise

  !
  !=================================================================================================================================
  !

  !>Find the component ID in the given field for the variable defined by the given variable ID in the provided CellML environment.
  !! This generic routine will be used to map variable ID's in CellML models to components in the various fields defined in the CellML models defined for the provided CellML environment.
  !! - may need to also provide a FIELD_VARIABLE_NUMBER (always 1?) for completeness
  !! - is the model ID also needed?
  !! - because the CellML fields should all be set up to allow direct use in the CellML code, the component number matches the index of the given variable in its associated array in the CellML generated code.
 SUBROUTINE CellML_FieldComponentGetC(cellML,modelIndex,cellMLFieldType,variableID,componentUserNumber,err,error,*)
   !Argument variables
   TYPE(CellMLType), POINTER :: cellML !<The CellML environment object from which to get the field component.
   INTEGER(INTG), INTENT(IN) :: modelIndex !<The index of the CellML model to map from.
   INTEGER(INTG), INTENT(IN) :: cellMLFieldType !<The type of CellML field type to get the component for. \see CELLML_FieldTypes,CmissCellML
   CHARACTER(LEN=*), INTENT(IN) :: variableID !<The ID of the model variable which needs to be located in the provided field.
   INTEGER(INTG), INTENT(OUT) :: componentUserNumber !<On return, the field component for the model variable defined by the given ID.
   INTEGER(INTG), INTENT(OUT) :: err !<The error code.
   TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
   !Local variables
   INTEGER(INTG) :: cellMLVariableIndex
   INTEGER(C_INT) :: errorC
   TYPE(CellMLModelType), POINTER :: cellMLModel
   TYPE(VARYING_STRING) :: localError
   CHARACTER(LEN=1,KIND=C_CHAR) :: cName(MAXSTRLEN)
   !INTEGER(INTG) :: cNameL

   ENTERS("CellML_FieldComponentGetC",err,error,*999)

#ifdef WITH_CELLML

   NULLIFY(cellMLModel)
   CALL CellML_CellMLModelGet(cellML,modelIndex,cellMLModel,err,error,*999)
   CALL CMISSF2CString(variableID,cName)
   errorC = CELLML_MODEL_DEFINITION_GET_VARIABLE_INDEX(cellMLModel%ptr,cName,cellMLVariableIndex)
   IF(errorC /= 0) THEN
     localError="Error "//TRIM(NumberToVString(errorC,"*",err,error))//". Failed to get the index for CellML variable: "// &
       & variableID
     CALL FlagError(localError,err,error,*999)
   ENDIF
   componentUserNumber=cellMLVariableIndex

#else

   CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*999)

#endif

   EXITS("CellML_FieldComponentGetC")
   RETURN
999 ERRORSEXITS("CellML_FieldComponentGetC",err,error)
   RETURN 1
   
 END SUBROUTINE CellML_FieldComponentGetC
  !
  !=================================================================================================================================
  !

  !>Find the component ID in the given field for the variable defined by the given variable ID in the provided CellML environment.
  !! This generic routine will be used to map variable ID's in CellML models to components in the various fields defined in the CellML models defined for the provided CellML environment.
  !! - may need to also provide a FIELD_VARIABLE_NUMBER (always 1?) for completeness
  !! - is the model ID also needed?
  SUBROUTINE CellML_FieldComponentGetVS(cellML,modelIndex,cellMLFieldType,variableID,componentUserNumber,err,error,*)
    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<The CellML environment object from which to get the field component.
    INTEGER(INTG), INTENT(IN) :: modelIndex !<The index of the CellML model to map from.
    INTEGER(INTG), INTENT(IN) :: cellMLFieldType !<The type of CellML field type to get the component for. \see CELLML_FieldTypes,CmissCellML
    TYPE(VARYING_STRING), INTENT(IN) :: variableID !<The ID of the model variable which needs to be located in the provided field.
    INTEGER(INTG), INTENT(OUT) :: componentUserNumber !<On return, the field component for the model variable defined by the given ID.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local variables

    ENTERS("CellML_FieldComponentGetVS",err,error,*999)

#ifdef WITH_CELLML

    CALL CellML_FieldComponentGet(cellML,modelIndex,cellMLFieldType,CHAR(variableID),componentUserNumber,err,error,*999)

#else

    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*999)

#endif

    EXITS("CellML_FieldComponentGetVS")
    RETURN
999 ERRORSEXITS("CellML_FieldComponentGetVS",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_FieldComponentGetVS

  !
  !=================================================================================================================================
  !

  !>Create a field used to store intermediate variables of interest.
  SUBROUTINE CellML_IntermediateFieldCreateStart(intermediateFieldUserNumber,cellML,intermediateField,err,error,*)
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: intermediateFieldUserNumber !<The unique identifier for the intermediate field to be created for the given CellML environment object.
    TYPE(CellMLType), POINTER :: cellML !<The CellML environment object from which to get the field component.
    TYPE(FieldType), POINTER :: intermediateField !<If associated on entry, a pointer to the user created intermediate field which has the same user number as the specified intermediate field user number. If not associated on entry, on exit, a pointer to the created intermediate field for the CellML environment.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local variables
    INTEGER(INTG) :: componentIdx,sourceMeshComponentNumber
    TYPE(CellMLFieldMapsType), POINTER :: cellMLFieldMaps
    TYPE(CellMLIntermediateFieldType), POINTER :: cellMLIntermediateField
    TYPE(DecompositionType), POINTER :: intermediateFieldDecomposition,sourceFieldDecomposition
    TYPE(DomainType), POINTER :: sourceFieldDomain
    TYPE(FieldType), POINTER :: field,intermediateGeometricField,sourceGeometricField
    TYPE(RegionType), POINTER :: region,intermediateFieldRegion
    TYPE(VARYING_STRING) :: localError

    ENTERS("CellML_IntermediateFieldCreateStart",err,error,*999)

#ifdef WITH_CELLML

    CALL CellML_AssertIsFinished(cellML,err,error,*999)
    NULLIFY(cellMLFieldMaps)
    CALL CellML_CellMLFieldMapsGet(cellML,cellMLFieldMaps,err,error,*999)
    CALL CellMLFieldMaps_AssertIsFinished(cellMLFieldMaps,err,error,*999)
    IF(ASSOCIATED(cellML%intermediateField)) &
      & CALL FlagError("The CellML environment intermediate field is already associated.",err,error,*999)
    
    NULLIFY(region)
    CALL CellML_RegionGet(cellML,region,err,error,*999)
    NULLIFY(sourceGeometricField)
    CALL CellMLFieldMaps_SourceGeometricFieldGet(cellMLFieldMaps,sourceGeometricField,err,error,*999)
    NULLIFY(sourceFieldDomain)
    CALL CellMLFieldMaps_SourceFieldDomainGet(cellMLFieldMaps,sourceFieldDomain,err,error,*999)
    NULLIFY(sourceFieldDecomposition)
    CALL Domain_DecompositionGet(sourceFieldDomain,sourceFieldDecomposition,err,error,*999)
    CALL Domain_MeshComponentNumberGet(sourceFieldDomain,sourceMeshComponentNumber,err,error,*999)
    
    IF(ASSOCIATED(intermediateField)) THEN
      !Check the field has been finished
      CALL Field_AssertIsFinished(intermediateField,err,error,*999)
      !Check the user numbers match
      IF(intermediateFieldUserNumber/=intermediateField%userNumber) THEN
        localError="The specified intermediate field user number of "// &
          & TRIM(NumberToVString(intermediateFieldUserNumber,"*",err,error))// &
          & " does not match the user number of the specified intermediate field of "// &
          & TRIM(NumberToVString(intermediateField%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      !Check that the state field has the same region as the CellML environment
      NULLIFY(intermediateFieldRegion)
      CALL Field_RegionGet(intermediateField,intermediateFieldRegion,err,error,*999)
      IF(intermediateFieldRegion%userNumber/=region%userNumber) THEN
        localError="Invalid region setup. The specified intermediate field has been created on region"// &
          & " number "//TRIM(NumberToVString(intermediateFieldRegion%userNumber,"*",err,error))// &
          & " and the specified CellML environment has been created on region number "// &
          & TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      !Check the specified intermediate field has the same geometric field as the source field
      NULLIFY(intermediateGeometricField)
      CALL Field_GeometricFieldGet(intermediateField,intermediateGeometricField,err,error,*999)
      IF(.NOT.ASSOCIATED(sourceGeometricField,intermediateGeometricField)) THEN
        CALL FlagError("The specified intermediate field does not have the same geometric field as the "// &
          & "geometric field for the specified CellML environment.",err,error,*999)
      ENDIF
      !Check the specified intermediate field has the same decomposition as the source field
      NULLIFY(intermediateFieldDecomposition)
      CALL Field_DecompositionGet(intermediateField,intermediateFieldDecomposition,err,error,*999)
      IF(.NOT.ASSOCIATED(sourceFieldDecomposition,intermediateFieldDecomposition)) THEN
        CALL FlagError("The specified intermediate field does not have the same decomposition as the source "// &
          & "domain decomposition for the specified CellML environment.",err,error,*999)
      ENDIF
    ELSE
      !Check the user number has not already been used for a field in this region.
      NULLIFY(field)
      CALL Field_UserNumberFind(intermediateFieldUserNumber,region,field,err,error,*999)
      IF(ASSOCIATED(field)) THEN
        localError="The specified intermediate field user number of "// &
          & TRIM(NumberToVString(intermediateFieldUserNumber,"*",err,error))// &
          & "has already been used to create a field on region number "// &
          & TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDIF
    CALL CellML_IntermediateFieldInitialise(cellML,err,error,*999)
    NULLIFY(cellMLIntermediateField)
    CALL CellML_CellMLIntermediateFieldGet(cellMl,cellMLIntermediateField,err,error,*999)
    IF(ASSOCIATED(intermediateField)) THEN
      !Now check the supplied field.
      CALL Field_DataTypeCheck(intermediateField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
      CALL Field_TypeCheck(intermediateField,FIELD_GENERAL_TYPE,err,error,*999)
      CALL Field_NumberOfVariablesCheck(intermediateField,1,err,error,*999)
      CALL Field_VariableTypesCheck(intermediateField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
      CALL Field_NumberOfComponentsCheck(intermediateField,FIELD_U_VARIABLE_TYPE,cellML%maximumNumberOfIntermediate,err,error,*999)
      DO componentIdx=1,cellML%maximumNumberOfIntermediate
        CALL Field_ComponentMeshComponentCheck(intermediateField,FIELD_U_VARIABLE_TYPE,componentIdx,sourceMeshComponentNumber, &
          & err,error,*999)
        CALL Field_ComponentInterpolationCheck(intermediateField,FIELD_U_VARIABLE_TYPE,componentIdx, &
          & cellMLFieldMaps%sourceFieldInterpolationType,err,error,*999)
      ENDDO !componentIdx
    ELSE
      cellMLIntermediateField%intermediateFieldAutoCreated=.TRUE.
      !Create the CellML environment intermediate field
      CALL Field_CreateStart(intermediateFieldUserNumber,region,cellMLIntermediateField%intermediateField,err,error,*999)
      CALL Field_DataTypeSetAndLock(cellMLIntermediateField%intermediateField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
      CALL Field_LabelSet(cellMLIntermediateField%intermediateField,"CellMLIntermediateField",err,error,*999)
      CALL Field_TypeSetAndLock(cellMLIntermediateField%intermediateField,FIELD_GENERAL_TYPE,err,error,*999)
      CALL Field_DecompositionSetAndLock(cellMLIntermediateField%intermediateField,sourceFieldDecomposition,err,error,*999)
      CALL Field_GeometricFieldSetAndLock(cellMLIntermediateField%intermediateField,sourceGeometricField,err,error,*999)
      CALL Field_NumberOfVariablesSetAndLock(cellMLIntermediateField%intermediateField,1,err,error,*999)
      CALL Field_VariableTypesSetAndLock(cellMLIntermediateField%intermediateField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
      CALL Field_VariableLabelSet(cellMLIntermediateField%intermediateField,FIELD_U_VARIABLE_TYPE,"IntermediateVariable", &
        & err,error,*999)
      CALL Field_DOFOrderTypeSet(cellMLIntermediateField%intermediateField,FIELD_U_VARIABLE_TYPE, &
        & FIELD_CONTIGUOUS_COMPONENT_DOF_ORDER,err,error,*999)
      CALL Field_NumberOfComponentsSetAndLock(cellMLIntermediateField%intermediateField,FIELD_U_VARIABLE_TYPE, &
        & cellML%maximumNumberOfIntermediate,err,error,*999)
      DO componentIdx=1,cellML%maximumNumberOfIntermediate
        CALL Field_ComponentMeshComponentSetAndLock(cellMLIntermediateField%intermediateField,FIELD_U_VARIABLE_TYPE, &
          & componentIdx,sourceMeshComponentNumber,err,error,*999)
        CALL Field_ComponentInterpolationSetAndLock(cellMLIntermediateField%intermediateField,FIELD_U_VARIABLE_TYPE, &
          & componentIdx,cellMLFieldMaps%sourceFieldInterpolationType,err,error,*999)
      ENDDO !componentIdx
    ENDIF
    !Set pointers
    IF(cellMLIntermediateField%intermediateFieldAutoCreated) THEN            
      intermediateField=>cellMLIntermediateField%intermediateField
    ELSE
      cellMLIntermediateField%intermediateField=>intermediateField
    ENDIF
   
#else

    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*999)

#endif

    EXITS("CellML_IntermediateFieldCreateStart")
    RETURN
999 ERRORSEXITS("CellML_IntermediateFieldCreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_IntermediateFieldCreateStart

  !
  !=================================================================================================================================
  !

  !>Finialse the creation of the intermediate field for the given CellML environment.
  !! Finish creating the intermediate variable field for the provided CellML environment.
  !! - check for valid variables being defined?
  !! - maybe delete the field if no valid variables/components exist?
  SUBROUTINE CellML_IntermediateFieldCreateFinish(cellML,err,error,*)
    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<The CellML environment object for which to finalise the creation of the intermediate field.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local variables
    TYPE(CellMLIntermediateFieldType), POINTER :: cellMLIntermediateField
    TYPE(CellMLModelsFieldType), POINTER :: cellMLModelsField
 
    ENTERS("CellML_IntermediateFieldCreateFinish",err,error,*999)

#ifdef WITH_CELLML

    NULLIFY(cellMLIntermediateField)
    CALL CellML_CellMLIntermediateFieldGet(cellML,cellMLIntermediateField,err,error,*999)
    CALL CellMLIntermediateField_AssertNotFinished(cellMLIntermediateField,err,error,*999)
    NULLIFY(cellMLModelsField)
    CALL CellML_CellMLModelsFieldGet(cellML,cellMLModelsField,err,error,*999)
    CALL CellMLModelsField_AssertIsFinished(cellMLModelsField,err,error,*999)
    CALL CellMLModelsField_Check(cellML%modelsField,err,error,*999)
    !Finish the intermediate field creation
    IF(cellMLIntermediateField%intermediateFieldAutoCreated) &
      & CALL Field_CreateFinish(cellMLIntermediateField%intermediateField,err,error,*999)
    !As the intermediate field is strictly output do not initialise the values.
    cellMLIntermediateField%intermediateFieldFinished=.TRUE.

#else

    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*999)

#endif

    EXITS("CellML_IntermediateFieldCreateFinish")
    RETURN
999 ERRORSEXITS("CellML_IntermediateFieldCreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_IntermediateFieldCreateFinish

 !
  !=================================================================================================================================
  !

  !>Finalise a CellML environment models field and deallocate all memory.
  SUBROUTINE CellML_IntermediateFieldFinalise(cellMLIntermediateField,err,error,*)
    !Argument variables
    TYPE(CellMLIntermediateFieldType), POINTER :: cellMLIntermediateField !<A pointer to the CellML environment intermediate field to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !< The error string
    !Local variables
    
    ENTERS("CellML_IntermediateFieldFinalise",err,error,*999)

#ifdef WITH_CELLML

    IF(ASSOCIATED(cellMLIntermediateField)) THEN
      DEALLOCATE(cellMLIntermediateField)
    ENDIF

#else

    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*999)

#endif

    EXITS("CellML_IntermediateFieldFinalise")
    RETURN
999 ERRORSEXITS("CellML_IntermediateFieldFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_IntermediateFieldFinalise

  !
  !=================================================================================================================================
  !

  !>Initialise a CellML environment intermediate field.
  SUBROUTINE CellML_IntermediateFieldInitialise(cellML,err,error,*)
    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML environment to initialise the intermediate field for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !< The error string
    !Local variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("CellML_IntermediateFieldInitialise",err,error,*998)

#ifdef WITH_CELLML

    IF(.NOT.ASSOCIATED(cellML)) CALL FlagError("CellML environment is not associated.",err,error,*998)
    IF(ASSOCIATED(cellML%intermediateField)) &
      & CALL FlagError("CellML environment intermediate field is already associated.",err,error,*998)
    
    ALLOCATE(cellML%intermediateField,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate CellML environment intermediate field.",err,error,*999)
    cellML%intermediateField%cellML=>cellML
    cellML%intermediateField%intermediateFieldFinished=.FALSE.
    cellML%intermediateField%intermediateFieldAutoCreated=.FALSE.
    NULLIFY(cellML%intermediateField%intermediateField)

#else

    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*998)

#endif

    EXITS("CellML_IntermediateFieldInitialise")
    RETURN
999 CALL CellML_IntermediateFieldFinalise(cellML%intermediateField,dummyErr,dummyError,*998)
998 ERRORSEXITS("CellML_IntermediateFieldInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_IntermediateFieldInitialise

  !
  !=================================================================================================================================
  !

  !>Start the creation of the parameters field for the given CellML environment.
  SUBROUTINE CellML_ParametersFieldCreateStart(parametersFieldUserNumber,CELLML,parametersField,err,error,*)
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: parametersFieldUserNumber !<The unique identifier for the parameters field to be created for the given CellML environment object.
    TYPE(CellMLType), POINTER :: cellML !<The CellML environment object for which we will be defining the parameters field.
    TYPE(FieldType), POINTER :: parametersField  !<If associated on entry, a pointer to the user created parameters field which has the same user number as the specified parameters field user number. If not associated on entry, on exit, a pointer to the created parameters field for the CellML environment.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local variables
    INTEGER(INTG) :: componentIdx,sourceMeshComponentNumber
    TYPE(CellMLFieldMapsType), POINTER :: cellMLFieldMaps
    TYPE(CellMLParametersFieldType), POINTER :: cellMLParametersField
    TYPE(DecompositionType), POINTER :: parametersFieldDecomposition,sourceFieldDecomposition
    TYPE(DomainType), POINTER :: sourceFieldDomain
    TYPE(FieldType), POINTER :: field,parametersGeometricField,sourceGeometricField
    TYPE(RegionType), POINTER :: region,parametersFieldRegion
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("CellML_ParametersFieldCreateStart",err,error,*999)

#ifdef WITH_CELLML

    CALL CellML_AssertIsFinished(cellML,err,error,*999)
    NULLIFY(cellMLFieldMaps)
    CALL CellML_CellMLFieldMapsGet(cellML,cellMLFieldMaps,err,error,*999)
    CALL CellMLFieldMaps_AssertIsFinished(cellMLFieldMaps,err,error,*999)    
    IF(ASSOCIATED(cellML%parametersField)) &
      & CALL FlagError("The CellML environment parameters field is already associated.",err,error,*999)
    
    NULLIFY(region)
    CALL CellML_RegionGet(cellML,region,err,error,*999)
    NULLIFY(sourceGeometricField)
    CALL CellMLFieldMaps_SourceGeometricFieldGet(cellMLFieldMaps,sourceGeometricField,err,error,*999)
    NULLIFY(sourceFieldDomain)
    CALL CellMLFieldMaps_SourceFieldDomainGet(cellMLFieldMaps,sourceFieldDomain,err,error,*999)
    NULLIFY(sourceFieldDecomposition)
    CALL Domain_DecompositionGet(sourceFieldDomain,sourceFieldDecomposition,err,error,*999)
    CALL Domain_MeshComponentNumberGet(sourceFieldDomain,sourceMeshComponentNumber,err,error,*999)
     
    IF(ASSOCIATED(parametersField)) THEN
      !Check the field has been finished
      CALL Field_AssertIsFinished(parametersField,err,error,*999)
      !Check the user numbers match
      IF(parametersFieldUserNumber/=parametersField%userNumber) THEN
        localError="The specified parameters field user number of "// &
          & TRIM(NumberToVString(parametersFieldUserNumber,"*",err,error))// &
          & " does not match the user number of the specified parameters field of "// &
          & TRIM(NumberToVString(parametersField%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      !Check the field is defined on the same region as the CellML region
      NULLIFY(parametersFieldRegion)
      CALL Field_RegionGet(parametersField,parametersFieldRegion,err,error,*999)
      IF(parametersFieldRegion%userNumber/=region%userNumber) THEN
        localError="Invalid region setup. The specified parameters field has been created on region number "// &
          & TRIM(NumberToVString(parametersFieldRegion%userNumber,"*",err,error))// &
          & " and the specified CellML environment has been created on region number "// &
          & TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      !Check the specified parameters field has the same geometric field as the source field
      NULLIFY(parametersGeometricField)
      CALL Field_GeometricFieldGet(parametersField,parametersGeometricField,err,error,*999)
      IF(.NOT.ASSOCIATED(sourceGeometricField,parametersGeometricField)) THEN
        CALL FlagError("The specified parameters field does not have the same geometric field as the "// &
          & "geometric field for the specified CellML environment.",err,error,*999)
      ENDIF
      !Check the specified parameters field has the same decomposition as the source field
      NULLIFY(parametersFieldDecomposition)
      CALL Field_DecompositionGet(parametersField,parametersFieldDecomposition,err,error,*999)
      IF(.NOT.ASSOCIATED(sourceFieldDecomposition,parametersFieldDecomposition)) THEN
        CALL FlagError("The specified parameters field does not have the same decomposition as the source "// &
          & "domain decomposition for the specified CellML environment.",err,error,*999)
      ENDIF
    ELSE
      !Check the user number has not already been used for a field in this region.
      NULLIFY(field)
      CALL Field_UserNumberFind(parametersFieldUserNumber,region,field,err,error,*999)
      IF(ASSOCIATED(field)) THEN
        localError="The specified parameters field user number of "// &
          & TRIM(NumberToVString(parametersFieldUserNumber,"*",err,error))// &
          & "has already been used to create a field on region number "// &
          & TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDIF
    CALL CellML_ParametersFieldInitialise(cellML,err,error,*999)
    NULLIFY(cellMLParametersField)
    CALL CellML_CellMLParametersFieldGet(cellML,cellMLParametersField,err,error,*999)
    IF(ASSOCIATED(parametersField)) THEN
      !Now check the supplied field.
      CALL Field_DataTypeCheck(parametersField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
      CALL Field_TypeCheck(parametersField,FIELD_GENERAL_TYPE,err,error,*999)
      CALL Field_NumberOfVariablesCheck(parametersField,1,err,error,*999)
      CALL Field_VariableTypesCheck(parametersField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
      CALL Field_NumberOfComponentsCheck(parametersField,FIELD_U_VARIABLE_TYPE,cellML%maximumNumberOfParameters,err,error,*999)
      DO componentIdx=1,cellML%maximumNumberOfParameters
        CALL Field_ComponentMeshComponentCheck(parametersField,FIELD_U_VARIABLE_TYPE,componentIdx,sourceMeshComponentNumber, &
          & err,error,*999)
        CALL Field_ComponentInterpolationCheck(parametersField,FIELD_U_VARIABLE_TYPE,componentIdx, &
          & cellMLFieldMaps%sourceFieldInterpolationType,err,error,*999)
      ENDDO !componentIdx
    ELSE
      cellMLParametersField%parametersFieldAutoCreated=.TRUE.
      !Create the CellML environment parameters field
      CALL Field_CreateStart(parametersFieldUserNumber,region,cellMLParametersField%parametersField,err,error,*999)
      CALL Field_DataTypeSetAndLock(cellMLParametersField%parametersField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
      CALL Field_LabelSet(cellMLParametersField%parametersField,"CellMLParametersField",err,error,*999)
      CALL Field_TypeSetAndLock(cellMLParametersField%parametersField,FIELD_GENERAL_TYPE,err,error,*999)
      CALL Field_DecompositionSetAndLock(cellMLParametersField%parametersField,sourceFieldDecomposition,err,error,*999)
      CALL Field_GeometricFieldSetAndLock(cellMLParametersField%parametersField,sourceGeometricField,err,error,*999)
      CALL Field_NumberOfVariablesSetAndLock(cellMLParametersField%parametersField,1,err,error,*999)
      CALL Field_VariableTypesSetAndLock(cellMLParametersField%parametersField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
      CALL Field_VariableLabelSet(cellMLParametersField%parametersField,FIELD_U_VARIABLE_TYPE,"ParametersVariable",err,error,*999)
      CALL Field_DOFOrderTypeSet(cellMLParametersField%parametersField,FIELD_U_VARIABLE_TYPE, &
        & FIELD_CONTIGUOUS_COMPONENT_DOF_ORDER,err,error,*999)
      CALL Field_NumberOfComponentsSetAndLock(cellMLParametersField%parametersField,FIELD_U_VARIABLE_TYPE, &
        & cellML%maximumNumberOfParameters,err,error,*999)
      DO componentIdx=1,cellML%maximumNumberOfParameters
        CALL Field_ComponentMeshComponentSetAndLock(cellMLParametersField%parametersField,FIELD_U_VARIABLE_TYPE,componentIdx, &
          & sourceMeshComponentNumber,err,error,*999)
        CALL Field_ComponentInterpolationSetAndLock(cellMLParametersField%parametersField,FIELD_U_VARIABLE_TYPE,componentIdx, &
          & cellMLFieldMaps%sourceFieldInterpolationType,err,error,*999)
      ENDDO !componentIdx
    ENDIF
    !Set pointers
    IF(cellMLParametersField%parametersFieldAutoCreated) THEN            
      parametersField=>cellMLParametersField%parametersField
    ELSE
      cellMLParametersField%parametersField=>parametersField
    ENDIF

#else

    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*999)

#endif

    EXITS("CellML_ParametersFieldCreateStart")
    RETURN
999 ERRORSEXITS("CellML_ParametersFieldCreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_ParametersFieldCreateStart

  !
  !=================================================================================================================================
  !

  !>Finish the creation of the parameters field for the given CellML environment.
  SUBROUTINE CellML_ParametersFieldCreateFinish(cellML,err,error,*)
    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<The CellML environment object for which to finalise the parameters field.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local variables
    INTEGER(INTG) :: cellMLVariableType,errorCode,modelIdx,modelsDOFIdx,parameterComponentIdx
    INTEGER(INTG), POINTER :: modelsData(:)
    REAL(DP) :: initialValue
    TYPE(CellMLFieldMapsType), POINTER :: cellMLFieldMaps
    TYPE(CellMLModelType), POINTER :: cellMLModel
    TYPE(CellMLModelsFieldType), POINTER :: cellMLModelsField
    TYPE(CellMLParametersFieldType), POINTER :: cellMLParametersField
    TYPE(FieldVariableType), POINTER :: modelsVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("CellML_ParametersFieldCreateFinish",err,error,*999)

#ifdef WITH_CELLML

    NULLIFY(cellMLParametersField)
    CALL CellML_CellMLParametersFieldGet(cellML,cellMLParametersField,err,error,*999)
    CALL CellMLParametersField_AssertNotFinished(cellMLParametersField,err,error,*999)
    NULLIFY(cellMLModelsField)
    CALL CellML_CellMLModelsFieldGet(cellML,cellMLModelsField,err,error,*999)
    CALL CellMLModelsField_AssertIsFinished(cellMLModelsField,err,error,*999)
    CALL CellMLModelsField_Check(cellMLModelsField,err,error,*999)
    !Finish the parameters field creation
    IF(cellMLParametersField%parametersFieldAutoCreated) &
      & CALL Field_CreateFinish(cellMLParametersField%parametersField,err,error,*999)
    IF(cellMLModelsField%onlyOneModelIndex/=0) THEN
      !There are CellML models on this rank
      IF(cellMLModelsField%onlyOneModelIndex/=CELLML_MODELS_FIELD_NOT_CONSTANT) THEN
        !Only one model so optimise
        NULLIFY(cellMLModel)
        CALL CellML_CellMLModelGet(cellML,cellMLModelsField%onlyOneModelIndex,cellMLModel,err,error,*999)
        DO parameterComponentIdx=1,cellMLModel%numberOfParameters
          cellMLVariableType=CellML_MapCellMLFieldTypeToVariableType(CELLML_PARAMETERS_FIELD,err,error)
          IF(err/=0) GOTO 999
          errorCode = CELLML_MODEL_DEFINITION_GET_INITIAL_VALUE_BY_INDEX(cellMLModel%ptr,cellMLVariableType,&
            & parameterComponentIdx,initialValue)
          IF(errorCode/=0) THEN
            !problem getting the initial value
            localError="Error "//TRIM(NumberToVString(errorCode,"*",err,error))// &
              & ". Failed to get an initial value for the parameter variable with component number "//&
              & TRIM(NumberToVString(parameterComponentIdx,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          CALL Field_ComponentValuesInitialise(cellMLParametersField%parametersField,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,parameterComponentIdx,initialValue,err,error,*999)
        ENDDO !parameterComponentIdx
      ELSE
        !Multiple models so go through each dof.
        NULLIFY(cellMLFieldMaps)
        CALL CellML_CellMLFieldMapsGet(cellML,cellMLFieldMaps,err,error,*999)
        NULLIFY(modelsVariable)
        CALL Field_VariableGet(cellMLModelsField%modelsField,FIELD_U_VARIABLE_TYPE,modelsVariable,err,error,*999)
        NULLIFY(modelsData)
        CALL FieldVariable_ParameterSetDataGet(modelsVariable,FIELD_VALUES_SET_TYPE,modelsData,err,error,*999)
        DO modelsDOFIdx=1,modelsVariable%numberOfDofs
          modelIdx=modelsData(modelsDOFIdx)
          IF(modelIdx/=0) THEN
            NULLIFY(cellMLModel)
            CALL CellML_CellMLModelGet(cellML,modelIdx,cellMLModel,err,error,*999)
            DO parameterComponentIdx=1,cellMLModel%numberOfParameters
              cellMLVariableType=CellML_MapCellMLFieldTypeToVariableType(CELLML_PARAMETERS_FIELD,err,error)
              IF(err/=0) GOTO 999
              errorCode = CELLML_MODEL_DEFINITION_GET_INITIAL_VALUE_BY_INDEX(cellMLModel%ptr,cellMLVariableType,&
                & parameterComponentIdx,initialValue)
              IF(errorCode/=0) THEN
                !problem getting the initial value
                localError="Error "//TRIM(NumberToVString(errorCode,"*",err,error))// &
                  & ". Failed to get an initial value for the parameter variable with component number "//&
                  & TRIM(NumberToVString(parameterComponentIdx,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              ENDIF
              CALL CellML_FieldModelDofSet(modelsVariable,modelsDOFIdx,cellMLParametersField%parametersField, &
                & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,parameterComponentIdx,initialValue,err,error,*999)
            ENDDO !parameterComponentIdx
          ENDIF
        ENDDO !modelsDOFIdx
        CALL FieldVariable_ParameterSetDataRestore(modelsVariable,FIELD_VALUES_SET_TYPE,modelsData,err,error,*999)
      ENDIF
    ENDIF
    cellMLParametersField%parametersFieldFinished=.TRUE.
    
#else

    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*999)

#endif

    EXITS("CellML_ParametersFieldCreateFinish")
    RETURN
999 ERRORSEXITS("CellML_ParametersFieldCreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_ParametersFieldCreateFinish

  !
  !=================================================================================================================================
  !

  !>Finalise a CellML environment parameters field and deallocate all memory.
  SUBROUTINE CellML_ParametersFieldFinalise(cellMLParametersField,err,error,*)
    !Argument variables
    TYPE(CellMLParametersFieldType), POINTER :: cellMLParametersField !<A pointer to the CellML environment parameters field to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !< The error string
    !Local variables
    
    ENTERS("CellML_ParametersFieldFinalise",err,error,*999)

#ifdef WITH_CELLML

    IF(ASSOCIATED(cellMLParametersField)) THEN
      DEALLOCATE(cellMLParametersField)
    ENDIF

#else

    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*999)

#endif

    EXITS("CellML_ParametersFieldFinalise")
    RETURN
999 ERRORSEXITS("CellML_ParametersFieldFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_ParametersFieldFinalise

  !
  !=================================================================================================================================
  !

  !>Initialise a CellML environment parameters field.
  SUBROUTINE CellML_ParametersFieldInitialise(cellML,err,error,*)
    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML environment to initialise the parameters field for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !< The error string
    !Local variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("CellML_ParametersFieldInitialise",err,error,*998)

#ifdef WITH_CELLML

    IF(ASSOCIATED(cellML)) CALL FlagError("CellML environment is not associated.",err,error,*998)
    IF(ASSOCIATED(cellML%parametersField)) &
      & CALL FlagError("CellML environment parameters field is already associated.",err,error,*998)
      
    ALLOCATE(cellML%parametersField,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate CellML environment parameters field.",err,error,*999)
    cellML%parametersField%cellML=>cellML
    cellML%parametersField%parametersFieldFinished=.FALSE.
    cellML%parametersField%parametersFieldAutoCreated=.FALSE.
    NULLIFY(cellML%parametersField%parametersField)

#else

    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*998)

#endif

    EXITS("CellML_ParametersFieldInitialise")
    RETURN
999 CALL CellML_ParametersFieldFinalise(cellML%parametersField,dummyErr,dummyError,*998)
998 ERRORSEXITS("CellML_ParametersFieldInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_ParametersFieldInitialise

  !
  !=================================================================================================================================
  !

  !>Validate and instantiate the specified CellML environment.
  !! Users should call this routine once they have set up the CellML environment to allow the CellML environment to be validated and instantiated into computable code.
  !! - generate and compile code.
  !! - distribute things amonst processors?
  !! - create the internal fields?
  !! - does this method need to get called or can it be done implicitly as needed by other solution procedures?
  SUBROUTINE CellML_Generate(cellML,err,error,*)
    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<The CellML environment object for which to generate computable code.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    !Local variables
 
    ENTERS("CellML_Generate",err,error,*999)

#ifdef WITH_CELLML

    CALL CellML_AssertIsFinished(cellML,err,error,*999)
    
    !Set the generated flag
    cellML%cellMLGenerated=.TRUE.

#else

    CALL FlagError("Must compile with WITH_CELLML ON to use CellML functionality.",err,error,*999)

#endif

    EXITS("CellML_Generate")
    RETURN
999 ERRORSEXITS("CellML_Generate",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_Generate

  !
  !=================================================================================================================================
  !

  !>Finalises the CellML environments and deallocates all memory.
  SUBROUTINE CellMLEnvironments_Finalise(cellMLEnvironments,err,error,*)

    !Argument variables
    TYPE(CellMLEnvironmentsType), POINTER :: cellMLEnvironments !<A pointer to the CellML environments to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: cellMLIdx

    ENTERS("CellMLEnvironments_Finalise",err,error,*999)

    IF(ASSOCIATED(cellMLEnvironments)) THEN
      IF(ALLOCATED(cellMLEnvironments%environments)) THEN
        DO cellMLIdx=1,SIZE(cellMLEnvironments%environments,1)
          CALL CellML_Finalise(cellMLEnvironments%environments(cellMLIdx)%ptr,err,error,*999)
        ENDDO !cellMLIdx
        DEALLOCATE(cellMLEnvironments%environments)
      ENDIF
      DEALLOCATE(cellMLEnvironments)
    ENDIF

    EXITS("CellMLEnvironments_Finalise")
    RETURN
999 ERRORSEXITS("CellMLEnvironments_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE CellMLEnvironments_Finalise
  
  !
  !=================================================================================================================================
  !

  !>Initialises the CellML environments.
  SUBROUTINE CellMLEnvironments_Initialise(region,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: region  !<A pointer to the region to initialise the CellML environments for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("CellMLEnvironments_Initialise",err,error,*998)

    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated.",err,error,*998)
    IF(ASSOCIATED(region%cellMLEnvironments)) CALL FlagError("Region CellML environments is already associated.",err,error,*998)
      
    ALLOCATE(region%cellMLEnvironments,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate region CellML environments.",err,error,*999)
    region%cellMLEnvironments%region=>region
    region%cellMLEnvironments%numberOfEnvironments=0

    EXITS("CellMLEnvironments_Initialise")
    RETURN
999 CALL CellMLEnvironments_Finalise(region%cellMLEnvironments,dummyErr,dummyError,*998)
998 ERRORSEXITS("CellMLEnvironments_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE CellMLEnvironments_Initialise

  !
  !================================================================================================================================
  !

  !> Maps a CellML variable type to a CellML field type
  FUNCTION CellML_MapCellMLVariableToFieldType(cellMLVariableType,err,error)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: cellMLVariableType !<The CellML variable type to map. \see CellML_FieldTypes,CmissCellML
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    INTEGER(INTG) :: CellML_MapCellMLVariableToFieldType
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ! CellML variable types in OpenCMISS(cellml)/CellMLModelDefinitionF.h:
    !   state=1; known=2; wanted=3.
    !   independent=4 - but untested and maybe not working yet.
    ENTERS("CellML_MapCellMLVariableToFieldType",err,error,*999)
    
    SELECT CASE(cellMLVariableType)
    CASE(1)
      CellML_MapCellMLVariableToFieldType=CELLML_STATE_FIELD
    CASE(2)
      CellML_MapCellMLVariableToFieldType=CELLML_PARAMETERS_FIELD
    CASE(3)
      CellML_MapCellMLVariableToFieldType=CELLML_INTERMEDIATE_FIELD
    CASE(4)
      localError="CellML variable type "//TRIM(NumberToVString(cellMLVariableType,"*",err,error))//&
      & " (independent variable) support not yet implemented"
      CALL FlagError(localError,err,error,*999)
    CASE DEFAULT
      localError="CellML variable type "//TRIM(NumberToVString(cellMLVariableType,"*",err,error))//&
      & " is invalid or not implemented."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("CellML_MapCellMLVariableToFieldType")
    RETURN
999 ERRORSEXITS("CellML_MapCellMLVariableToFieldType",err,error)
    RETURN

  END FUNCTION CellML_MapCellMLVariableToFieldType

  !
  !================================================================================================================================
  !

  !> Maps a CellML field type to a CellML variable type
  FUNCTION CellML_MapCellMLFieldTypeToVariableType(cellMLFieldType,err,error)
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: cellMLFieldType !<The CellML field type to map. \see CellML_FieldTypes,CmissCellML
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    INTEGER(INTG) :: CellML_MapCellMLFieldTypeToVariableType
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ! CellML variable types in OpenCMISS(cellml)/CellMLModelDefinitionF.h:
    !   state=1; known=2; wanted=3.
    !   independent=4 - but untested and maybe not working yet.
    ENTERS("CellML_MapCellMLFieldTypeToVariableType",err,error,*999)
    SELECT CASE(cellMLFieldType)
    CASE(CELLML_STATE_FIELD)
      CellML_MapCellMLFieldTypeToVariableType=1
    CASE(CELLML_PARAMETERS_FIELD)
      CellML_MapCellMLFieldTypeToVariableType=2
    CASE(CELLML_INTERMEDIATE_FIELD)
      CellML_MapCellMLFieldTypeToVariableType=3
    !CASE(4)
    !  localError="CellML variable type "//TRIM(NumberToVString(cellMLVariableType,"*",err,error))//&
    !  & " (independent variable) support not yet implemented"
    !  CALL FlagError(localError,err,error,*999)
    CASE DEFAULT
      localError="CellML field type "//TRIM(NumberToVString(cellMLFieldType,"*",err,error))//&
      & " is invalid or not implemented"
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("CellML_MapCellMLFieldTypeToVariableType")
    RETURN
999 ERRORSEXITS("CellML_MapCellMLFieldTypeToVariableType",err,error)
    RETURN

  END FUNCTION CellML_MapCellMLFieldTypeToVariableType

END MODULE CmissCellML
