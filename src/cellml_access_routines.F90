!> \file
!> \author Chris Bradley
!> \brief This module contains all CellML access method routines.
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

!> \addtogroup OpenCMISS_CellML OpenCMISS::Iron::CellML
!> This module contains all CellML access method routines.
MODULE CellMLAccessRoutines
  
  USE BaseRoutines
  USE Kinds
  USE ISO_VARYING_STRING
  USE Strings
  USE Types
  
#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup CellML_Types OpenCMISS::Iron::CellML::Constants
  !> \brief CellML constants
  !>@{
  !> \addtogroup CellML_FieldTypes OpenCMISS::Iron::CellML::Constants::FieldTypes
  !> \brief CellML field parameter types
  !> \see CmissCellML,OpenCMISS_CellMLFieldTypes
  !> CellML model variables being mapped to or from fields will have an initial type of CellML_UNKNOWN_FIELD. This will be set to
  !> the appropriate type once the model is instantiated and the type can be correctly determined.
  !> \todo Link to appropriate methods for instantiation.
  !>@{
  INTEGER(INTG), PARAMETER :: CELLML_MODELS_FIELD = 1 !<CellML models field \see CellML_FieldTypes,CmissCellML
  INTEGER(INTG), PARAMETER :: CELLML_STATE_FIELD = 2 !<CellML state field \see CellML_FieldTypes,CmissCellML
  INTEGER(INTG), PARAMETER :: CELLML_INTERMEDIATE_FIELD = 3 !<CellML intermediate field \see CellML_FieldTypes,CmissCellML
  INTEGER(INTG), PARAMETER :: CELLML_PARAMETERS_FIELD = 4 !<CellML parameters field \see CellML_FieldTypes,CmissCellML
  !>@}

  !> \addtogroup CellML_FieldMappingTypes OpenCMISS::Iron::CellML::Constants::FieldMappingTypes
  !> \brief CellML field parameter types
  !> \see CmissCellML,OpenCMISS_CellMLFieldTypes
  !>@{
  INTEGER(INTG), PARAMETER :: CELLML_MAP_TO_FIELD_TYPE = 1 !<A CellML to field mapping type \see CellML_FieldMappingTypes,CmissCellML
  INTEGER(INTG), PARAMETER :: CELLML_MAP_FROM_FIELD_TYPE = 2 !<A field to CellML mapping type \see CellML_FieldMappingTypes,CmissCellML
  !>@}

  !> \addtogroup CellML_ModelsFieldTypes OpenCMISS::Iron::CellML::Constants::ModelsFieldTypes
  !> \brief CellML field parameter types
  !> \see CmissCellML,OpenCMISS_CellMLFieldTypes
  !>@{
  INTEGER(INTG), PARAMETER :: CELLML_MODELS_FIELD_NOT_CHECKED = -2 !<The CellML environment models field has not been checked. \see CellML_ModelsFieldTypes,CmissCellML
  INTEGER(INTG), PARAMETER :: CELLML_MODELS_FIELD_NOT_CONSTANT =-1 !<The CellML environement models field is not constant. \see CellML_ModelsFieldTypes,CmissCellML
  !>@}
  !>@}

  !Module types

  !Module variables

  !Interfaces

  PUBLIC CELLML_MODELS_FIELD,CELLML_STATE_FIELD,CELLML_INTERMEDIATE_FIELD,CELLML_PARAMETERS_FIELD

  PUBLIC CELLML_MAP_TO_FIELD_TYPE,CELLML_MAP_FROM_FIELD_TYPE
  
  PUBLIC CELLML_MODELS_FIELD_NOT_CHECKED,CELLML_MODELS_FIELD_NOT_CONSTANT

  PUBLIC CellML_AssertIsFinished,CellML_AssertNotFinished

  PUBLIC CellML_CellMLEnvironmentsGet

  PUBLIC CellML_CellMLFieldMapsGet

  PUBLIC CellML_CellMLIntermediateFieldExists

  PUBLIC CellML_CellMLIntermediateFieldGet

  PUBLIC CellML_CellMLModelGet
  
  PUBLIC CellML_CellMLModelsFieldExists

  PUBLIC CellML_CellMLModelsFieldGet

  PUBLIC CellML_CellMLParametersFieldExists
  
  PUBLIC CellML_CellMLParametersFieldGet
  
  PUBLIC CellML_CellMLStateFieldExists

  PUBLIC CellML_CellMLStateFieldGet

  PUBLIC CellML_IntermediateFieldExists

  PUBLIC CellML_IntermediateFieldGet

  PUBLIC CellML_MaximumNumberOfIntermediateGet

  PUBLIC CellML_MaximumNumberOfParametersGet

  PUBLIC CellML_MaximumNumberOfStateGet

  PUBLIC CellML_ModelsFieldExists

  PUBLIC CellML_ModelsFieldGet

  PUBLIC CellML_NumberOfModelsGet

  PUBLIC CellML_ParametersFieldExists

  PUBLIC CellML_ParametersFieldGet

  PUBLIC CellML_RegionGet

  PUBLIC CellML_StateFieldExists
  
  PUBLIC CellML_StateFieldGet
  
  PUBLIC CellML_UserNumberFind

  PUBLIC CellMLFieldMaps_AssertIsFinished,CellMLFieldMaps_AssertNotFinished

  PUBLIC CellMLFieldMaps_CellMLModelMapsGet

  PUBLIC CellMLFieldMaps_SourceFieldDomainGet

  PUBLIC CellMLFieldMaps_SourceGeometricFieldGet
  
  PUBLIC CellMLIntermediateField_AssertIsFinished,CellMLIntermediateField_AssertNotFinished

  PUBLIC CellMLIntermediateField_CellMLGet

  PUBLIC CellMLIntermediateField_IntermediateFieldGet

  PUBLIC CellMLModel_CellMLGet
  
  PUBLIC CellMLModel_NumberOfIntermediateGet
  
  PUBLIC CellMLModel_NumberOfParametersGet
  
  PUBLIC CellMLModel_NumberOfStateGet
  
  PUBLIC CellMLModelMaps_CellMLModelFromMapGet
  
  PUBLIC CellMLModelMaps_CellMLModelToMapGet

  PUBLIC CellMLModelsField_AssertIsFinished,CellMLModelsField_AssertNotFinished

  PUBLIC CellMLModelsField_CellMLGet

  PUBLIC CellMLModelsField_ModelsFieldGet

  PUBLIC CellMLModelsField_OnlyOneModelIndexGet
  
  PUBLIC CellMLParametersField_AssertIsFinished,CellMLParametersField_AssertNotFinished

  PUBLIC CellMLParametersField_CellMLGet

  PUBLIC CellMLParametersField_ParametersFieldGet
  
  PUBLIC CellMLStateField_AssertIsFinished,CellMLStateField_AssertNotFinished

  PUBLIC CellMLStateField_CellMLGet

  PUBLIC CellMLStateField_StateFieldGet
  

CONTAINS

  !
  !=================================================================================================================================
  !

  !>Assert that a cellml has been finished
  SUBROUTINE CellML_AssertIsFinished(cellml,err,error,*)

    !Argument Variables
    TYPE(CellMLType), POINTER, INTENT(IN) :: cellml !<The cellml to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("CellML_AssertIsFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(cellml)) CALL FlagError("CellML is not associated.",err,error,*999)
#endif    

    IF(.NOT.cellml%cellmlFinished) THEN
      localError="CellML user number "//TRIM(NumberToVString(cellml%userNumber,"*",err,error))
      IF(ASSOCIATED(cellML%region)) &
        & localError=localError//" of region number "//TRIM(NumberToVString(cellML%region%userNumber,"*",err,error))
      localError=localError//" has not been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("CellML_AssertIsFinished")
    RETURN
999 ERRORSEXITS("CellML_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_AssertIsFinished

  !
  !=================================================================================================================================
  !

  !>Assert that a cellml has not been finished
  SUBROUTINE CellML_AssertNotFinished(cellml,err,error,*)

    !Argument Variables
    TYPE(CellMLType), POINTER, INTENT(IN) :: cellml !<The cellml to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("CellML_AssertNotFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(cellml)) CALL FlagError("CellML is not associated.",err,error,*999)
#endif    

    IF(cellml%cellmlFinished) THEN
      localError="CellML user number "//TRIM(NumberToVString(cellml%userNumber,"*",err,error))
      IF(ASSOCIATED(cellML%region)) &
        & localError=localError//" of region number "//TRIM(NumberToVString(cellML%region%userNumber,"*",err,error))
      localError=localError//" has already been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("CellML_AssertNotFinished")
    RETURN
999 ERRORSEXITS("CellML_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_AssertNotFinished
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the CellML environments for the specified CellML environment.
  SUBROUTINE CellML_CellMLEnvironmentsGet(cellML,cellMLEnvironments,err,error,*)

    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML to get the environments for
    TYPE(CellMLEnvironmentsType), POINTER :: cellMLEnvironments  !<On exit, a pointer to environments for the CellML environment. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("CellML_CellMLEnvironmentsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(cellMLEnvironments)) CALL FlagError("CellML environments is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cellML)) CALL FlagError("CellML is not associated.",err,error,*999)
#endif    
    
    cellMLEnvironments=>cellML%environments

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(cellMLEnvironments)) THEN
      localError="The environments for CellML user number "//TRIM(NumberToVString(cellML%userNumber,"*",err,error))
      IF(ASSOCIATED(cellML%region)) &
        & localError=localError//" of region number "//TRIM(NumberToVString(cellML%region%userNumber,"*",err,error))
      localError=localError//" is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    EXITS("CellML_CellMLEnvironmentsGet")
    RETURN
999 NULLIFY(cellMLEnvironments)
998 ERRORSEXITS("CellML_CellMLEnvironmentsGet",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_CellMLEnvironmentsGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the field maps for the specified CellML environment.
  SUBROUTINE CellML_CellMLFieldMapsGet(cellML,cellMLFieldMaps,err,error,*)

    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML to get the field maps for
    TYPE(CellMLFieldMapsType), POINTER :: cellMLFieldMaps  !<On exit, a pointer to field maps for the CellML environment. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("CellML_CellMLFieldMapsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(cellMLFieldMaps)) CALL FlagError("CellML field maps is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cellML)) CALL FlagError("CellML environment is not associated.",err,error,*999)
#endif    
    
    cellMLFieldMaps=>cellML%fieldMaps

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(cellMLFieldMaps)) THEN
      localError="The field maps for CellML user number "//TRIM(NumberToVString(cellML%userNumber,"*",err,error))
      IF(ASSOCIATED(cellML%region)) &
        & localError=localError//" of region number "//TRIM(NumberToVString(cellML%region%userNumber,"*",err,error))
      localError=localError//" are not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    EXITS("CellML_CellMLFieldMapsGet")
    RETURN
999 NULLIFY(cellMLFieldMaps)
998 ERRORSEXITS("CellML_CellMLFieldMapsGet",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_CellMLFieldMapsGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the intermediate field information if it exists for the specified CellML environment.
  SUBROUTINE CellML_CellMLIntermediateFieldExists(cellML,cellMLIntermediateField,err,error,*)

    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML to check the intermediate field information exists for
    TYPE(CellMLIntermediateFieldType), POINTER :: cellMLIntermediateField  !<On exit, a pointer to intermediate field information for the CellML environment if it exists. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("CellML_CellMLIntermediateFieldExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(cellMLIntermediateField)) CALL FlagError("CellML intermediate field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cellML)) CALL FlagError("CellML environment is not associated.",err,error,*999)
#endif    
    
    cellMLIntermediateField=>cellML%intermediateField

    EXITS("CellML_CellMLIntermediateFieldExists")
    RETURN
999 NULLIFY(cellMLIntermediateField)
998 ERRORSEXITS("CellML_CellMLIntermediateFieldExists",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_CellMLIntermediateFieldExists

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the intermediate field information for the specified CellML environment.
  SUBROUTINE CellML_CellMLIntermediateFieldGet(cellML,cellMLIntermediateField,err,error,*)

    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML to get the intermediate field information for
    TYPE(CellMLIntermediateFieldType), POINTER :: cellMLIntermediateField  !<On exit, a pointer to intermediate field information for the CellML environment. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("CellML_CellMLIntermediateFieldGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(cellMLIntermediateField)) CALL FlagError("CellML intermediate field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cellML)) CALL FlagError("CellML environment is not associated.",err,error,*999)
#endif    
    
    cellMLIntermediateField=>cellML%intermediateField

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(cellMLIntermediateField)) THEN
      localError="The intermediate field information for CellML user number "// &
        & TRIM(NumberToVString(cellML%userNumber,"*",err,error))
      IF(ASSOCIATED(cellML%region)) &
        & localError=localError//" of region number "//TRIM(NumberToVString(cellML%region%userNumber,"*",err,error))
      localError=localError//" are not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    EXITS("CellML_CellMLIntermediateFieldGet")
    RETURN
999 NULLIFY(cellMLIntermediateField)
998 ERRORSEXITS("CellML_CellMLIntermediateFieldGet",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_CellMLIntermediateFieldGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to a CellML model for the specified index in a CellML environment.
  SUBROUTINE CellML_CellMLModelGet(cellML,modelIndex,cellMLModel,err,error,*)

    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML to get the model for
    INTEGER(INTG), INTENT(IN) :: modelIndex !<The model index to get
    TYPE(CellMLModelType), POINTER :: cellMLModel !<On exit, a pointer to specified model in the CellML environment. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("CellML_CellMLModelGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(cellMLModel)) CALL FlagError("CellML model is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cellML)) CALL FlagError("CellML environment is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(cellML%models)) THEN
      localError="Models is not allocated for CellML user number "//TRIM(NumberToVString(cellML%userNumber,"*",err,error))
      IF(ASSOCIATED(cellML%region)) &
        & localError=localError//" of region number "//TRIM(NumberToVString(cellML%region%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(modelIndex<1.OR.modelIndex>cellML%numberOfModels) THEN
      localError="The specified model index of "//TRIM(NumberToVString(modelIndex,"*",err,error))// &
        & " for CellML user number "//TRIM(NumberToVString(cellML%userNumber,"*",err,error))
      IF(ASSOCIATED(cellML%region)) &
        & localError=localError//" of region number "//TRIM(NumberToVString(cellML%region%userNumber,"*",err,error))
      localError=localError//". The model index should be >= 1 and <= "// &
        & TRIM(NumberToVString(cellML%numberOfModels,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    cellMLModel=>cellML%models(modelIndex)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(cellMLModel)) THEN
      localError="The CellML model for model index "//TRIM(NumberToVString(modelIndex,"*",err,error))// &
        & " for CellML user number "//TRIM(NumberToVString(cellML%userNumber,"*",err,error))
      IF(ASSOCIATED(cellML%region)) &
        & localError=localError//" of region number "//TRIM(NumberToVString(cellML%region%userNumber,"*",err,error))
      localError=localError//" is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    EXITS("CellML_CellMLModelGet")
    RETURN
999 NULLIFY(cellMLModel)
998 ERRORSEXITS("CellML_CellMLModelGet",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_CellMLModelGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the models field information if it exists for the specified CellML environment.
  SUBROUTINE CellML_CellMLModelsFieldExists(cellML,cellMLModelsField,err,error,*)

    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML to check if the models field information exists
    TYPE(CellMLModelsFieldType), POINTER :: cellMLModelsField  !<On exit, a pointer to models field information for the CellML environment if it exists. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("CellML_CellMLModelsFieldExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(cellMLModelsField)) CALL FlagError("CellML models field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cellML)) CALL FlagError("CellML environment is not associated.",err,error,*999)
#endif    
    
    cellMLModelsField=>cellML%modelsField

    EXITS("CellML_CellMLModelsFieldExists")
    RETURN
999 NULLIFY(cellMLModelsField)
998 ERRORSEXITS("CellML_CellMLModelsFieldExists",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_CellMLModelsFieldExists

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the models field information for the specified CellML environment.
  SUBROUTINE CellML_CellMLModelsFieldGet(cellML,cellMLModelsField,err,error,*)

    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML to get the models field information for
    TYPE(CellMLModelsFieldType), POINTER :: cellMLModelsField  !<On exit, a pointer to models field information for the CellML environment. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("CellML_CellMLModelsFieldGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(cellMLModelsField)) CALL FlagError("CellML models field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cellML)) CALL FlagError("CellML environment is not associated.",err,error,*999)
#endif    
    
    cellMLModelsField=>cellML%modelsField

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(cellMLModelsField)) THEN
      localError="The models field information for CellML user number "//TRIM(NumberToVString(cellML%userNumber,"*",err,error))
      IF(ASSOCIATED(cellML%region)) &
        & localError=localError//" of region number "//TRIM(NumberToVString(cellML%region%userNumber,"*",err,error))
      localError=localError//" are not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    EXITS("CellML_CellMLModelsFieldGet")
    RETURN
999 NULLIFY(cellMLModelsField)
998 ERRORSEXITS("CellML_CellMLModelsFieldGet",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_CellMLModelsFieldGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the parameters field information if it exists for the specified CellML environment.
  SUBROUTINE CellML_CellMLParametersFieldExists(cellML,cellMLParametersField,err,error,*)

    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML to check the parameters field information exists for
    TYPE(CellMLParametersFieldType), POINTER :: cellMLParametersField  !<On exit, a pointer to parameters field information for the CellML environment if it exists. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("CellML_CellMLParametersFieldExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(cellMLParametersField)) CALL FlagError("CellML parameters field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cellML)) CALL FlagError("CellML environment is not associated.",err,error,*999)
#endif
    
    cellMLParametersField=>cellML%parametersField

    EXITS("CellML_CellMLParametersFieldExists")
    RETURN
999 NULLIFY(cellMLParametersField)
998 ERRORSEXITS("CellML_CellMLParametersFieldExists",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_CellMLParametersFieldExists

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the parameters field information for the specified CellML environment.
  SUBROUTINE CellML_CellMLParametersFieldGet(cellML,cellMLParametersField,err,error,*)

    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML to get the parameters field information for
    TYPE(CellMLParametersFieldType), POINTER :: cellMLParametersField  !<On exit, a pointer to parameters field information for the CellML environment. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("CellML_CellMLParametersFieldGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(cellMLParametersField)) CALL FlagError("CellML parameters field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cellML)) CALL FlagError("CellML environment is not associated.",err,error,*999)
#endif
    
    cellMLParametersField=>cellML%parametersField

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(cellMLParametersField)) THEN
      localError="The parameters field information for CellML user number "//TRIM(NumberToVString(cellML%userNumber,"*",err,error))
      IF(ASSOCIATED(cellML%region)) &
        & localError=localError//" of region number "//TRIM(NumberToVString(cellML%region%userNumber,"*",err,error))
      localError=localError//" are not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    EXITS("CellML_CellMLParametersFieldGet")
    RETURN
999 NULLIFY(cellMLParametersField)
998 ERRORSEXITS("CellML_CellMLParametersFieldGet",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_CellMLParametersFieldGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the state field information if it exists for the specified CellML environment.
  SUBROUTINE CellML_CellMLStateFieldExists(cellML,cellMLStateField,err,error,*)

    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML to check if the state field information exists for
    TYPE(CellMLStateFieldType), POINTER :: cellMLStateField  !<On exit, a pointer to state field information for the CellML environment if it exists. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("CellML_CellMLStateFieldExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(cellMLStateField)) CALL FlagError("CellML state field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cellML)) CALL FlagError("CellML environment is not associated.",err,error,*999)
#endif    
    
    cellMLStateField=>cellML%stateField

    EXITS("CellML_CellMLStateFieldExists")
    RETURN
999 NULLIFY(cellMLStateField)
998 ERRORSEXITS("CellML_CellMLStateFieldExists",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_CellMLStateFieldExists

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the state field information for the specified CellML environment.
  SUBROUTINE CellML_CellMLStateFieldGet(cellML,cellMLStateField,err,error,*)

    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML to get the state field information for
    TYPE(CellMLStateFieldType), POINTER :: cellMLStateField  !<On exit, a pointer to state field information for the CellML environment. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("CellML_CellMLStateFieldGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(cellMLStateField)) CALL FlagError("CellML state field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cellML)) CALL FlagError("CellML environment is not associated.",err,error,*999)
#endif    
    
    cellMLStateField=>cellML%stateField

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(cellMLStateField)) THEN
      localError="The state field information for CellML user number "//TRIM(NumberToVString(cellML%userNumber,"*",err,error))
      IF(ASSOCIATED(cellML%region)) &
        & localError=localError//" of region number "//TRIM(NumberToVString(cellML%region%userNumber,"*",err,error))
      localError=localError//" are not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    EXITS("CellML_CellMLStateFieldGet")
    RETURN
999 NULLIFY(cellMLStateField)
998 ERRORSEXITS("CellML_CellMLStateFieldGet",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_CellMLStateFieldGet

  !
  !================================================================================================================================
  !

  !>Checks if an intermediate field for the specified CellML environment exists.
  SUBROUTINE CellML_IntermediateFieldExists(cellML,intermediateField,err,error,*)

    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML to check the existance of the intermediate field for
    TYPE(FieldType), POINTER :: intermediateField  !<On exit, a pointer to intermediate field for the CellML environment if it exists. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("CellML_IntermediateFieldExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(intermediateField)) CALL FlagError("Intermediate field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cellML)) CALL FlagError("CellML environment is not associated.",err,error,*999)
#endif
    
    IF(ASSOCIATED(cellML%intermediateField)) THEN
      CALL CellMLIntermediateField_AssertIsFinished(cellML%intermediateField,err,error,*999)    
      intermediateField=>cellML%intermediateField%intermediateField
    ELSE
      NULLIFY(intermediateField)
    ENDIF

    EXITS("CellML_IntermediateFieldExists")
    RETURN
999 NULLIFY(intermediateField)
998 ERRORSEXITS("CellML_IntermediateFieldExists",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_IntermediateFieldExists

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the intermediate field for the specified CellML environment.
  SUBROUTINE CellML_IntermediateFieldGet(cellML,intermediateField,err,error,*)

    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML to get the intermediate field for
    TYPE(FieldType), POINTER :: intermediateField  !<On exit, a pointer to intermediate field for the CellML environment. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("CellML_IntermediateFieldGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(intermediateField)) CALL FlagError("Intermediate field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cellML)) CALL FlagError("CellML environment is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(cellML%intermediateField)) THEN
      localError="The CellML intermediate field information is not associated for CellML user number "// &
        & TRIM(NumberToVString(cellML%userNumber,"*",err,error))
      IF(ASSOCIATED(cellML%region)) &
        & localError=localError//" of region number "//TRIM(NumberToVString(cellML%region%userNumber,"*",err,error))
      localError=localError//" is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    CALL CellMLIntermediateField_AssertIsFinished(cellML%intermediateField,err,error,*999)
    
    intermediateField=>cellML%intermediateField%intermediateField

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(intermediateField)) THEN
      localError="The intermediate field is not associated for the intermediate field information for CellML user number "// &
        & TRIM(NumberToVString(cellML%userNumber,"*",err,error))
      IF(ASSOCIATED(cellML%region)) &
        & localError=localError//" of region number "//TRIM(NumberToVString(cellML%region%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    EXITS("CellML_IntermediateFieldGet")
    RETURN
999 NULLIFY(intermediateField)
998 ERRORSEXITS("CellML_IntermediateFieldGet",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_IntermediateFieldGet

  !
  !================================================================================================================================
  !

  !>Returns the maximum number of intermediate variables across all models for the specified CellML environment.
  SUBROUTINE CellML_MaximumNumberOfIntermediateGet(cellML,maximumNumberOfIntermediate,err,error,*)

    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML to get the maximum number of intermediate variables for
    INTEGER(INTG), INTENT(OUT) :: maximumNumberOfIntermediate  !<On exit, the maximum number of intermidate variables across all models for the CellML environment.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("CellML_MaximumNumberOfIntermediateGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(cellML)) CALL FlagError("CellML environment is not associated.",err,error,*999)
#endif    
    
    maximumNumberOfIntermediate=cellML%maximumNumberOfIntermediate

    EXITS("CellML_MaximumNumberOfIntermediateGet")
    RETURN
999 ERRORSEXITS("CellML_MaximumNumberOfIntermediateGet",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_MaximumNumberOfIntermediateGet

  !
  !================================================================================================================================
  !

  !>Returns the maximum number of parameters variables across all models for the specified CellML environment.
  SUBROUTINE CellML_MaximumNumberOfParametersGet(cellML,maximumNumberOfParameters,err,error,*)

    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML to get the maximum number of parameters variables for
    INTEGER(INTG), INTENT(OUT) :: maximumNumberOfParameters  !<On exit, the maximum number of intermidate variables across all models for the CellML environment.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("CellML_MaximumNumberOfParametersGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(cellML)) CALL FlagError("CellML environment is not associated.",err,error,*999)
#endif    
    
    maximumNumberOfParameters=cellML%maximumNumberOfParameters

    EXITS("CellML_MaximumNumberOfParametersGet")
    RETURN
999 ERRORSEXITS("CellML_MaximumNumberOfParametersGet",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_MaximumNumberOfParametersGet

  !
  !================================================================================================================================
  !

  !>Returns the maximum number of state variables across all models for the specified CellML environment.
  SUBROUTINE CellML_MaximumNumberOfStateGet(cellML,maximumNumberOfState,err,error,*)

    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML to get the maximum number of state variables for
    INTEGER(INTG), INTENT(OUT) :: maximumNumberOfState  !<On exit, the maximum number of intermidate variables across all models for the CellML environment.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("CellML_MaximumNumberOfStateGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(cellML)) CALL FlagError("CellML environment is not associated.",err,error,*999)
#endif    
    
    maximumNumberOfState=cellML%maximumNumberOfState

    EXITS("CellML_MaximumNumberOfStateGet")
    RETURN
999 ERRORSEXITS("CellML_MaximumNumberOfStateGet",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_MaximumNumberOfStateGet

  !
  !================================================================================================================================
  !

  !>Checks if the models field for the specified CellML environment exists.
  SUBROUTINE CellML_ModelsFieldExists(cellML,modelsField,err,error,*)

    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML to check the existance of the models field for
    TYPE(FieldType), POINTER :: modelsField  !<On exit, a pointer to models field for the CellML environment if it exists. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("CellML_ModelsFieldGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(modelsField)) CALL FlagError("Models field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cellML)) CALL FlagError("CellML environment is not associated.",err,error,*999)
#endif

    IF(ASSOCIATED(cellML%modelsField)) THEN
      CALL CellMLModelsField_AssertIsFinished(cellML%modelsField,err,error,*999)
      modelsField=>cellML%modelsField%modelsField
    ELSE
      NULLIFY(modelsField)
    ENDIF

    EXITS("CellML_ModelsFieldExists")
    RETURN
999 NULLIFY(modelsField)
998 ERRORSEXITS("CellML_ModelsFieldExists",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_ModelsFieldExists

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the models field for the specified CellML environment.
  SUBROUTINE CellML_ModelsFieldGet(cellML,modelsField,err,error,*)

    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML to get the models field for
    TYPE(FieldType), POINTER :: modelsField  !<On exit, a pointer to models field for the CellML environment. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("CellML_ModelsFieldGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(modelsField)) CALL FlagError("Models field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cellML)) CALL FlagError("CellML environment is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(cellML%modelsField)) THEN
      localError="The CellML models field information is not associated for CellML user number "// &
        & TRIM(NumberToVString(cellML%userNumber,"*",err,error))
      IF(ASSOCIATED(cellML%region)) &
        & localError=localError//" of region number "//TRIM(NumberToVString(cellML%region%userNumber,"*",err,error))
      localError=localError//" is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    CALL CellMLModelsField_AssertIsFinished(cellML%modelsField,err,error,*999)
    
    modelsField=>cellML%modelsField%modelsField

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(modelsField)) THEN
      localError="The models field is not associated for the models field information for CellML user number "// &
        & TRIM(NumberToVString(cellML%userNumber,"*",err,error))
      IF(ASSOCIATED(cellML%region)) &
        & localError=localError//" of region number "//TRIM(NumberToVString(cellML%region%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    EXITS("CellML_ModelsFieldGet")
    RETURN
999 NULLIFY(modelsField)
998 ERRORSEXITS("CellML_ModelsFieldGet",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_ModelsFieldGet

  !
  !================================================================================================================================
  !

  !>Returns the number of models for the specified CellML environment.
  SUBROUTINE CellML_NumberOfModelsGet(cellML,numberOfModels,err,error,*)

    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML to get the number of models for
    INTEGER(INTG), INTENT(OUT) :: numberOfModels  !<On exit, the number of models for the CellML environment.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("CellML_NumberOfModelsGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(cellML)) CALL FlagError("CellML environment is not associated.",err,error,*999)
#endif    
    
    numberOfModels=cellML%numberOfModels

    EXITS("CellML_NumberOfModelsGet")
    RETURN
999 ERRORSEXITS("CellML_NumberOfModelsGet",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_NumberOfModelsGet

  !
  !================================================================================================================================
  !

  !>Checks if a parameters field for the specified CellML environment exists.
  SUBROUTINE CellML_ParametersFieldExists(cellML,parametersField,err,error,*)

    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML to check the existance of the parameters field for
    TYPE(FieldType), POINTER :: parametersField  !<On exit, a pointer to parameters field for the CellML environment if it exists. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("CellML_ParametersFieldExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(parametersField)) CALL FlagError("Parameters field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cellML)) CALL FlagError("CellML environment is not associated.",err,error,*999)
#endif

    IF(ASSOCIATED(cellML%parametersField)) THEN
      CALL CellMLParametersField_AssertIsFinished(cellML%parametersField,err,error,*999)    
      parametersField=>cellML%parametersField%parametersField
    ELSE
      NULLIFY(parametersField)
    ENDIF

    EXITS("CellML_ParametersFieldExists")
    RETURN
999 NULLIFY(parametersField)
998 ERRORSEXITS("CellML_ParametersFieldExists",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_ParametersFieldExists

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the parameters field for the specified CellML environment.
  SUBROUTINE CellML_ParametersFieldGet(cellML,parametersField,err,error,*)

    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML to get the parameters field for
    TYPE(FieldType), POINTER :: parametersField  !<On exit, a pointer to parameters field for the CellML environment. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("CellML_ParametersFieldGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(parametersField)) CALL FlagError("Parameters field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cellML)) CALL FlagError("CellML environment is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(cellML%parametersField)) THEN
      localError="The CellML parameters field information is not associated for CellML user number "// &
        & TRIM(NumberToVString(cellML%userNumber,"*",err,error))
      IF(ASSOCIATED(cellML%region)) &
        & localError=localError//" of region number "//TRIM(NumberToVString(cellML%region%userNumber,"*",err,error))
      localError=localError//" is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    CALL CellMLParametersField_AssertIsFinished(cellML%parametersField,err,error,*999)
    
    parametersField=>cellML%parametersField%parametersField

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(parametersField)) THEN
      localError="The parameters field is not associated for the parameters field information for CellML user number "// &
        & TRIM(NumberToVString(cellML%userNumber,"*",err,error))
      IF(ASSOCIATED(cellML%region)) &
        & localError=localError//" of region number "//TRIM(NumberToVString(cellML%region%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    EXITS("CellML_ParametersFieldGet")
    RETURN
999 NULLIFY(parametersField)
998 ERRORSEXITS("CellML_ParametersFieldGet",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_ParametersFieldGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the region for the specified CellML environment.
  SUBROUTINE CellML_RegionGet(cellML,region,err,error,*)

    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML to get the region for
    TYPE(RegionType), POINTER :: region  !<On exit, a pointer to region for the CellML environment. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("CellML_RegionGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(region)) CALL FlagError("Region is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cellML)) CALL FlagError("CellML environment is not associated.",err,error,*999)
#endif    
    
    region=>cellML%region

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(region)) THEN
      localError="The region for CellML user number "//TRIM(NumberToVString(cellML%userNumber,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    EXITS("CellML_RegionGet")
    RETURN
999 NULLIFY(region)
998 ERRORSEXITS("CellML_RegionGet",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_RegionGet

  !
  !================================================================================================================================
  !

  !>Checks if the state field for the specified CellML environment exists.
  SUBROUTINE CellML_StateFieldExists(cellML,stateField,err,error,*)

    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML to check the existance of the state field for
    TYPE(FieldType), POINTER :: stateField  !<On exit, a pointer to state field for the CellML environment if it exists. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("CellML_StateFieldExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(stateField)) CALL FlagError("State field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cellML)) CALL FlagError("CellML environment is not associated.",err,error,*999)
#endif

    IF(ASSOCIATED(cellML%stateField)) THEN
      CALL CellMLStateField_AssertIsFinished(cellML%stateField,err,error,*999)    
      stateField=>cellML%stateField%stateField
    ELSE
      NULLIFY(stateField)
    ENDIF

    EXITS("CellML_StateFieldExists")
    RETURN
999 NULLIFY(stateField)
998 ERRORSEXITS("CellML_StateFieldExists",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_StateFieldExists

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the state field for the specified CellML environment.
  SUBROUTINE CellML_StateFieldGet(cellML,stateField,err,error,*)

    !Argument variables
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML to get the state field for
    TYPE(FieldType), POINTER :: stateField  !<On exit, a pointer to state field for the CellML environment. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("CellML_StateFieldGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(stateField)) CALL FlagError("State field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cellML)) CALL FlagError("CellML environment is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(cellML%stateField)) THEN
      localError="The CellML state field information is not associated for CellML user number "// &
        & TRIM(NumberToVString(cellML%userNumber,"*",err,error))
      IF(ASSOCIATED(cellML%region)) &
        & localError=localError//" of region number "//TRIM(NumberToVString(cellML%region%userNumber,"*",err,error))
      localError=localError//" is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    CALL CellMLStateField_AssertIsFinished(cellML%stateField,err,error,*999)
    
    stateField=>cellML%stateField%stateField

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(stateField)) THEN
      localError="The state field is not associated for the state field information for CellML user number "// &
        & TRIM(NumberToVString(cellML%userNumber,"*",err,error))
      IF(ASSOCIATED(cellML%region)) &
        & localError=localError//" of region number "//TRIM(NumberToVString(cellML%region%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    EXITS("CellML_StateFieldGet")
    RETURN
999 NULLIFY(stateField)
998 ERRORSEXITS("CellML_StateFieldGet",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_StateFieldGet

  !
  !=================================================================================================================================
  !

  !>Finds and returns a pointer to the CellML environment identified by a user number on a region. If no CellML environment with that user number exists cellml is left nullified.
  SUBROUTINE CellML_UserNumberFind(userNumber,region,cellml,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number to find.
    TYPE(RegionType), POINTER :: region !<A pointer to the region to find the CellML user number.
    TYPE(CellMLType), POINTER :: cellml !<On return a pointer to the CellML environment with the given user number. If no CellML environment with that user number exists then the pointer is returned as NULL. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: cellmlIdx
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("CellML_UserNumberFind",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated.",err,error,*999)
    IF(ASSOCIATED(cellml)) CALL FlagError("CellML is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(region%cellMLEnvironments)) CALL FlagError("Region CellML environments is not associated.",err,error,*999)
#endif    

    NULLIFY(cellml)
    IF(ALLOCATED(region%cellMLEnvironments%environments)) THEN
      DO cellmlIdx=1,region%cellMLEnvironments%numberOfEnvironments
#ifdef WITH_PRECHECKS        
        IF(.NOT.ASSOCIATED(region%cellMLEnvironments%environments(cellmlIdx)%ptr)) THEN
          localError="The CellML pointer in the CellML environments is not associated for CellML index "// &
            & TRIM(NumberToVString(cellmlIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
#endif        
        IF(region%cellMLEnvironments%environments(cellmlIdx)%ptr%userNumber==userNumber) THEN
          cellml=>region%cellMLEnvironments%environments(cellmlIdx)%ptr
          EXIT
        ENDIF
      ENDDO !cellmlIdx
    ENDIF

    EXITS("CellML_UserNumberFind")
    RETURN
999 ERRORSEXITS("CellML_UserNumberFind",err,error)
    RETURN 1
    
  END SUBROUTINE CellML_UserNumberFind

  !
  !=================================================================================================================================
  !

  !>Assert that a cellml field maps has been finished
  SUBROUTINE CellMLFieldMaps_AssertIsFinished(cellmlFieldMaps,err,error,*)

    !Argument Variables
    TYPE(CellMLFieldMapsType), POINTER, INTENT(IN) :: cellmlFieldMaps !<The cellml field map to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("CellMLFieldMaps_AssertIsFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(cellmlFieldMaps)) CALL FlagError("CellML field maps is not associated.",err,error,*999)
#endif    

    IF(.NOT.cellmlFieldMaps%cellmlFieldMapsFinished) THEN
      localError="The CellML field maps"
      IF(ASSOCIATED(cellMLFieldMaps%cellML)) THEN
        localError=localError//" for CellML user number "//TRIM(NumberToVString(cellMLFieldMaps%cellML%userNumber,"*",err,error))
        IF(ASSOCIATED(cellMLFieldMaps%cellML%region)) &
          & localError=localError//" of region number "// &
          & TRIM(NumberToVString(cellMLFieldMaps%cellML%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//" has not been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("CellMLFieldMaps_AssertIsFinished")
    RETURN
999 ERRORSEXITS("CellMLFieldMaps_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE CellMLFieldMaps_AssertIsFinished

  !
  !=================================================================================================================================
  !

  !>Assert that a cellml field maps has not been finished
  SUBROUTINE CellMLFieldMaps_AssertNotFinished(cellMLFieldMaps,err,error,*)

    !Argument Variables
    TYPE(CellMLFieldMapsType), POINTER, INTENT(IN) :: cellmlFieldMaps !<The cellml field maps to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("CellMLFieldMaps_AssertNotFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(cellMLFieldMaps)) CALL FlagError("CellML field maps is not associated.",err,error,*999)
#endif    

    IF(cellMLFieldMaps%cellMLFieldMapsFinished) THEN
      localError="The CellML field maps"
      IF(ASSOCIATED(cellMLFieldMaps%cellML)) THEN
        localError=localError//" for CellML user number "//TRIM(NumberToVString(cellMLFieldMaps%cellML%userNumber,"*",err,error))
        IF(ASSOCIATED(cellMLFieldMaps%cellML%region)) &
          & localError=localError//" of region number "// &
          & TRIM(NumberToVString(cellMLFieldMaps%cellML%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//" has already been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("CellMLFieldMaps_AssertNotFinished")
    RETURN
999 ERRORSEXITS("CellMLFieldMaps_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE CellMLFieldMaps_AssertNotFinished
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the model maps for the specified model index in the CellML field maps.
  SUBROUTINE CellMLFieldMaps_CellMLModelMapsGet(cellMLFieldMaps,modelIndex,cellMLModelMaps,err,error,*)

    !Argument variables
    TYPE(CellMLFieldMapsType), POINTER :: cellMLFieldMaps !<A pointer to the CellML field maps to get the model maps for
    INTEGER(INTG), INTENT(IN) :: modelIndex !<The model index of the model maps to get
    TYPE(CellMLModelMapsType), POINTER :: cellMLModelMaps  !<On exit, a pointer to model maps for the CellML field maps. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("CellMLFieldMaps_CellMLModelMapsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(cellMLModelMaps)) CALL FlagError("CellML model maps is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cellMLFieldMaps)) CALL FlagError("CellML field maps is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(cellMLFieldMaps%modelMaps)) THEN
      localError="The CellML field maps models maps"
      IF(ASSOCIATED(cellMLFieldMaps%cellML)) THEN
        localError=localError//" for CellML user number "//TRIM(NumberToVString(cellMLFieldMaps%cellML%userNumber,"*",err,error))
        IF(ASSOCIATED(cellMLFieldMaps%cellML%region)) &
          & localError=localError//" of region number "// &
          & TRIM(NumberToVString(cellMLFieldMaps%cellML%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//" are not allocated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(modelIndex<1.OR.modelIndex>SIZE(cellMLFieldMaps%modelMaps,1)) THEN
      localError="The specified model index of "//TRIM(NumberToVString(modelIndex,"*",err,error))// &
        & " is invalid for the CellML field maps"
      IF(ASSOCIATED(cellMLFieldMaps%cellML)) THEN
        localError=localError//" for CellML user number "//TRIM(NumberToVString(cellMLFieldMaps%cellML%userNumber,"*",err,error))
        IF(ASSOCIATED(cellMLFieldMaps%cellML%region)) &
          & localError=localError//" of region number "// &
          & TRIM(NumberToVString(cellMLFieldMaps%cellML%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//". The model index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(cellMLFieldMaps%modelMaps,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    cellMLModelMaps=>cellMLFieldMaps%modelMaps(modelIndex)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(cellMLModelMaps)) THEN
      localError="The model maps for model index "//TRIM(NumberToVString(modelIndex,"*",err,error))// &
        & " of the CellML field maps"
      IF(ASSOCIATED(cellMLFieldMaps%cellML)) THEN
        localError=localError//" for CellML user number "//TRIM(NumberToVString(cellMLFieldMaps%cellML%userNumber,"*",err,error))
        IF(ASSOCIATED(cellMLFieldMaps%cellML%region)) &
          & localError=localError//" of region number "// &
          & TRIM(NumberToVString(cellMLFieldMaps%cellML%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//" is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    EXITS("CellMLFieldMaps_CellMLModelMapsGet")
    RETURN
999 NULLIFY(cellMLModelMaps)
998 ERRORSEXITS("CellMLFieldMaps_CellMLModelMapsGet",err,error)
    RETURN 1
    
  END SUBROUTINE CellMLFieldMaps_CellMLModelMapsGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the source field domain in the CellML field maps.
  SUBROUTINE CellMLFieldMaps_SourceFieldDomainGet(cellMLFieldMaps,sourceFieldDomain,err,error,*)

    !Argument variables
    TYPE(CellMLFieldMapsType), POINTER :: cellMLFieldMaps !<A pointer to the CellML field maps to get the source field domain for
    TYPE(DomainType), POINTER :: sourceFieldDomain  !<On exit, a pointer to source field domain for the CellML field maps. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif
    
    ENTERS("CellMLFieldMaps_SourceFieldDomainGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(sourceFieldDomain)) CALL FlagError("Source field domain is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cellMLFieldMaps)) CALL FlagError("CellML field maps is not associated.",err,error,*999)
#endif    
    
    sourceFieldDomain=>cellMLFieldMaps%sourceFieldDomain

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(sourceFieldDomain)) THEN
      localError="The source field domain of the CellML field maps"
      IF(ASSOCIATED(cellMLFieldMaps%cellML)) THEN
        localError=localError//" for CellML user number "//TRIM(NumberToVString(cellMLFieldMaps%cellML%userNumber,"*",err,error))
        IF(ASSOCIATED(cellMLFieldMaps%cellML%region)) &
          & localError=localError//" of region number "// &
          & TRIM(NumberToVString(cellMLFieldMaps%cellML%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//" is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    EXITS("CellMLFieldMaps_SourceFieldDomainGet")
    RETURN
999 NULLIFY(sourceFieldDomain)
998 ERRORSEXITS("CellMLFieldMaps_SourceFieldDomainGet",err,error)
    RETURN 1
    
  END SUBROUTINE CellMLFieldMaps_SourceFieldDomainGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the source geometric field in the CellML field maps.
  SUBROUTINE CellMLFieldMaps_SourceGeometricFieldGet(cellMLFieldMaps,sourceGeometricField,err,error,*)

    !Argument variables
    TYPE(CellMLFieldMapsType), POINTER :: cellMLFieldMaps !<A pointer to the CellML field maps to get the source geometric field for
    TYPE(FieldType), POINTER :: sourceGeometricField  !<On exit, a pointer to source geometric field for the CellML field maps. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("CellMLFieldMaps_SourceGeometricFieldGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(sourceGeometricField)) CALL FlagError("Source geometric field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cellMLFieldMaps)) CALL FlagError("CellML field maps is not associated.",err,error,*999)
#endif    
    
    sourceGeometricFIeld=>cellMLFieldMaps%sourceGeometricField

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(sourceGeometricField)) THEN
      localError="The source geometric field of the CellML field maps"
      IF(ASSOCIATED(cellMLFieldMaps%cellML)) THEN
        localError=localError//" for CellML user number "//TRIM(NumberToVString(cellMLFieldMaps%cellML%userNumber,"*",err,error))
        IF(ASSOCIATED(cellMLFieldMaps%cellML%region)) &
          & localError=localError//" of region number "// &
          & TRIM(NumberToVString(cellMLFieldMaps%cellML%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//" is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    EXITS("CellMLFieldMaps_SourceGeometricFieldGet")
    RETURN
999 NULLIFY(sourceGeometricField)
998 ERRORSEXITS("CellMLFieldMaps_SourceGeometricFieldGet",err,error)
    RETURN 1
    
  END SUBROUTINE CellMLFieldMaps_SourceGeometricFieldGet

  !
  !=================================================================================================================================
  !

  !>Assert that a CellML intermediate field has been finished
  SUBROUTINE CellMLIntermediateField_AssertIsFinished(cellMLIntermediateField,err,error,*)

    !Argument Variables
    TYPE(CellMLIntermediateFieldType), POINTER, INTENT(IN) :: cellMLIntermediateField !<The CellML intermediate field to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("CellMLIntermediateField_AssertIsFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(cellMLIntermediateField)) CALL FlagError("CellML intermediate field is not associated.",err,error,*999)
#endif    

    IF(.NOT.cellMLIntermediateField%intermediateFieldFinished) THEN
      localError="The CellML intermediate field"
      IF(ASSOCIATED(cellMLIntermediateField%cellML)) THEN
        localError=localError//" for CellML user number "// &
          & TRIM(NumberToVString(cellMLIntermediateField%cellML%userNumber,"*",err,error))
        IF(ASSOCIATED(cellMLIntermediateField%cellML%region)) &
          & localError=localError//" of region number "// &
          & TRIM(NumberToVString(cellMLIntermediateField%cellML%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//" has not been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("CellMLIntermediateField_AssertIsFinished")
    RETURN
999 ERRORSEXITS("CellMLIntermediateField_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE CellMLIntermediateField_AssertIsFinished

  !
  !=================================================================================================================================
  !

  !>Assert that a CellML intermediate field has not been finished
  SUBROUTINE CellMLIntermediateField_AssertNotFinished(cellMLIntermediateField,err,error,*)

    !Argument Variables
    TYPE(CellMLIntermediateFieldType), POINTER, INTENT(IN) :: cellMLIntermediateField !<The CellML intermediate field to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("CellMLIntermediateField_AssertNotFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(cellMLIntermediateField)) CALL FlagError("CellML intermediate field is not associated.",err,error,*999)
#endif    

    IF(cellMLIntermediateField%intermediateFieldFinished) THEN
      localError="The CellML intermediate field"
      IF(ASSOCIATED(cellMLIntermediateField%cellML)) THEN
        localError=localError//" for CellML user number "// &
          & TRIM(NumberToVString(cellMLIntermediateField%cellML%userNumber,"*",err,error))
        IF(ASSOCIATED(cellMLIntermediateField%cellML%region)) &
          & localError=localError//" of region number "// &
          & TRIM(NumberToVString(cellMLIntermediateField%cellML%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//" has already been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("CellMLIntermediateField_AssertNotFinished")
    RETURN
999 ERRORSEXITS("CellMLIntermediateField_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE CellMLIntermediateField_AssertNotFinished
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the CellML environment for the specified CellML intermediate field information.
  SUBROUTINE CellMLIntermediateField_CellMLGet(cellMLIntermediateField,cellML,err,error,*)

    !Argument variables
    TYPE(CellMLIntermediateFieldType), POINTER :: cellMLIntermediateField  !<A pointer to CellML intermediate field information to get the CellML environment from
    TYPE(CellMLType), POINTER :: cellML !<On return, the CellML environment. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("CellMLIntermediateField_CellMLGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(cellML)) CALL FlagError("CellML environment is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cellMLIntermediateField)) CALL FlagError("CellML intermediate field is not associated.",err,error,*999)
#endif    
    
    cellML=>cellMLIntermediateField%cellML

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(cellML)) &
      & CALL FlagError("The CellML environment for the intermediate field information is not associated.",err,error,*999)
#endif    

    EXITS("CellMLIntermediateField_CellMLGet")
    RETURN
999 NULLIFY(cellML)
998 ERRORSEXITS("CellMLIntermediateField_CellMLGet",err,error)
    RETURN 1
    
  END SUBROUTINE CellMLIntermediateField_CellMLGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the intermediate field for the specified CellML intermediate field information.
  SUBROUTINE CellMLIntermediateField_IntermediateFieldGet(cellMLIntermediateField,intermediateField,err,error,*)

    !Argument variables
    TYPE(CellMLIntermediateFieldType), POINTER :: cellMLIntermediateField  !<A pointer to CellML intermediate field information to get the intermediate field from
    TYPE(FieldType), POINTER :: intermediateField !<On return, the intermediate field. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("CellMLIntermediateField_IntermediateFieldGet",err,error,*998)

#ifdef WITH_PRECHECKS
    IF(ASSOCIATED(intermediateField)) CALL FlagError("Intermediate field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cellMLIntermediateField)) CALL FlagError("CellML intermediate field is not associated.",err,error,*999)
#endif    
    
    intermediateField=>cellMLIntermediateField%intermediateField

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(intermediateField)) THEN
      localError="The intermediate field for the intermediate field information"
      IF(ASSOCIATED(cellMLIntermediateField%cellML)) THEN
        localError=localError//" for CellML user number "// &
          & TRIM(NumberToVString(cellMLIntermediateField%cellML%userNumber,"*",err,error))
        IF(ASSOCIATED(cellMLIntermediateField%cellML%region)) &
          & localError=localError//" of region number "// &
          & TRIM(NumberToVString(cellMLIntermediateField%cellML%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//" is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    EXITS("CellMLIntermediateField_IntermediateFieldGet")
    RETURN
999 NULLIFY(intermediateField)
998 ERRORS("CellMLIntermediateField_IntermediateFieldGet",err,error)
    EXITS("CellMLIntermediateField_IntermediateFieldGet")
    RETURN 1
    
  END SUBROUTINE CellMLIntermediateField_IntermediateFieldGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the CellML environment for the specified CellML model.
  SUBROUTINE CellMLModel_CellMLGet(cellMLModel,cellML,err,error,*)

    !Argument variables
    TYPE(CellMLModelType), POINTER :: cellMLModel  !<A pointer to CellML model to get the CellML environment from
    TYPE(CellMLType), POINTER :: cellML !<On return, the CellML environment for the CellML model. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("CellMLModel_CellMLGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(cellML)) CALL FlagError("CellML environment is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cellMLModel)) CALL FlagError("CellML model is not associated.",err,error,*999)
#endif    
    
    cellML=>cellMLModel%cellML

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(cellML)) &
      & CALL FlagError("The CellML environment for the CellML model is not associated.",err,error,*999)
#endif    

    EXITS("CellMLModel_CellMLGet")
    RETURN
999 NULLIFY(cellML)
998 ERRORSEXITS("CellMLModel_CellMLGet",err,error)
    RETURN 1
    
  END SUBROUTINE CellMLModel_CellMLGet

  !
  !================================================================================================================================
  !

  !>Returns the number of intermediate variables for the specified CellML model.
  SUBROUTINE CellMLModel_NumberOfIntermediateGet(cellMLModel,numberOfIntermediate,err,error,*)

    !Argument variables
    TYPE(CellMLModelType), POINTER :: cellMLModel  !<A pointer to CellML model to get the number of intermediate variables for
    INTEGER(INTG), INTENT(OUT) :: numberOfIntermediate !<On return, the number of intermediate variables for the CellML model.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("CellMLModel_NumberOfIntermediateGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(cellMLModel)) CALL FlagError("CellML model is not associated.",err,error,*999)
#endif    
    
    numberOfIntermediate=cellMLModel%numberOfIntermediate

    EXITS("CellMLModel_NumberOfIntermediateGet")
    RETURN
999 ERRORSEXITS("CellMLModel_NumberOfIntermediateGet",err,error)
    RETURN 1
    
  END SUBROUTINE CellMLModel_NumberOfIntermediateGet

  !
  !================================================================================================================================
  !

  !>Returns the number of parameters variables for the specified CellML model.
  SUBROUTINE CellMLModel_NumberOfParametersGet(cellMLModel,numberOfParameters,err,error,*)

    !Argument variables
    TYPE(CellMLModelType), POINTER :: cellMLModel  !<A pointer to CellML model to get the number of parameter variables for
    INTEGER(INTG), INTENT(OUT) :: numberOfParameters !<On return, the number of parameters variables for the CellML model.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("CellMLModel_NumberOfParametersGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(cellMLModel)) CALL FlagError("CellML model is not associated.",err,error,*999)
#endif    
    
    numberOfParameters=cellMLModel%numberOfParameters

    EXITS("CellMLModel_NumberOfParametersGet")
    RETURN
999 ERRORSEXITS("CellMLModel_NumberOfParametersGet",err,error)
    RETURN 1
    
  END SUBROUTINE CellMLModel_NumberOfParametersGet

  !
  !================================================================================================================================
  !

  !>Returns the number of state variables for the specified CellML model.
  SUBROUTINE CellMLModel_NumberOfStateGet(cellMLModel,numberOfState,err,error,*)

    !Argument variables
    TYPE(CellMLModelType), POINTER :: cellMLModel  !<A pointer to CellML model to get the number of state variables for
    INTEGER(INTG), INTENT(OUT) :: numberOfState !<On return, the number of state variables for the CellML model.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("CellMLModel_NumberOfStateGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(cellMLModel)) CALL FlagError("CellML model is not associated.",err,error,*999)
#endif    
    
    numberOfState=cellMLModel%numberOfState

    EXITS("CellMLModel_NumberOfStateGet")
    RETURN
999 ERRORSEXITS("CellMLModel_NumberOfStateGet",err,error)
    RETURN 1
    
  END SUBROUTINE CellMLModel_NumberOfStateGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the from model map for the specified from map index in the CellML model maps.
  SUBROUTINE CellMLModelMaps_CellMLModelFromMapGet(cellMLModelMaps,fromMapIndex,cellMLModelFromMap,err,error,*)

    !Argument variables
    TYPE(CellMLModelMapsType), POINTER :: cellMLModelMaps !<A pointer to the CellML model maps to get the from model map for
    INTEGER(INTG), INTENT(IN) :: fromMapIndex !<The from map index of the model from map to get
    TYPE(CellMLModelMapType), POINTER :: cellMLModelFromMap  !<On exit, a pointer to model from map for the CellML models maps. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("CellMLModelMaps_CellMLModelFromMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(cellMLModelFromMap)) CALL FlagError("CellML model from map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cellMLModelMaps)) CALL FlagError("CellML model maps is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(cellMLModelMaps%fieldsMappedFrom)) THEN
      localError="The CellML fields mapped from model map"
      IF(ASSOCIATED(cellMLModelMaps%cellMLFieldMaps)) THEN
        IF(ASSOCIATED(cellMLModelMaps%cellMLFieldMaps%cellML)) THEN
          localError=localError//" for CellML user number "// &
            & TRIM(NumberToVString(cellMLModelMaps%cellMLFieldMaps%cellML%userNumber,"*",err,error))
          IF(ASSOCIATED(cellMLModelMaps%cellMLFieldMaps%cellML%region)) &
            & localError=localError//" of region number "// &
            & TRIM(NumberToVString(cellMLModelMaps%cellMLFieldMaps%cellML%region%userNumber,"*",err,error))
        ENDIF
      ENDIF
      localError=localError//" is not allocated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(fromMapIndex<1.OR.fromMapIndex>cellMLModelMaps%numberOfFieldsMappedFrom) THEN
      localError="The specified from map index of "//TRIM(NumberToVString(fromMapIndex,"*",err,error))// &
        & " is invalid for the model maps of the CellML field maps"
      IF(ASSOCIATED(cellMLModelMaps%cellMLFieldMaps)) THEN
        IF(ASSOCIATED(cellMLModelMaps%cellMLFieldMaps%cellML)) THEN
          localError=localError//" for CellML user number "// &
            & TRIM(NumberToVString(cellMLModelMaps%cellMLFieldMaps%cellML%userNumber,"*",err,error))
          IF(ASSOCIATED(cellMLModelMaps%cellMLFieldMaps%cellML%region)) &
            & localError=localError//" of region number "// &
            & TRIM(NumberToVString(cellMLModelMaps%cellMLFieldMaps%cellML%region%userNumber,"*",err,error))
        ENDIF
      ENDIF
      localError=localError//". The from map index should be >= 1 and <= "// &
        & TRIM(NumberToVString(cellMLModelMaps%numberOfFieldsMappedFrom,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    cellMLModelFromMap=>cellMLModelMaps%fieldsMappedFrom(fromMapIndex)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(cellMLModelFromMap)) THEN
      localError="The from model map for from map index "//TRIM(NumberToVString(fromMapIndex,"*",err,error))// &
        & " of the model maps of the CellML field maps"
      IF(ASSOCIATED(cellMLModelMaps%cellMLFieldMaps)) THEN
        IF(ASSOCIATED(cellMLModelMaps%cellMLFieldMaps%cellML)) THEN
          localError=localError//" for CellML user number "// &
            & TRIM(NumberToVString(cellMLModelMaps%cellMLFieldMaps%cellML%userNumber,"*",err,error))
          IF(ASSOCIATED(cellMLModelMaps%cellMLFieldMaps%cellML%region)) &
            & localError=localError//" of region number "// &
            & TRIM(NumberToVString(cellMLModelMaps%cellMLFieldMaps%cellML%region%userNumber,"*",err,error))
        ENDIF
      ENDIF
      localError=localError//" is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    EXITS("CellMLModelMaps_CellMLModelFromMapGet")
    RETURN
999 NULLIFY(cellMLModelFromMap)
998 ERRORSEXITS("CellMLModelMaps_CellMLModelFromMapGet",err,error)
    RETURN 1
    
  END SUBROUTINE CellMLModelMaps_CellMLModelFromMapGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the to model map for the specified to map index in the CellML model maps.
  SUBROUTINE CellMLModelMaps_CellMLModelToMapGet(cellMLModelMaps,toMapIndex,cellMLModelToMap,err,error,*)

    !Argument variables
    TYPE(CellMLModelMapsType), POINTER :: cellMLModelMaps !<A pointer to the CellML model maps to get the to model map for
    INTEGER(INTG), INTENT(IN) :: toMapIndex !<The to map index of the model to map to get
    TYPE(CellMLModelMapType), POINTER :: cellMLModelToMap !<On exit, a pointer to model to map for the CellML models maps. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("CellMLModelMaps_CellMLModelToMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(cellMLModelToMap)) CALL FlagError("CellML model to map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cellMLModelMaps)) CALL FlagError("CellML model maps is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(cellMLModelMaps%fieldsMappedTo)) THEN
      localError="The CellML fields mapped to model map"
      IF(ASSOCIATED(cellMLModelMaps%cellMLFieldMaps)) THEN
        IF(ASSOCIATED(cellMLModelMaps%cellMLFieldMaps%cellML)) THEN
          localError=localError//" for CellML user number "// &
            & TRIM(NumberToVString(cellMLModelMaps%cellMLFieldMaps%cellML%userNumber,"*",err,error))
          IF(ASSOCIATED(cellMLModelMaps%cellMLFieldMaps%cellML%region)) &
            & localError=localError//" of region number "// &
            & TRIM(NumberToVString(cellMLModelMaps%cellMLFieldMaps%cellML%region%userNumber,"*",err,error))
        ENDIF
      ENDIF
      localError=localError//" is not allocated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(toMapIndex<1.OR.toMapIndex>cellMLModelMaps%numberOfFieldsMappedTo) THEN
      localError="The specified to map index of "//TRIM(NumberToVString(toMapIndex,"*",err,error))// &
        & " is invalid for the model maps of the CellML field maps"
      IF(ASSOCIATED(cellMLModelMaps%cellMLFieldMaps)) THEN
        IF(ASSOCIATED(cellMLModelMaps%cellMLFieldMaps%cellML)) THEN
          localError=localError//" for CellML user number "// &
            & TRIM(NumberToVString(cellMLModelMaps%cellMLFieldMaps%cellML%userNumber,"*",err,error))
          IF(ASSOCIATED(cellMLModelMaps%cellMLFieldMaps%cellML%region)) &
            & localError=localError//" of region number "// &
            & TRIM(NumberToVString(cellMLModelMaps%cellMLFieldMaps%cellML%region%userNumber,"*",err,error))
        ENDIF
      ENDIF
      localError=localError//". The to map index should be >= 1 and <= "// &
        & TRIM(NumberToVString(cellMLModelMaps%numberOfFieldsMappedTo,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    cellMLModelToMap=>cellMLModelMaps%fieldsMappedTo(toMapIndex)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(cellMLModelToMap)) THEN
      localError="The to model map for to map index "//TRIM(NumberToVString(toMapIndex,"*",err,error))// &
        & " of the model maps of the CellML field maps"
      IF(ASSOCIATED(cellMLModelMaps%cellMLFieldMaps)) THEN
        IF(ASSOCIATED(cellMLModelMaps%cellMLFieldMaps%cellML)) THEN
          localError=localError//" for CellML user number "// &
            & TRIM(NumberToVString(cellMLModelMaps%cellMLFieldMaps%cellML%userNumber,"*",err,error))
          IF(ASSOCIATED(cellMLModelMaps%cellMLFieldMaps%cellML%region)) &
            & localError=localError//" of region number "// &
            & TRIM(NumberToVString(cellMLModelMaps%cellMLFieldMaps%cellML%region%userNumber,"*",err,error))
        ENDIF
      ENDIF
      localError=localError//" is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    EXITS("CellMLModelMaps_CellMLModelToMapGet")
    RETURN
999 NULLIFY(cellMLModelToMap)
998 ERRORSEXITS("CellMLModelMaps_CellMLModelToMapGet",err,error)
    RETURN 1
    
  END SUBROUTINE CellMLModelMaps_CellMLModelToMapGet

  !
  !=================================================================================================================================
  !

  !>Assert that a CellML models field has been finished
  SUBROUTINE CellMLModelsField_AssertIsFinished(cellMLModelsField,err,error,*)

    !Argument Variables
    TYPE(CellMLModelsFieldType), POINTER, INTENT(IN) :: cellMLModelsField !<The CellML models field to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("CellMLModelsField_AssertIsFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(cellMLModelsField)) CALL FlagError("CellML models field is not associated.",err,error,*999)
#endif    

    IF(.NOT.cellMLModelsField%modelsFieldFinished) THEN
      localError="The CellML models field"
      IF(ASSOCIATED(cellMLModelsField%cellML)) THEN
        localError=localError//" for CellML user number "//TRIM(NumberToVString(cellMLModelsField%cellML%userNumber,"*",err,error))
        IF(ASSOCIATED(cellMLModelsField%cellML%region)) &
          & localError=localError//" of region number "// &
          & TRIM(NumberToVString(cellMLModelsField%cellML%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//" has not been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("CellMLModelsField_AssertIsFinished")
    RETURN
999 ERRORSEXITS("CellMLModelsField_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE CellMLModelsField_AssertIsFinished

  !
  !=================================================================================================================================
  !

  !>Assert that a CellML models field has not been finished
  SUBROUTINE CellMLModelsField_AssertNotFinished(cellMLModelsField,err,error,*)

    !Argument Variables
    TYPE(CellMLModelsFieldType), POINTER, INTENT(IN) :: cellMLModelsField !<The CellML models field to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("CellMLModelsField_AssertNotFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(cellMLModelsField)) CALL FlagError("CellML models field is not associated.",err,error,*999)
#endif    

    IF(cellMLModelsField%modelsFieldFinished) THEN
      localError="The CellML models field"
      IF(ASSOCIATED(cellMLModelsField%cellML)) THEN
        localError=localError//" for CellML user number "//TRIM(NumberToVString(cellMLModelsField%cellML%userNumber,"*",err,error))
        IF(ASSOCIATED(cellMLModelsField%cellML%region)) &
          & localError=localError//" of region number "// &
          & TRIM(NumberToVString(cellMLModelsField%cellML%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//" has already been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("CellMLModelsField_AssertNotFinished")
    RETURN
999 ERRORSEXITS("CellMLModelsField_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE CellMLModelsField_AssertNotFinished
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the CellML environment for the specified CellML models field information.
  SUBROUTINE CellMLModelsField_CellMLGet(cellMLModelsField,cellML,err,error,*)

    !Argument variables
    TYPE(CellMLModelsFieldType), POINTER :: cellMLModelsField  !<A pointer to CellML models field information to get the CellML environment from
    TYPE(CellMLType), POINTER :: cellML !<On return, the CellML environment. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("CellMLModelsField_CellMLGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(cellML)) CALL FlagError("CellML environment is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cellMLModelsField)) CALL FlagError("CellML models field is not associated.",err,error,*999)
#endif    
    
    cellML=>cellMLModelsField%cellML

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(cellML)) &
      & CALL FlagError("The CellML environment for the models field information is not associated.",err,error,*999)
#endif    

    EXITS("CellMLModelsField_CellMLGet")
    RETURN
999 NULLIFY(cellML)
998 ERRORSEXITS("CellMLModelsField_CellMLGet",err,error)
    RETURN 1
    
  END SUBROUTINE CellMLModelsField_CellMLGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the models field for the specified CellML models field information.
  SUBROUTINE CellMLModelsField_ModelsFieldGet(cellMLModelsField,modelsField,err,error,*)

    !Argument variables
    TYPE(CellMLModelsFieldType), POINTER :: cellMLModelsField  !<A pointer to CellML models field information to get the models field from
    TYPE(FieldType), POINTER :: modelsField !<On return, the models field. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("CellMLModelsField_ModelsFieldGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(modelsField)) CALL FlagError("Models field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cellMLModelsField)) CALL FlagError("CellML models field is not associated.",err,error,*999)
#endif    
    
    modelsField=>cellMLModelsField%modelsField

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(modelsField)) THEN
      localError="The models field for the models field information"
      IF(ASSOCIATED(cellMLModelsField%cellML)) THEN
        localError=localError//" for CellML user number "//TRIM(NumberToVString(cellMLModelsField%cellML%userNumber,"*",err,error))
        IF(ASSOCIATED(cellMLModelsField%cellML%region)) &
          & localError=localError//" of region number "// &
          & TRIM(NumberToVString(cellMLModelsField%cellML%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//" is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    EXITS("CellMLModelsField_ModelsFieldGet")
    RETURN
999 NULLIFY(modelsField)
998 ERRORSEXITS("CellMLModelsField_ModelsFieldGet",err,error)
    RETURN 1
    
  END SUBROUTINE CellMLModelsField_ModelsFieldGet

  !
  !================================================================================================================================
  !

  !>Returns the only one model index for the specified CellML models field information.
  SUBROUTINE CellMLModelsField_OnlyOneModelIndexGet(cellMLModelsField,onlyOneModelIndex,err,error,*)

    !Argument variables
    TYPE(CellMLModelsFieldType), POINTER :: cellMLModelsField  !<A pointer to CellML models field information to get the only one model index from
    INTEGER(INTG), INTENT(OUT) :: onlyOneModelIndex !<On return, the only one model index in the CellML models field information.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("CellMLModelsField_OnlyOneModelIndexGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(cellMLModelsField)) CALL FlagError("CellML models field is not associated.",err,error,*999)
#endif    
    
    onlyOneModelIndex=cellMLModelsField%onlyOneModelIndex

    EXITS("CellMLModelsField_OnlyOneModelIndexGet")
    RETURN
999 ERRORSEXITS("CellMLModelsField_OnlyOneModelIndexGet",err,error)
    RETURN 1
    
  END SUBROUTINE CellMLModelsField_OnlyOneModelIndexGet

  !
  !=================================================================================================================================
  !

  !>Assert that a CellML parameters field has been finished
  SUBROUTINE CellMLParametersField_AssertIsFinished(cellMLParametersField,err,error,*)

    !Argument Variables
    TYPE(CellMLParametersFieldType), POINTER, INTENT(IN) :: cellMLParametersField !<The CellML parameters field to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("CellMLParametersField_AssertIsFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(cellMLParametersField)) CALL FlagError("CellML parameters field is not associated.",err,error,*999)
#endif    

    IF(.NOT.cellMLParametersField%parametersFieldFinished) THEN
      localError="The CellML parameters field"
      IF(ASSOCIATED(cellMLParametersField%cellML)) THEN
        localError=localError//" for CellML user number "// &
          & TRIM(NumberToVString(cellMLParametersField%cellML%userNumber,"*",err,error))
        IF(ASSOCIATED(cellMLParametersField%cellML%region)) &
          & localError=localError//" of region number "// &
          & TRIM(NumberToVString(cellMLParametersField%cellML%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//" has not been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("CellMLParametersField_AssertIsFinished")
    RETURN
999 ERRORSEXITS("CellMLParametersField_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE CellMLParametersField_AssertIsFinished

  !
  !=================================================================================================================================
  !

  !>Assert that a CellML parameters field has not been finished
  SUBROUTINE CellMLParametersField_AssertNotFinished(cellMLParametersField,err,error,*)

    !Argument Variables
    TYPE(CellMLParametersFieldType), POINTER, INTENT(IN) :: cellMLParametersField !<The CellML parameters field to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("CellMLParametersField_AssertNotFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(cellMLParametersField)) CALL FlagError("CellML parameters field is not associated.",err,error,*999)
#endif    

    IF(cellMLParametersField%parametersFieldFinished) THEN
      localError="The CellML parameters field"
      IF(ASSOCIATED(cellMLParametersField%cellML)) THEN
        localError=localError//" for CellML user number "// &
          & TRIM(NumberToVString(cellMLParametersField%cellML%userNumber,"*",err,error))
        IF(ASSOCIATED(cellMLParametersField%cellML%region)) &
          & localError=localError//" of region number "// &
          & TRIM(NumberToVString(cellMLParametersField%cellML%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//" has already been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("CellMLParametersField_AssertNotFinished")
    RETURN
999 ERRORSEXITS("CellMLParametersField_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE CellMLParametersField_AssertNotFinished
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the CellML environment for the specified CellML parameters field information.
  SUBROUTINE CellMLParametersField_CellMLGet(cellMLParametersField,cellML,err,error,*)

    !Argument variables
    TYPE(CellMLParametersFieldType), POINTER :: cellMLParametersField  !<A pointer to CellML parameters field information to get the CellML environment from
    TYPE(CellMLType), POINTER :: cellML !<On return, the CellML environment. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("CellMLParametersField_CellMLGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(cellML)) CALL FlagError("CellML environment is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cellMLParametersField)) CALL FlagError("CellML parameters field is not associated.",err,error,*999)
#endif    
    
    cellML=>cellMLParametersField%cellML

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(cellML)) &
      & CALL FlagError("The CellML environment for the parameters field information is not associated.",err,error,*999)
#endif    

    EXITS("CellMLParametersField_CellMLGet")
    RETURN
999 NULLIFY(cellML)
998 ERRORSEXITS("CellMLParametersField_CellMLGet",err,error)
    RETURN 1
    
  END SUBROUTINE CellMLParametersField_CellMLGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the parameters field for the specified CellML parameters field information.
  SUBROUTINE CellMLParametersField_ParametersFieldGet(cellMLParametersField,parametersField,err,error,*)

    !Argument variables
    TYPE(CellMLParametersFieldType), POINTER :: cellMLParametersField  !<A pointer to CellML parameters field information to get the parameters field from
    TYPE(FieldType), POINTER :: parametersField !<On return, the parameters field. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("CellMLParametersField_ParametersFieldGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(parametersField)) CALL FlagError("Parameters field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cellMLParametersField)) CALL FlagError("CellML parameters field is not associated.",err,error,*999)
#endif    
    
    parametersField=>cellMLParametersField%parametersField

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(parametersField)) THEN
      localError="The parameters field for the parameters field information"
      IF(ASSOCIATED(cellMLParametersField%cellML)) THEN
        localError=localError//" for CellML user number "// &
          & TRIM(NumberToVString(cellMLParametersField%cellML%userNumber,"*",err,error))
        IF(ASSOCIATED(cellMLParametersField%cellML%region)) &
          & localError=localError//" of region number "// &
          & TRIM(NumberToVString(cellMLParametersField%cellML%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//" is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    EXITS("CellMLParametersField_ParametersFieldGet")
    RETURN
999 NULLIFY(parametersField)
998 ERRORSEXITS("CellMLParametersField_ParametersFieldGet",err,error)
    RETURN 1
    
  END SUBROUTINE CellMLParametersField_ParametersFieldGet

  !
  !=================================================================================================================================
  !

  !>Assert that a CellML state field has been finished
  SUBROUTINE CellMLStateField_AssertIsFinished(cellMLStateField,err,error,*)

    !Argument Variables
    TYPE(CellMLStateFieldType), POINTER, INTENT(IN) :: cellMLStateField !<The CellML state field to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("CellMLStateField_AssertIsFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(cellMLStateField)) CALL FlagError("CellML state field is not associated.",err,error,*999)
#endif    

    IF(.NOT.cellMLStateField%stateFieldFinished) THEN
      localError="The CellML state field"
      IF(ASSOCIATED(cellMLStateField%cellML)) THEN
        localError=localError//" for CellML user number "//TRIM(NumberToVString(cellMLStateField%cellML%userNumber,"*",err,error))
        IF(ASSOCIATED(cellMLStateField%cellML%region)) &
          & localError=localError//" of region number "// &
          & TRIM(NumberToVString(cellMLStateField%cellML%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//" has not been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("CellMLStateField_AssertIsFinished")
    RETURN
999 ERRORSEXITS("CellMLStateField_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE CellMLStateField_AssertIsFinished

  !
  !=================================================================================================================================
  !

  !>Assert that a CellML state field has not been finished
  SUBROUTINE CellMLStateField_AssertNotFinished(cellMLStateField,err,error,*)

    !Argument Variables
    TYPE(CellMLStateFieldType), POINTER, INTENT(IN) :: cellMLStateField !<The CellML state field to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("CellMLStateField_AssertNotFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(cellMLStateField)) CALL FlagError("CellML state field is not associated.",err,error,*999)
#endif    

    IF(cellMLStateField%stateFieldFinished) THEN
      localError="The CellML state field"
      IF(ASSOCIATED(cellMLStateField%cellML)) THEN
        localError=localError//" for CellML user number "//TRIM(NumberToVString(cellMLStateField%cellML%userNumber,"*",err,error))
        IF(ASSOCIATED(cellMLStateField%cellML%region)) &
          & localError=localError//" of region number "// &
          & TRIM(NumberToVString(cellMLStateField%cellML%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//" has already been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("CellMLStateField_AssertNotFinished")
    RETURN
999 ERRORSEXITS("CellMLStateField_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE CellMLStateField_AssertNotFinished
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the CellML environment for the specified CellML state field information.
  SUBROUTINE CellMLStateField_CellMLGet(cellMLStateField,cellML,err,error,*)

    !Argument variables
    TYPE(CellMLStateFieldType), POINTER :: cellMLStateField  !<A pointer to CellML state field information to get the CellML environment from
    TYPE(CellMLType), POINTER :: cellML !<On return, the CellML environment. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("CellMLStateField_CellMLGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(cellML)) CALL FlagError("CellML environment is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cellMLStateField)) CALL FlagError("CellML state field is not associated.",err,error,*999)
#endif    
    
    cellML=>cellMLStateField%cellML

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(cellML)) &
      & CALL FlagError("The CellML environment for the state field information is not associated.",err,error,*999)
#endif    

    EXITS("CellMLStateField_CellMLGet")
    RETURN
999 NULLIFY(cellML)
998 ERRORSEXITS("CellMLStateField_CellMLGet",err,error)
    RETURN 1
    
  END SUBROUTINE CellMLStateField_CellMLGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the state field for the specified CellML state field information.
  SUBROUTINE CellMLStateField_StateFieldGet(cellMLStateField,stateField,err,error,*)

    !Argument variables
    TYPE(CellMLStateFieldType), POINTER :: cellMLStateField  !<A pointer to CellML state field information to get the state field from
    TYPE(FieldType), POINTER :: stateField !<On return, the state field. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("CellMLStateField_StateFieldGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(stateField)) CALL FlagError("State field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cellMLStateField)) CALL FlagError("CellML state field is not associated.",err,error,*999)
#endif    
    
    stateField=>cellMLStateField%stateField

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(stateField)) THEN
      localError="The state field for the state field information"
      IF(ASSOCIATED(cellMLStateField%cellML)) THEN
        localError=localError//" for CellML user number "//TRIM(NumberToVString(cellMLStateField%cellML%userNumber,"*",err,error))
        IF(ASSOCIATED(cellMLStateField%cellML%region)) &
          & localError=localError//" of region number "// &
          & TRIM(NumberToVString(cellMLStateField%cellML%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//" is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    EXITS("CellMLStateField_StateFieldGet")
    RETURN
999 NULLIFY(stateField)
998 ERRORSEXITS("CellMLStateField_StateFieldGet",err,error)
    RETURN 1
    
  END SUBROUTINE CellMLStateField_StateFieldGet

  !
  !================================================================================================================================
  !

END MODULE CellMLAccessRoutines
