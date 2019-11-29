!> \file
!> \author Chris Bradley
!> \brief This module contains all region access method routines.
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

!> This module contains all region access method routines.
MODULE RegionAccessRoutines
  
  USE BaseRoutines
  USE CellMLAccessRoutines
  USE DataPointAccessRoutines
  USE DecompositionAccessRoutines
  USE EquationsSetAccessRoutines
  USE FieldAccessRoutines
  USE GeneratedMeshAccessRoutines
  USE InterfaceAccessRoutines
  USE ISO_VARYING_STRING
  USE Kinds
  USE MeshAccessRoutines
  USE Strings
  USE Types
  
#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  INTERFACE REGION_COORDINATE_SYSTEM_GET
    MODULE PROCEDURE Region_CoordinateSystemGet
  END INTERFACE REGION_COORDINATE_SYSTEM_GET

  INTERFACE Region_LabelGet
    MODULE PROCEDURE Region_LabelGetC
    MODULE PROCEDURE Region_LabelGetVS
  END INTERFACE Region_LabelGet
  
  INTERFACE REGION_NODES_GET
    MODULE PROCEDURE Region_NodesGet
  END INTERFACE REGION_NODES_GET

  INTERFACE REGION_USER_NUMBER_FIND
    MODULE PROCEDURE Region_UserNumberFind
  END INTERFACE REGION_USER_NUMBER_FIND

  PUBLIC Region_AssertIsFinished,Region_AssertNotFinished

  PUBLIC Region_CellMLGet

  PUBLIC Region_ContextGet

  PUBLIC Region_CoordinateSystemGet

  PUBLIC REGION_COORDINATE_SYSTEM_GET
  
  PUBLIC Region_DataPointsGet

  PUBLIC Region_DecomposerGet
  
  PUBLIC Region_DecomposersGet
  
  PUBLIC Region_EquationsSetGet

  PUBLIC Region_FieldGet

  PUBLIC Region_FieldsGet

  PUBLIC Region_GeneratedMeshGet

  PUBLIC Region_Get

  PUBLIC Region_InterfaceGet

  PUBLIC Region_IsSubRegion

  PUBLIC Region_LabelGet

  PUBLIC Region_MeshGet

  PUBLIC Region_MeshesGet

  PUBLIC Region_NodesGet

  PUBLIC REGION_NODES_GET

  PUBLIC Region_RegionsGet

  PUBLIC Region_UserNumberFind

  PUBLIC REGION_USER_NUMBER_FIND

  PUBLIC Region_UserNumberGet

  PUBLIC Regions_WorldRegionGet

CONTAINS

  !
  !=================================================================================================================================
  !

  !>Assert that a region has been finished
  SUBROUTINE Region_AssertIsFinished(region,err,error,*)

    !Argument Variables
    TYPE(RegionType), POINTER, INTENT(INOUT) :: region !<The region to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Region_AssertIsFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated.",err,error,*999)

    IF(.NOT.region%regionFinished) THEN
      localError="Region number "//TRIM(NumberToVString(region%userNumber,"*",err,error))// &
        & " has not been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Region_AssertIsFinished")
    RETURN
999 ERRORSEXITS("Region_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE Region_AssertIsFinished

  !
  !=================================================================================================================================
  !

  !>Assert that a region has not been finished
  SUBROUTINE Region_AssertNotFinished(region,err,error,*)

    !Argument Variables
    TYPE(RegionType), POINTER, INTENT(INOUT) :: region !<The region to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Region_AssertNotFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated.",err,error,*999)

    IF(region%regionFinished) THEN
      localError="Region number "//TRIM(NumberToVString(region%userNumber,"*",err,error))// &
        & " has already been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Region_AssertNotFinished")
    RETURN
999 ERRORSEXITS("Region_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE Region_AssertNotFinished

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the cellml for a given user number in a region. \see OPENCMISS::Iron::cmfe_Region_CellMLGet
  SUBROUTINE Region_CellMLGet(region,userNumber,cellml,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: region !<A pointer to the region to get the cellml for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the cellml to get.
    TYPE(CELLML_TYPE), POINTER :: cellml !<On exit, a pointer to the cellml for the region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Region_CellMLGet",err,error,*998)

    IF(ASSOCIATED(cellml)) CALL FlagError("CellML is already associated.",err,error,*998)
    CALL Region_AssertIsFinished(region,err,error,*998)
     
    NULLIFY(cellml)
    CALL CellML_UserNumberFind(userNumber,region,cellml,err,error,*999)
    IF(.NOT.ASSOCIATED(cellml)) THEN
      localError="A cellml environment with a user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " do not exist on region number "//TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Region_CellMLGet")
    RETURN
999 NULLIFY(cellml)
998 ERRORSEXITS("Region_CellMLGet",err,error)
    RETURN 1
    
  END SUBROUTINE Region_CellMLGet
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the context for a region.
  SUBROUTINE Region_ContextGet(region,context,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: region !<A pointer to the region to get the context for
    TYPE(ContextType), POINTER :: context !<On exit, a pointer to the context for the region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Region_ContextGet",err,error,*998)

    IF(ASSOCIATED(context)) CALL FlagError("Context is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated.",err,error,*999)

    IF(.NOT.ASSOCIATED(region%regions)) THEN
      localError="Regions is not associated for region number "// &
        & TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    context=>region%regions%context
    IF(.NOT.ASSOCIATED(context)) THEN
      localError="The context is not associated for the regions for region number "// &
        & TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Region_ContextGet")
    RETURN
999 NULLIFY(context)
998 ERRORSEXITS("Region_ContextGet",err,error)
    RETURN 1
    
  END SUBROUTINE Region_ContextGet

  !
  !================================================================================================================================
  !

  !>Returns the coordinate system of region. \see OPENCMISS::Iron::cmfe_RegionCoordinateSystemGet
  SUBROUTINE Region_CoordinateSystemGet(region,coordinateSystem,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: region !<A pointer to the region to get the coordinate system for
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem !<On exit, the coordinate system for the specified region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Region_CoordinateSystemGet",err,error,*998)

    !Check input arguments
    IF(ASSOCIATED(coordinateSystem)) CALL FlagError("Coordinate system is already associated.",err,error,*998)
    CALL Region_AssertIsFinished(region,err,error,*998)

    coordinateSystem=>region%coordinateSystem
    IF(.NOT.ASSOCIATED(coordinateSystem)) THEN
      localError="The coordinate system for region number "//TRIM(NumberToVString(region%userNumber,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Region_CoordinateSystemGet")
    RETURN
999 NULLIFY(coordinateSystem)
998 ERRORSEXITS("Region_CoordinateSystemGet",err,error)    
    RETURN 1
    
  END SUBROUTINE Region_CoordinateSystemGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the data points for a given user number in a region. \see OPENCMISS::Iron::cmfe_Region_DataPointsGet
  SUBROUTINE Region_DataPointsGet(region,userNumber,dataPoints,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: region !<A pointer to the region to get the data points for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the data points to get.
    TYPE(DataPointsType), POINTER :: dataPoints !<On exit, a pointer to the data points for the region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Region_DataPointsGet",err,error,*998)

    IF(ASSOCIATED(dataPoints)) CALL FlagError("Data points is already associated.",err,error,*998)
    CALL Region_AssertIsFinished(region,err,error,*998)
    IF(.NOT.ASSOCIATED(region%dataPointSets)) THEN
      localError="Region data point sets is not associated for region number "// &
        & TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    NULLIFY(dataPoints)
    CALL DataPointSets_UserNumberFind(region%dataPointSets,userNumber,dataPoints,err,error,*999)
    IF(.NOT.ASSOCIATED(dataPoints)) THEN
      localError="Data points with a user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " does not exist on region number "//TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("Region_DataPointsGet")
    RETURN
999 NULLIFY(dataPoints)
998 ERRORSEXITS("Region_DataPointsGet",err,error)
    RETURN 1
    
  END SUBROUTINE Region_DataPointsGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the decomposer for a given user number in a region. \see OPENCMISS::Iron::cmfe_Region_DecomposerGet
  SUBROUTINE Region_DecomposerGet(region,userNumber,decomposer,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: region !<A pointer to the region to get the decomposer for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the decomposer to get.
    TYPE(DecomposerType), POINTER :: decomposer !<On exit, a pointer to the decomposer for the region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Region_DecomposerGet",err,error,*998)

    IF(ASSOCIATED(decomposer)) CALL FlagError("Decomposer is already associated.",err,error,*998)
    CALL Region_AssertIsFinished(region,err,error,*999)
    
    NULLIFY(decomposer)
    CALL Decomposer_UserNumberFind(userNumber,region,decomposer,err,error,*999)
    IF(.NOT.ASSOCIATED(decomposer)) THEN
      localError="A decomposer with a user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " does not exist on region number "//TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Region_DecomposerGet")
    RETURN
999 NULLIFY(decomposer)
998 ERRORSEXITS("Region_DecomposerGet",err,error)
    RETURN 1
    
  END SUBROUTINE Region_DecomposerGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the decomposers for a region. 
  SUBROUTINE Region_DecomposersGet(region,decomposers,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: region !<A pointer to the region to get the decomposers for
    TYPE(DecomposersType), POINTER :: decomposers !<On exit, a pointer to the decomposers for the region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Region_DecomposersGet",err,error,*998)

    IF(ASSOCIATED(decomposers)) CALL FlagError("Decomposers is already associated.",err,error,*998)
    CALL Region_AssertIsFinished(region,err,error,*998)

    decomposers=>region%decomposers
    IF(.NOT.ASSOCIATED(decomposers)) THEN
      localError="The decomposers for region number "//TRIM(NumberToVString(region%userNumber,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("Region_DecomposersGet")
    RETURN
999 NULLIFY(decomposers)
998 ERRORSEXITS("Region_DecomposersGet",err,error)
    RETURN 1
    
  END SUBROUTINE Region_DecomposersGet
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the equations set for a given user number in a region. \see OPENCMISS::Iron::cmfe_Region_EquationsSetGet
  SUBROUTINE Region_EquationsSetGet(region,userNumber,equationsSet,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: region !<A pointer to the region to get the equationsSet for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the equations set to get.
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<On exit, a pointer to the equations set for the region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Region_EquationsSetGet",err,error,*998)

    IF(ASSOCIATED(equationsSet)) CALL FlagError("Equations set is already associated.",err,error,*998)
    CALL Region_AssertIsFinished(region,err,error,*998)
    
    NULLIFY(equationsSet)
    CALL EquationsSet_UserNumberFind(userNumber,region,equationsSet,err,error,*999)
    IF(.NOT.ASSOCIATED(equationsSet)) THEN
      localError="An equations set with a user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " does not exist on region number "//TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Region_EquationsSetGet")
    RETURN
999 NULLIFY(equationsSet)
998 ERRORSEXITS("Region_EquationsSetGet",err,error)
    RETURN 1
    
  END SUBROUTINE Region_EquationsSetGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the field for a given user number in a region. \see OPENCMISS::Iron::cmfe_Region_FieldGet
  SUBROUTINE Region_FieldGet(region,userNumber,field,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: region !<A pointer to the region to get the field for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the field to get.
    TYPE(FieldType), POINTER :: field !<On exit, a pointer to the field for the region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Region_FieldGet",err,error,*998)

    IF(ASSOCIATED(field)) CALL FlagError("Field is already associated.",err,error,*998)
    CALL Region_AssertIsFinished(region,err,error,*998)
    
    NULLIFY(field)
    CALL Field_UserNumberFind(userNumber,region,field,err,error,*999)
    IF(.NOT.ASSOCIATED(field)) THEN
      localError="A field with a user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " does not exist on region number "//TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Region_FieldGet")
    RETURN
999 NULLIFY(field)
998 ERRORSEXITS("Region_FieldGet",err,error)
    RETURN 1
    
  END SUBROUTINE Region_FieldGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the fields for a given region.
  SUBROUTINE Region_FieldsGet(region,fields,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: region !<A pointer to the region to get the fields for
    TYPE(FieldsType), POINTER :: fields !<On exit, a pointer to the fields for the region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Region_FieldsGet",err,error,*998)

    IF(ASSOCIATED(fields)) CALL FlagError("Fields is already associated.",err,error,*998)
    CALL Region_AssertIsFinished(region,err,error,*998)

    fields=>region%fields
    IF(.NOT.ASSOCIATED(fields)) THEN
      localError="The fields for region number "//TRIM(NumberToVString(region%userNumber,"*",err,error))// &
        & " are not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Region_FieldsGet")
    RETURN
999 NULLIFY(fields)
998 ERRORSEXITS("Region_FieldsGet",err,error)
    RETURN 1
    
  END SUBROUTINE Region_FieldsGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the generated mesh for a given user number in a region. \see OPENCMISS::Iron::cmfe_Region_GeneratedMeshGet
  SUBROUTINE Region_GeneratedMeshGet(region,userNumber,generatedMesh,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: region !<A pointer to the region to get the generated mesh for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the generated mesh to get.
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<On exit, a pointer to the generated mesh for the region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Region_GeneratedMeshGet",err,error,*998)

    IF(ASSOCIATED(generatedMesh)) CALL FlagError("Generated mesh is already associated.",err,error,*998)
    CALL Region_AssertIsFinished(region,err,error,*998)
    
    NULLIFY(generatedMesh)
    CALL GeneratedMesh_UserNumberFind(userNumber,region,generatedMesh,err,error,*999)
    IF(.NOT.ASSOCIATED(generatedMesh)) THEN
      localError="A generated mesh with a user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " does not exist on region number "//TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Region_GeneratedMeshGet")
    RETURN
999 NULLIFY(generatedMesh)
998 ERRORSEXITS("Region_GeneratedMeshGet",err,error)
    RETURN 1
    
  END SUBROUTINE Region_GeneratedMeshGet

   !
  !================================================================================================================================
  !

  !>Finds and returns a pointer to the region with the given user number. 
  SUBROUTINE Region_Get(regions,userNumber,region,err,error,*)

    !Argument variables
    TYPE(RegionsType), POINTER :: regions !<The regions to get the user number for.
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the region to find
    TYPE(RegionType), POINTER :: region !<On exit, a pointer to the region with the specified user number if it exists. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Region_Get",err,error,*999)

    CALL Region_UserNumberFind(regions,userNumber,region,err,error,*999)
    IF(.NOT.ASSOCIATED(region)) THEN
      localError="A region with an user number of "//TRIM(NumberToVString(userNumber,"*",err,error))//" does not exist."
      CALL FlagError(localError,err,error,*999)
    ENDIF
  
    EXITS("Region_Get")
    RETURN
999 ERRORSEXITS("Region_Get",err,error)
    RETURN 1
    
  END SUBROUTINE Region_Get

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the interface for a given user number in a (parent) region. \see OPENCMISS::Iron::cmfe_Region_InterfaceGet
  SUBROUTINE Region_InterfaceGet(region,userNumber,interface,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: region !<A pointer to the parent region to get the interface for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the interface to get.
    TYPE(InterfaceType), POINTER :: interface !<On exit, a pointer to the interface for the parent region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Region_InterfaceGet",err,error,*998)

    IF(ASSOCIATED(INTERFACE)) CALL FlagError("Interface is already associated.",err,error,*998)
    CALL Region_AssertIsFinished(region,err,error,*998)
     
    NULLIFY(interface)
    CALL Interface_UserNumberFind(userNumber,region,interface,err,error,*999)
    IF(.NOT.ASSOCIATED(interface)) THEN
      localError="An interface with a user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " does not exist on region number "//TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Region_InterfaceGet")
    RETURN
999 NULLIFY(interface)
998 ERRORSEXITS("Region_InterfaceGet",err,error)
    RETURN 1
    
  END SUBROUTINE Region_InterfaceGet

  !
  !================================================================================================================================
  !

  !>Returns if a region is a sub-region (or the same region) of another region.
  SUBROUTINE Region_IsSubRegion(region,subRegion,isSubRegion,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: region !<A pointer to the region to check the sub-region for
    TYPE(RegionType), POINTER :: subRegion !<A pointer to the sub-region to region for
    LOGICAL, INTENT(OUT) :: isSubRegion !<On return, is .TRUE. if sub-region is either the same as region or a sub-region of region, .FALSE. if not.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(RegionType), POINTER :: parentRegion
 
    ENTERS("Region_IsSubRegion",err,error,*999)

    isSubRegion=.FALSE.
    CALL Region_AssertIsFinished(region,err,error,*999)
    CALL Region_AssertIsFinished(subRegion,err,error,*999)
    
    parentRegion=>subRegion
    DO WHILE(ASSOCIATED(parentRegion))
      IF(ASSOCIATED(region,parentRegion)) THEN
        isSubRegion=.TRUE.
        EXIT
      ELSE
        parentRegion=>parentRegion%parentRegion
      ENDIF
    ENDDO
    
    EXITS("Region_IsSubRegion")
    RETURN
999 ERRORSEXITS("Region_IsSubRegion",err,error)
    RETURN 1
    
  END SUBROUTINE Region_IsSubRegion

  !
  !================================================================================================================================
  !

  !>Returns the label of a region. \see OPENCMISS::Iron::cmfe_Region_LabelGet
  SUBROUTINE Region_LabelGetC(region,label,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: region !<A pointer to the region to get the label for
    CHARACTER(LEN=*), INTENT(OUT) :: label !<On return the region label.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: cLength,vsLength

    ENTERS("Region_LabelGetC",err,error,*999)

    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated.",err,error,*999)
    
    cLength=LEN(label)
    vsLength=LEN_TRIM(region%label)
    IF(cLength>vsLength) THEN
      label=CHAR(region%label,vsLength)
    ELSE
      label=CHAR(region%label,cLength)
    ENDIF
    
    EXITS("Region_LabelGetC")
    RETURN
999 ERRORSEXITS("Region_LabelGetC",err,error)
    RETURN 1
    
  END SUBROUTINE Region_LabelGetC

   !
  !================================================================================================================================
  !

  !>Returns the label of a region. \see OPENCMISS::Iron::cmfe_Region_LabelGet
  SUBROUTINE Region_LabelGetVS(region,label,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: region !<A pointer to the region to get the label for
    TYPE(VARYING_STRING), INTENT(OUT) :: label !<On return the region label.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Region_LabelGetVS",err,error,*999)

    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated.",err,error,*999)
    
    !CPB 20/2/07 The following line crashes the AIX compiler unless it has a VAR_STR(CHAR()) around it
    label=VAR_STR(CHAR(region%label))
     
    EXITS("Region_LabelGetVS")
    RETURN
999 ERRORSEXITS("Region_LabelGetVS",err,error)
    RETURN 1
    
  END SUBROUTINE Region_LabelGetVS

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the mesh for a given user number in a region. \see OPENCMISS::Iron::cmfe_Region_MeshGet
  SUBROUTINE Region_MeshGet(region,userNumber,mesh,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: region !<A pointer to the region to get the mesh for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the mesh to get.
    TYPE(MeshType), POINTER :: mesh !<On exit, a pointer to the mesh for the region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Region_MeshGet",err,error,*998)

    IF(ASSOCIATED(mesh)) CALL FlagError("Mesh is already associated.",err,error,*998)
    CALL Region_AssertIsFinished(region,err,error,*998)
    
    NULLIFY(mesh)
    CALL Mesh_UserNumberFind(userNumber,region,mesh,err,error,*999)
    IF(.NOT.ASSOCIATED(mesh)) THEN
      localError="A mesh with a user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " does not exist on region number "//TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Region_MeshGet")
    RETURN
999 NULLIFY(mesh)
998 ERRORSEXITS("Region_MeshGet",err,error)
    RETURN 1
    
  END SUBROUTINE Region_MeshGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the meshes for a region.
  SUBROUTINE Region_MeshesGet(region,meshes,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: region !<A pointer to the region to get the meshes for
    TYPE(MeshesType), POINTER :: meshes !<On exit, a pointer to the meshes for the region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Region_MeshesGet",err,error,*998)

    IF(ASSOCIATED(meshes)) CALL FlagError("Meshes is already associated.",err,error,*998)
    CALL Region_AssertIsFinished(region,err,error,*999)

    meshes=>region%meshes
    IF(.NOT.ASSOCIATED(meshes)) THEN
      localError="The meshes for region number "//TRIM(NumberToVString(region%userNumber,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("Region_MeshesGet")
    RETURN
999 NULLIFY(meshes)
998 ERRORSEXITS("Region_MeshesGet",err,error)
    RETURN 1
    
  END SUBROUTINE Region_MeshesGet
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the nodes for a region. \see OPENCMISS::Iron::cmfe_Region_NodesGet
  SUBROUTINE Region_NodesGet(region,nodes,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: region !<A pointer to the region to get the nodes for
    TYPE(NodesType), POINTER :: nodes !<On exit, a pointer to the nodes for the region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Region_NodesGet",err,error,*998)

    IF(ASSOCIATED(nodes)) CALL FlagError("Nodes is already associated.",err,error,*998)
    CALL Region_AssertIsFinished(region,err,error,*998)

    nodes=>region%nodes
    IF(.NOT.ASSOCIATED(nodes)) THEN
      localError="The nodes for region number "//TRIM(NumberToVString(region%userNumber,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("Region_NodesGet")
    RETURN
999 NULLIFY(nodes)
998 ERRORSEXITS("Region_NodesGet",err,error)
    RETURN 1
    
  END SUBROUTINE Region_NodesGet
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the regions for a region. 
  SUBROUTINE Region_RegionsGet(region,regions,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: region !<A pointer to the region to get the regions for
    TYPE(RegionsType), POINTER :: regions !<On exit, a pointer to the regions for the region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Region_RegionsGet",err,error,*998)

    IF(ASSOCIATED(regions)) CALL FlagError("Regions is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated.",err,error,*999)
    
    regions=>region%regions   
    IF(.NOT.ASSOCIATED(regions)) THEN
      localError="Regions is not associated for region number "// &
        & TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("Region_RegionsGet")
    RETURN
999 NULLIFY(regions)
998 ERRORSEXITS("Region_RegionsGet",err,error)
    RETURN 1
    
  END SUBROUTINE Region_RegionsGet
  
  !
  !================================================================================================================================
  !

  !>Finds and returns a pointer to the region with the given user number. If no region with that number exists region is left nullified.
  SUBROUTINE Region_UserNumberFind(regions,userNumber,region,err,error,*)

    !Argument variables
    TYPE(RegionsType), POINTER :: regions !<The regions to find the user number for.
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the region to find
    TYPE(RegionType), POINTER :: region !<On exit, a pointer to the region with the specified user number if it exists. If no region exists with the specified user number a NULL pointer is returned. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: regionIdx
    TYPE(RegionType), POINTER :: worldRegion
    
    ENTERS("Region_UserNumberFind",err,error,*999)

    IF(ASSOCIATED(region)) CALL FlagError("Region is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(regions)) CALL FlagError("Regions is not associated.",err,error,*999)
    
    worldRegion=>regions%worldRegion
    IF(.NOT.ASSOCIATED(worldRegion)) CALL FlagError("World region is not associated.",err,error,*999)
    
    NULLIFY(region)
    IF(userNumber==0) THEN
      region=>worldRegion
    ELSE
      IF(ASSOCIATED(worldRegion%subRegions)) THEN
        DO regionIdx=1,worldRegion%numberOfSubRegions
          CALL Region_UserNumberFindPtr(userNumber,worldRegion%subRegions(regionIdx)%ptr,region,err,error,*999)
          IF(ASSOCIATED(region)) EXIT
        ENDDO !regionIdx
      ENDIF
    ENDIF
  
    EXITS("Region_UserNumberFind")
    RETURN
999 ERRORSEXITS("Region_UserNumberFind",err,error)
    RETURN 1
    
  END SUBROUTINE Region_UserNumberFind

  !
  !================================================================================================================================
  !

  !>Finds and returns a pointer to the region with the given user number starting from the given start region and searching all sub-regions under the start region. If no region with that number exists region is left nullified.
  RECURSIVE SUBROUTINE Region_UserNumberFindPtr(userNumber,startRegion,region,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number to find
    TYPE(RegionType), POINTER :: startRegion !<A pointer to the region to start the search from
    TYPE(RegionType), POINTER :: region !<On exit, a pointer to the region with the specified user number if it exists. If no region exists with the specified user number a NULL pointer is returned.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: regionIdx

    ENTERS("Region_UserNumberFindPtr",err,error,*999)

    IF(.NOT.ASSOCIATED(startRegion)) CALL FlagError("Start region is not associated",err,error,*999)
    IF(ASSOCIATED(region)) CALL FlagError("Region is already associated.",err,error,*999)

    NULLIFY(region)
    IF(startRegion%userNumber==userNumber) THEN
      region=>startRegion
    ELSE
      IF(ASSOCIATED(startRegion%subRegions)) THEN
        DO regionIdx=1,startRegion%numberOfSubRegions
          CALL Region_UserNumberFindPtr(userNumber,startRegion%subRegions(regionIdx)%ptr,region,err,error,*999)
          IF(ASSOCIATED(region)) EXIT
        ENDDO !regionIdx
      ENDIF
    ENDIF
    
    EXITS("Region_UserNumberFindPtr")
    RETURN
999 ERRORSEXITS("Region_UserNumberFindPtr",err,error)
    RETURN 1
    
  END SUBROUTINE Region_UserNumberFindPtr

  !
  !================================================================================================================================
  !

  !>Returns the user number for a region.
  SUBROUTINE Region_UserNumberGet(region,userNumber,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: region !<A pointer to the region to get the user number for
    INTEGER(INTG), INTENT(OUT) :: userNumber !<On exit, the user number of the region.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Region_UserNumberGet",err,error,*999)

    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated.",err,error,*999)

    userNumber=region%userNumber
  
    EXITS("Region_UserNumberGet")
    RETURN
999 ERRORSEXITS("Region_UserNumberGet",err,error)
    RETURN 1
    
  END SUBROUTINE Region_UserNumberGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the world region for a regions.
  SUBROUTINE Regions_WorldRegionGet(regions,worldRegion,err,error,*)

    !Argument variables
    TYPE(RegionsType), POINTER :: regions !<A pointer to the regions to get the world region for
    TYPE(RegionType), POINTER :: worldRegion !<On exit, a pointer to the world region for the regions. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Regions_WorldRegionGet",err,error,*998)

    IF(ASSOCIATED(worldRegion)) CALL FlagError("World region is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(regions)) CALL FlagError("Regions is not associated.",err,error,*998)

    worldRegion=>regions%worldREgion
    IF(.NOT.ASSOCIATED(worldRegion)) CALL FlagError("Regions world region is not associated.",err,error,*999)
       
    EXITS("Regions_WorldRegionGet")
    RETURN
999 NULLIFY(worldRegion)
998 ERRORSEXITS("Regions_WorldRegionGet",err,error)
    RETURN 1
    
  END SUBROUTINE Regions_WorldRegionGet
  
  !
  !================================================================================================================================
  !

END MODULE RegionAccessRoutines
