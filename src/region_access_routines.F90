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

!> This module contains all region access method routines.
MODULE RegionAccessRoutines
  
  USE BaseRoutines
  USE CellMLAccessRoutines
  USE DataPointAccessRoutines
  USE EquationsSetAccessRoutines
  USE FieldAccessRoutines
  USE GeneratedMeshAccessRoutines
  USE InterfaceAccessRoutines
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

  TYPE(REGIONS_TYPE) :: regions
  
  !Interfaces

  INTERFACE REGION_COORDINATE_SYSTEM_GET
    MODULE PROCEDURE Region_CoordinateSystemGet
  END INTERFACE REGION_COORDINATE_SYSTEM_GET

  INTERFACE REGION_NODES_GET
    MODULE PROCEDURE Region_NodesGet
  END INTERFACE REGION_NODES_GET

  INTERFACE REGION_USER_NUMBER_FIND
    MODULE PROCEDURE Region_UserNumberFind
  END INTERFACE REGION_USER_NUMBER_FIND

  PUBLIC regions

  PUBLIC Region_CellMLGet

  PUBLIC Region_CoordinateSystemGet

  PUBLIC REGION_COORDINATE_SYSTEM_GET
  
  PUBLIC Region_DataPointsGet

  PUBLIC Region_EquationsSetGet

  PUBLIC Region_FieldGet

  PUBLIC Region_GeneratedMeshGet

  PUBLIC Region_Get

  PUBLIC Region_InterfaceGet

  PUBLIC Region_MeshGet

  PUBLIC Region_NodesGet

  PUBLIC REGION_NODES_GET

  PUBLIC Region_UserNumberFind

  PUBLIC REGION_USER_NUMBER_FIND

  PUBLIC Region_UserNumberGet

CONTAINS

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the cellml for a given user number in a region. \see OPENCMISS::Iron::cmfe_Region_CellMLGet
  SUBROUTINE Region_CellMLGet(region,userNumber,cellml,err,error,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: region !<A pointer to the region to get the cellml for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the cellml to get.
    TYPE(CELLML_TYPE), POINTER :: cellml !<On exit, a pointer to the cellml for the region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Region_CellMLGet",err,error,*998)

    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated.",err,error,*998)
    IF(.NOT.region%REGION_FINISHED) CALL FlagError("Region has not been finished.",err,error,*998)
    IF(ASSOCIATED(cellml)) CALL FlagError("CellML is already associated.",err,error,*998)
    
    NULLIFY(cellml)
    CALL CellML_UserNumberFind(userNumber,region,cellml,err,error,*999)
    IF(.NOT.ASSOCIATED(cellml)) THEN
      localError="A cellml environment with a user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " do not exist on region number "//TRIM(NumberToVString(region%USER_NUMBER,"*",err,error))//"."
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

  !>Returns the coordinate system of region. \see OPENCMISS::Iron::cmfe_RegionCoordinateSystemGet
  SUBROUTINE Region_CoordinateSystemGet(region,coordinateSystem,err,error,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: region !<A pointer to the region to get the coordinate system for
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: coordinateSystem !<On exit, the coordinate system for the specified region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Region_CoordinateSystemGet",ERR,ERROR,*999)

    !Check input arguments
    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated.",err,error,*999)
    IF(.NOT.region%REGION_FINISHED) CALL FlagError("Region has not been finished.",err,error,*999)
    IF(ASSOCIATED(coordinateSystem)) CALL FlagError("Coordinate system is already associated.",ERR,ERROR,*999)

    coordinateSystem=>region%COORDINATE_SYSTEM
    IF(.NOT.ASSOCIATED(coordinateSystem)) THEN
      localError="The coordinate system for region number "//TRIM(NumberToVString(region%USER_NUMBER,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Region_CoordinateSystemGet")
    RETURN
999 ERRORSEXITS("Region_CoordinateSystemGet",err,error)
    RETURN 1
    
  END SUBROUTINE Region_CoordinateSystemGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the data points for a given user number in a region. \see OPENCMISS::Iron::cmfe_Region_DataPointsGet
  SUBROUTINE Region_DataPointsGet(region,userNumber,dataPoints,err,error,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: region !<A pointer to the region to get the data points for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the data points to get.
    TYPE(DataPointsType), POINTER :: dataPoints !<On exit, a pointer to the data points for the region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Region_DataPointsGet",err,error,*998)

    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated.",err,error,*998)
    IF(.NOT.region%REGION_FINISHED) CALL FlagError("Region data point sets is not associated.",err,error,*998)
    IF(ASSOCIATED(dataPoints)) CALL FlagError("Data points is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(region%dataPointSets)) CALL FlagError("Region data point sets is not associated.",err,error,*998)
    
    NULLIFY(dataPoints)
    CALL DataPointSets_UserNumberFind(region%dataPointSets,userNumber,dataPoints,err,error,*999)
    IF(.NOT.ASSOCIATED(dataPoints)) THEN
      localError="Data points with a user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " does not exist on region number "//TRIM(NumberToVString(region%USER_NUMBER,"*",err,error))//"."
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

  !>Returns a pointer to the equations set for a given user number in a region. \see OPENCMISS::Iron::cmfe_Region_EquationsSetGet
  SUBROUTINE Region_EquationsSetGet(region,userNumber,equationsSet,err,error,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: region !<A pointer to the region to get the equationsSet for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the equations set to get.
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<On exit, a pointer to the equations set for the region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Region_EquationsSetGet",err,error,*998)

    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated.",err,error,*998)
    IF(.NOT.region%REGION_FINISHED) CALL FlagError("Region has not been finished.",err,error,*998)
    IF(ASSOCIATED(equationsSet)) CALL FlagError("Equations set is already associated.",err,error,*998)
    
    NULLIFY(equationsSet)
    CALL EquationsSet_UserNumberFind(userNumber,region,equationsSet,err,error,*999)
    IF(.NOT.ASSOCIATED(equationsSet)) THEN
      localError="An equations set with a user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " does not exist on region number "//TRIM(NumberToVString(region%USER_NUMBER,"*",err,error))//"."
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
    TYPE(REGION_TYPE), POINTER :: region !<A pointer to the region to get the field for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the field to get.
    TYPE(FIELD_TYPE), POINTER :: field !<On exit, a pointer to the field for the region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Region_FieldGet",err,error,*998)

    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated.",err,error,*998)
    IF(.NOT.region%REGION_FINISHED) CALL FlagError("Region has not been finished.",err,error,*998)
    IF(ASSOCIATED(field)) CALL FlagError("Field is already associated.",err,error,*998)
    
    NULLIFY(field)
    CALL Field_UserNumberFind(userNumber,region,field,err,error,*999)
    IF(.NOT.ASSOCIATED(field)) THEN
      localError="A field with a user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " does not exist on region number "//TRIM(NumberToVString(region%USER_NUMBER,"*",err,error))//"."
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

  !>Returns a pointer to the generated mesh for a given user number in a region. \see OPENCMISS::Iron::cmfe_Region_GeneratedMeshGet
  SUBROUTINE Region_GeneratedMeshGet(region,userNumber,generatedMesh,err,error,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: region !<A pointer to the region to get the generated mesh for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the generated mesh to get.
    TYPE(GENERATED_MESH_TYPE), POINTER :: generatedMesh !<On exit, a pointer to the generated mesh for the region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Region_GeneratedMeshGet",err,error,*998)

    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated.",err,error,*998)
    IF(.NOT.region%REGION_FINISHED) CALL FlagError("Region has not been finished.",err,error,*998)
    IF(ASSOCIATED(generatedMesh)) CALL FlagError("Generated mesh is already associated.",err,error,*998)
    
    NULLIFY(generatedMesh)
    CALL GeneratedMesh_UserNumberFind(userNumber,region,generatedMesh,err,error,*999)
    IF(.NOT.ASSOCIATED(generatedMesh)) THEN
      localError="A generated mesh with a user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " does not exist on region number "//TRIM(NumberToVString(region%USER_NUMBER,"*",err,error))//"."
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
  SUBROUTINE Region_Get(userNumber,region,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the region to find
    TYPE(REGION_TYPE), POINTER :: region !<On exit, a pointer to the region with the specified user number if it exists. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Region_Get",err,error,*999)

    CALL Region_UserNumberFind(userNumber,region,err,error,*999)
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
    TYPE(REGION_TYPE), POINTER :: region !<A pointer to the parent region to get the interface for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the interface to get.
    TYPE(INTERFACE_TYPE), POINTER :: interface !<On exit, a pointer to the interface for the parent region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Region_InterfaceGet",err,error,*998)

    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated.",err,error,*998)
    IF(.NOT.region%REGION_FINISHED) CALL FlagError("Region has not been finished.",err,error,*998)
    IF(ASSOCIATED(interface)) CALL FlagError("Interface is already associated.",err,error,*998)
    
    NULLIFY(interface)
    CALL Interface_UserNumberFind(userNumber,region,interface,err,error,*999)
    IF(.NOT.ASSOCIATED(interface)) THEN
      localError="An interface with a user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " does not exist on region number "//TRIM(NumberToVString(region%USER_NUMBER,"*",err,error))//"."
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

  !>Returns a pointer to the mesh for a given user number in a region. \see OPENCMISS::Iron::cmfe_Region_MeshGet
  SUBROUTINE Region_MeshGet(region,userNumber,mesh,err,error,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: region !<A pointer to the region to get the mesh for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the mesh to get.
    TYPE(MESH_TYPE), POINTER :: mesh !<On exit, a pointer to the mesh for the region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Region_MeshGet",err,error,*998)

    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated.",err,error,*998)
    IF(.NOT.region%REGION_FINISHED) CALL FlagError("Region has not been finished.",err,error,*998)
    IF(ASSOCIATED(mesh)) CALL FlagError("Mesh is already associated.",err,error,*998)
    
    NULLIFY(mesh)
    CALL Mesh_UserNumberFind(userNumber,region,mesh,err,error,*999)
    IF(.NOT.ASSOCIATED(mesh)) THEN
      localError="A mesh with a user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " does not exist on region number "//TRIM(NumberToVString(region%USER_NUMBER,"*",err,error))//"."
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

  !>Returns a pointer to the nodes for a region. \see OPENCMISS::Iron::cmfe_RegionNodesGet
  SUBROUTINE Region_NodesGet(region,nodes,err,error,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: region !<A pointer to the region to get the nodes for
    TYPE(NodesType), POINTER :: nodes !<On exit, a pointer to the nodes for the region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Region_NodesGet",err,error,*998)

    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated.",err,error,*998)
    IF(.NOT.region%REGION_FINISHED) CALL FlagError("Region has not been finished.",err,error,*998)
    IF(ASSOCIATED(NODES)) CALL FlagError("Nodes is already associated.",err,error,*998)

    nodes=>region%nodes
    IF(.NOT.ASSOCIATED(nodes)) CALL FlagError("Region nodes is not associated.",err,error,*999)
       
    EXITS("Region_NodesGet")
    RETURN
999 NULLIFY(nodes)
998 ERRORSEXITS("Region_NodesGet",err,error)
    RETURN 1
    
  END SUBROUTINE Region_NodesGet
  
  !
  !================================================================================================================================
  !

  !>Finds and returns a pointer to the region with the given user number. If no region with that number exists region is left nullified.
  SUBROUTINE Region_UserNumberFind(userNumber,region,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the region to find
    TYPE(REGION_TYPE), POINTER :: region !<On exit, a pointer to the region with the specified user number if it exists. If no region exists with the specified user number a NULL pointer is returned. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: regionIdx
    TYPE(REGION_TYPE), POINTER :: worldRegion
    
    ENTERS("Region_UserNumberFind",err,error,*999)

    IF(ASSOCIATED(region)) CALL FlagError("Region is already associated.",err,error,*999)
    worldRegion=>regions%WORLD_REGION
    IF(.NOT.ASSOCIATED(worldRegion)) CALL FlagError("World region is not associated.",err,error,*999)
    
    NULLIFY(region)
    IF(userNumber==0) THEN
      region=>worldRegion
    ELSE
      IF(ASSOCIATED(worldRegion%SUB_REGIONS)) THEN
        DO regionIdx=1,worldRegion%NUMBER_OF_SUB_REGIONS
          CALL Region_UserNumberFindPtr(userNumber,worldRegion%SUB_REGIONS(regionIdx)%ptr,region,err,error,*999)
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
    TYPE(REGION_TYPE), POINTER :: startRegion !<A pointer to the region to start the search from
    TYPE(REGION_TYPE), POINTER :: region !<On exit, a pointer to the region with the specified user number if it exists. If no region exists with the specified user number a NULL pointer is returned.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: regionIdx

    ENTERS("Region_UserNumberFindPtr",err,error,*999)

    IF(.NOT.ASSOCIATED(startRegion)) CALL FlagError("Start region is not associated",err,error,*999)
    IF(ASSOCIATED(region)) CALL FlagError("Region is already associated.",err,error,*999)

    NULLIFY(region)
    IF(startRegion%USER_NUMBER==userNumber) THEN
      region=>startRegion
    ELSE
      IF(ASSOCIATED(startRegion%SUB_REGIONS)) THEN
        DO regionIdx=1,startRegion%NUMBER_OF_SUB_REGIONS
          CALL Region_UserNumberFindPtr(userNumber,startRegion%SUB_REGIONS(regionIdx)%ptr,region,err,error,*999)
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
    TYPE(REGION_TYPE), POINTER :: region !<A pointer to the region to get the user number for
    INTEGER(INTG), INTENT(OUT) :: userNumber !<On exit, the user number of the region.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Region_UserNumberGet",err,error,*999)

    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated.",err,error,*999)

    userNumber=region%USER_NUMBER
  
    EXITS("Region_UserNumberGet")
    RETURN
999 ERRORSEXITS("Region_UserNumberGet",err,error)
    RETURN 1
    
  END SUBROUTINE Region_UserNumberGet

  !
  !================================================================================================================================
  !

END MODULE RegionAccessRoutines
