!> \file
!> \author Chris Bradley
!> \brief This module contains all interface access method routines.
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

!> This module contains all interface access method routines.
MODULE InterfaceAccessRoutines
  
  USE BaseRoutines
  USE DataPointAccessRoutines
  USE FieldAccessRoutines
  USE InterfaceConditionAccessRoutines
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

  INTERFACE INTERFACE_COORDINATE_SYSTEM_GET
    MODULE PROCEDURE Interface_CoordinateSystemGet
  END INTERFACE INTERFACE_COORDINATE_SYSTEM_GET

  INTERFACE INTERFACE_DATA_POINTS_GET
    MODULE PROCEDURE Interface_DataPointsGet
  END INTERFACE INTERFACE_DATA_POINTS_GET

  INTERFACE INTERFACE_NODES_GET
    MODULE PROCEDURE Interface_NodesGet
  END INTERFACE INTERFACE_NODES_GET

  INTERFACE INTERFACE_USER_NUMBER_FIND
    MODULE PROCEDURE Interface_UserNumberFind
  END INTERFACE INTERFACE_USER_NUMBER_FIND

  PUBLIC Interface_AssertIsFinished,Interface_AssertNotFinished

  PUBLIC Interface_CoordinateSystemGet

  PUBLIC INTERFACE_COORDINATE_SYSTEM_GET

  PUBLIC Interface_CoupledMeshGet

  PUBLIC Interface_DataPointsGet

  PUBLIC INTERFACE_DATA_POINTS_GET

  PUBLIC Interface_FieldGet

  PUBLIC Interface_InterfaceConditionGet

  PUBLIC Interface_InterfacesGet

  PUBLIC Interface_MeshGet

  PUBLIC Interface_MeshConnectivityGet

  PUBLIC Interface_MeshesGet

  PUBLIC Interface_NodesGet

  PUBLIC INTERFACE_NODES_GET

  PUBLIC Interface_ParentRegionGet

  PUBLIC Interface_PointsConnectivityGet

  PUBLIC Interface_UserNumberFind

  PUBLIC INTERFACE_USER_NUMBER_FIND

  PUBLIC InterfaceMeshConnectivity_AssertIsFinished,InterfaceMeshConnectivity_AssertNotFinished

  PUBLIC InterfaceMeshConnectivity_BasisGet

  PUBLIC InterfaceMeshConnectivity_InterfaceGet

  PUBLIC InterfaceMeshConnectivity_InterfaceMeshGet

  PUBLIC InterfacePointsConnectivity_AssertIsFinished,InterfacePointsConnectivity_AssertNotFinished

  PUBLIC InterfacePointsConnectivity_DataPointsGet

  PUBLIC InterfacePointsConnectivity_InterfaceGet

  PUBLIC InterfacePointsConnectivity_InterfaceMeshGet

CONTAINS

  !
  !=================================================================================================================================
  !

  !>Assert that a interface has been finished
  SUBROUTINE Interface_AssertIsFinished(interface,err,error,*)

    !Argument Variables
    TYPE(InterfaceType), POINTER, INTENT(INOUT) :: interface !<The work group to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Interface_AssertIsFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(interface)) CALL FlagError("Interface is not associated.",err,error,*999)

    IF(.NOT.interface%interfaceFinished) THEN
      localError="Interface number "//TRIM(NumberToVString(interface%userNumber,"*",err,error))
      IF(ASSOCIATED(INTERFACE%parentRegion)) localError=localError// &
        & " of parent region number "//TRIM(NumberToVString(INTERFACE%parentRegion%userNumber,"*",err,error))
       localError=localError//" has not been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Interface_AssertIsFinished")
    RETURN
999 ERRORSEXITS("Interface_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE Interface_AssertIsFinished

  !
  !=================================================================================================================================
  !

  !>Assert that a interface has not been finished
  SUBROUTINE Interface_AssertNotFinished(interface,err,error,*)

    !Argument Variables
    TYPE(InterfaceType), POINTER, INTENT(INOUT) :: interface !<The work group to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Interface_AssertNotFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(interface)) CALL FlagError("Interface is not associated.",err,error,*999)

    IF(interface%interfaceFinished) THEN
      localError="Interface number "//TRIM(NumberToVString(interface%userNumber,"*",err,error))
      IF(ASSOCIATED(INTERFACE%parentRegion)) localError=localError// &
        & " of parent region number "//TRIM(NumberToVString(INTERFACE%parentRegion%userNumber,"*",err,error))
      localError=localError//" has already been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Interface_AssertNotFinished")
    RETURN
999 ERRORSEXITS("Interface_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE Interface_AssertNotFinished

  !
  !================================================================================================================================
  !
  
  !>Returns the coordinate system of an interface. \see OPENCMISS::Iron::cmfe_Interface_CoordinateSystemGet
  SUBROUTINE Interface_CoordinateSystemGet(interface,coordinateSystem,err,error,*)

    !Argument variables
    TYPE(InterfaceType), POINTER :: interface !<A pointer to the interface to get the coordinate system for
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem !<On exit, the coordinate system for the specified interface. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Interface_CoordinateSystemGet",err,error,*999)

    IF(.NOT.ASSOCIATED(INTERFACE)) CALL FlagError("Interface is not associated.",err,error,*999)
    IF(.NOT.INTERFACE%interfaceFinished) CALL FlagError("Interface has not been finished.",err,error,*999)
    IF(ASSOCIATED(coordinateSystem)) CALL FlagError("Coordinate system is already associated.",ERR,ERROR,*999)

    coordinateSystem=>interface%coordinateSystem
    IF(.NOT.ASSOCIATED(coordinateSystem)) THEN
      localError="Coordinate system is not associated for interface number "// &
        & TRIM(NumberToVString(interface%userNumber,"*",err,error))
      IF(ASSOCIATED(INTERFACE%parentRegion)) localError=localError// &
        & " of parent region number "//TRIM(NumberToVString(INTERFACE%parentRegion%userNumber,"*",err,error))
      localError=localError//"."
     CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Interface_CoordinateSystemGet")
    RETURN
999 ERRORSEXITS("Interface_CoordinateSystemGet",err,error)
    RETURN 1
    
  END SUBROUTINE Interface_CoordinateSystemGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to a coupled mesh in an interface. 
  SUBROUTINE Interface_CoupledMeshGet(INTERFACE,coupledMeshIndex,coupledMesh,err,error,*)

    !Argument variables
    TYPE(InterfaceType), POINTER :: interface !<A pointer to the interface to get the coupled mesh for
    INTEGER(INTG), INTENT(IN) :: coupledMeshIndex !<The index of the coupled mesh to get.
    TYPE(MeshType), POINTER :: coupledMesh !<On exit, a pointer to the specified coupled mesh for the interface. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Interface_CoupledMeshGet",err,error,*998)

    IF(ASSOCIATED(coupledMesh)) CALL FlagError("Coupled mesh is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interface)) CALL FlagError("Interface is not associated.",err,error,*999)
    IF(coupledMeshIndex<1.OR.coupledMeshIndex>interface%numberOfCoupledMeshes) THEN
      localError="The specified coupled mesh index of "//TRIM(NumberToVString(coupledMeshIndex,"*",err,error))// &
        & " is invalid. The coupled mesh index should be >=1 and <= "// &
        & TRIM(NumberToVString(INTERFACE%numberOfCoupledMeshes,"*",err,error))// &
        & " for interface number "//TRIM(NumberToVString(INTERFACE%userNumber,"*",err,error))
      IF(ASSOCIATED(INTERFACE%parentRegion)) localError=localError// &
        & " of parent region number "//TRIM(NumberToVString(INTERFACE%parentRegion%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(INTERFACE%coupledMeshes)) THEN
      localError="Coupled meshes is not allocated for interface number "// &
        & TRIM(NumberToVString(INTERFACE%userNumber,"*",err,error))
      IF(ASSOCIATED(INTERFACE%parentRegion)) localError=localError// &
        & " of parent region number "//TRIM(NumberToVString(INTERFACE%parentRegion%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    coupledMesh=>interface%coupledMeshes(coupledMeshIndex)%ptr
    IF(.NOT.ASSOCIATED(coupledMesh)) THEN
      localError="The coupled mesh for coupled mesh index "//TRIM(NumberToVString(coupledMeshIndex,"*",err,error))// &
        & " is not associated for interface number "//TRIM(NumberToVString(INTERFACE%userNumber,"*",err,error))
      IF(ASSOCIATED(INTERFACE%parentRegion)) localError=localError// &
        & " of parent region number "//TRIM(NumberToVString(INTERFACE%parentRegion%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
      
    EXITS("Interface_CoupledMeshGet")
    RETURN
999 NULLIFY(coupledMesh)
998 ERRORSEXITS("Interface_CoupledMeshGet",err,error)
    RETURN 1
    
  END SUBROUTINE Interface_CoupledMeshGet

  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the data points for a given user number in an interface. \see OPENCMISS::Iron::cmfe_Interface_DataPointsGet
  SUBROUTINE Interface_DataPointsGet(interface,userNumber,dataPoints,err,error,*)

    !Argument variables
    TYPE(InterfaceType), POINTER :: interface !<A pointer to the interface to get the data points for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the data points to get.
    TYPE(DataPointsType), POINTER :: dataPoints !<On exit, a pointer to the data points for the interface. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Interface_DataPointsGet",err,error,*998)

    IF(.NOT.ASSOCIATED(interface)) CALL FlagError("Interface is not associated.",err,error,*998)
    IF(.NOT.interface%interfaceFinished) CALL FlagError("Interface has not been finished.",err,error,*998)
    IF(ASSOCIATED(dataPoints)) CALL FlagError("Data points is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interface%dataPointSets)) CALL FlagError("Interface data point sets is not associated.",err,error,*998)

    NULLIFY(dataPoints)
    CALL DataPointSets_UserNumberFind(INTERFACE%dataPointSets,userNumber,dataPoints,err,error,*999)
    IF(.NOT.ASSOCIATED(dataPoints)) THEN
      localError="Data points with a user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " do not exist on interface number "//TRIM(NumberToVString(interface%userNumber,"*",err,error))
      IF(ASSOCIATED(INTERFACE%parentRegion)) localError=localError// &
        & " of parent region number "//TRIM(NumberToVString(INTERFACE%parentRegion%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
      
    EXITS("Interface_DataPointsGet")
    RETURN
999 NULLIFY(dataPoints)
998 ERRORSEXITS("Interface_DataPointsGet",err,error)
    RETURN 1
    
  END SUBROUTINE Interface_DataPointsGet    

  !
  !================================================================================================================================
  !

  !>Returns a pointer to a field for a given user number in an interface. \see OPENCMISS::Iron::cmfe_Interface_FieldGet
  SUBROUTINE Interface_FieldGet(interface,userNumber,field,err,error,*)

    !Argument variables
    TYPE(InterfaceType), POINTER :: INTERFACE !<A pointer to the interface to get the field for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the field to get.
    TYPE(FieldType), POINTER :: field !<On exit, a pointer to the field for the interface. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Interface_FieldGet",err,error,*998)

    IF(.NOT.ASSOCIATED(interface)) CALL FlagError("Interface is not associated.",err,error,*998)
    IF(.NOT.interface%interfaceFinished) CALL FlagError("Interface has not been finished.",err,error,*998)
    IF(ASSOCIATED(field)) CALL FlagError("Field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interface%fields)) CALL FlagError("Interface fields is not associated.",err,error,*998)

    NULLIFY(field)
    CALL Field_UserNumberFind(userNumber,interface,field,err,error,*999)
    IF(.NOT.ASSOCIATED(field)) THEN
      localError="A field with a user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " do not exist on interface number "//TRIM(NumberToVString(interface%userNumber,"*",err,error))
      IF(ASSOCIATED(INTERFACE%parentRegion)) localError=localError// &
        & " of parent region number "//TRIM(NumberToVString(INTERFACE%parentRegion%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("Interface_FieldGet")
    RETURN
999 NULLIFY(field)
998 ERRORSEXITS("Interface_FieldGet",err,error)
    RETURN 1
    
  END SUBROUTINE Interface_FieldGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to a interface condition for a given user number in an interface. \see OPENCMISS::Iron::cmfe_Interface_InterfaceConditionGet
  SUBROUTINE Interface_InterfaceConditionGet(interface,userNumber,interfaceCondition,err,error,*)

    !Argument variables
    TYPE(InterfaceType), POINTER :: INTERFACE !<A pointer to the interface to get the interface condition for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the interface condition to get.
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<On exit, a pointer to the interface condition for the interface. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Interface_InterfaceConditionGet",err,error,*998)

    IF(.NOT.ASSOCIATED(interface)) CALL FlagError("Interface is not associated.",err,error,*998)
    IF(.NOT.interface%interfaceFinished) CALL FlagError("Interface has not been finished.",err,error,*998)
    IF(ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interface%interfaceConditions)) &
      & CALL FlagError("Interface interface conditions is not associated.",err,error,*998)

    NULLIFY(interfaceCondition)
    CALL InterfaceCondition_UserNumberFind(userNumber,interface,interfaceCondition,err,error,*999)
    IF(.NOT.ASSOCIATED(interfaceCondition)) THEN
      localError="An interface Condition with a user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " do not exist on interface number "//TRIM(NumberToVString(interface%userNumber,"*",err,error))
      IF(ASSOCIATED(INTERFACE%parentRegion)) localError=localError// &
        & " of parent region number "//TRIM(NumberToVString(INTERFACE%parentRegion%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("Interface_InterfaceConditionGet")
    RETURN
999 NULLIFY(interfaceCondition)
998 ERRORSEXITS("Interface_InterfaceConditionGet",err,error)
    RETURN 1
    
  END SUBROUTINE Interface_InterfaceConditionGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the interfaces for a given interface. 
  SUBROUTINE Interface_InterfacesGet(INTERFACE,interfaces,err,error,*)

    !Argument variables
    TYPE(InterfaceType), POINTER :: interface !<A pointer to the interface to get the mesh for
    TYPE(InterfacesType), POINTER :: interfaces !<On return, a pointer to the interfaces for the interface. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Interface_InterfacesGet",err,error,*998)

    IF(ASSOCIATED(interfaces)) CALL FlagError("Interfaces is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interface)) CALL FlagError("Interface is not associated.",err,error,*998)

    interfaces=>interface%interfaces
    IF(.NOT.ASSOCIATED(interfaces)) THEN
      localError="Interfaces is not associated for interface number "// &
        & TRIM(NumberToVString(interface%userNumber,"*",err,error))
      IF(ASSOCIATED(INTERFACE%parentRegion)) localError=localError// &
        & " of parent region number "//TRIM(NumberToVString(INTERFACE%parentRegion%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Interface_InterfacesGet")
    RETURN
999 NULLIFY(interfaces)
998 ERRORSEXITS("Interface_InterfacesGet",err,error)
    RETURN 1
    
  END SUBROUTINE Interface_InterfacesGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the mesh for a given user number in a interface. \see OPENCMISS::Iron::cmfe_Interface_MeshGet
  SUBROUTINE Interface_MeshGet(interface,userNumber,mesh,err,error,*)

    !Argument variables
    TYPE(InterfaceType), POINTER :: interface !<A pointer to the interface to get the mesh for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the mesh to get.
    TYPE(MeshType), POINTER :: mesh !<On exit, a pointer to the mesh for the interface. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Interface_MeshGet",err,error,*998)

    IF(.NOT.ASSOCIATED(interface)) CALL FlagError("Interface is not associated.",err,error,*998)
    IF(.NOT.interface%interfaceFinished) CALL FlagError("Interface has not been finished.",err,error,*998)
    IF(ASSOCIATED(mesh)) CALL FlagError("Mesh is already associated.",err,error,*998)
    
    NULLIFY(mesh)
    CALL Mesh_UserNumberFind(userNumber,interface,mesh,err,error,*999)
    IF(.NOT.ASSOCIATED(mesh)) THEN
      localError="A mesh with a user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " does not exist on interface number "//TRIM(NumberToVString(interface%userNumber,"*",err,error))
      IF(ASSOCIATED(INTERFACE%parentRegion)) localError=localError// &
        & " of parent region number "//TRIM(NumberToVString(INTERFACE%parentRegion%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Interface_MeshGet")
    RETURN
999 NULLIFY(mesh)
998 ERRORSEXITS("Interface_MeshGet",err,error)
    RETURN 1
    
  END SUBROUTINE Interface_MeshGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the mesh connectivity for a given interface.
  SUBROUTINE Interface_MeshConnectivityGet(INTERFACE,meshConnectivity,err,error,*)

    !Argument variables
    TYPE(InterfaceType), POINTER :: interface !<A pointer to the interface to get the mesh connectivity for
    TYPE(InterfaceMeshConnectivityType), POINTER :: meshConnectivity !<On exit, a pointer to the mesh connectivity for the interface. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Interface_MeshConnectivityGet",err,error,*998)

    IF(ASSOCIATED(meshConnectivity)) CALL FlagError("Mesh connectivity is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interface)) CALL FlagError("Interface is not associated.",err,error,*998)
    IF(.NOT.interface%interfaceFinished) CALL FlagError("Interface has not been finished.",err,error,*998)
    
    meshConnectivity=>interface%meshConnectivity
    IF(.NOT.ASSOCIATED(meshConnectivity)) THEN
      localError="Mesh connectivity is not associated on intereface number "// &
        & TRIM(NumberToVString(interface%userNumber,"*",err,error))
      IF(ASSOCIATED(INTERFACE%parentRegion)) localError=localError// &
        & " of parent region number "//TRIM(NumberToVString(INTERFACE%parentRegion%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Interface_MeshConnectivityGet")
    RETURN
999 NULLIFY(meshConnectivity)
998 ERRORSEXITS("Interface_MeshConnectivityGet",err,error)
    RETURN 1
    
  END SUBROUTINE Interface_MeshConnectivityGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the meshes for a interface. 
  SUBROUTINE Interface_MeshesGet(INTERFACE,meshes,err,error,*)

    !Argument variables
    TYPE(InterfaceType), POINTER :: interface !<A pointer to the interface to get the meshes for
    TYPE(MeshesType), POINTER :: meshes !<On exit, a pointer to the meshes for the interface. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Interface_MeshesGet",err,error,*998)

    IF(ASSOCIATED(meshes)) CALL FlagError("Meshes is already associated.",err,error,*998)
    CALL Interface_AssertIsFinished(interface,err,error,*998)
 
    meshes=>interface%meshes
    IF(.NOT.ASSOCIATED(meshes)) THEN
      localError="Meshes is not associated for interface number "//TRIM(NumberToVString(INTERFACE%userNumber,"*",err,error))
      IF(ASSOCIATED(INTERFACE%parentRegion)) localError=localError// &
        & " of parent region number "//TRIM(NumberToVString(INTERFACE%parentRegion%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("Interface_MeshesGet")
    RETURN
999 NULLIFY(meshes)
998 ERRORSEXITS("Interface_MeshesGet",err,error)
    RETURN 1
    
  END SUBROUTINE Interface_MeshesGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the nodes for a interface. \see OPENCMISS::Iron::cmfe_InterfaceNodesGet
  SUBROUTINE Interface_NodesGet(interface,nodes,err,error,*)

    !Argument variables
    TYPE(InterfaceType), POINTER :: interface !<A pointer to the interface to get the nodes for
    TYPE(NodesType), POINTER :: nodes !<On exit, a pointer to the nodes for the interface. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Interface_NodesGet",err,error,*998)

    IF(.NOT.ASSOCIATED(interface)) CALL FlagError("Interface is not associated.",err,error,*998)
    IF(.NOT.interface%interfaceFinished) CALL FlagError("Interface has not been finished.",err,error,*998)
    IF(ASSOCIATED(nodes)) CALL FlagError("Nodes is already associated.",err,error,*998)

    nodes=>interface%nodes
    IF(.NOT.ASSOCIATED(nodes)) THEN
      localError="Nodes is not associated for interface number "//TRIM(NumberToVString(INTERFACE%userNumber,"*",err,error))
      IF(ASSOCIATED(INTERFACE%parentRegion)) localError=localError// &
        & " of parent region number "//TRIM(NumberToVString(INTERFACE%parentRegion%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("Interface_NodesGet")
    RETURN
999 NULLIFY(nodes)
998 ERRORSEXITS("Interface_NodesGet",err,error)
    RETURN 1
    
  END SUBROUTINE Interface_NodesGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the parent region for an interface. 
  SUBROUTINE Interface_ParentRegionGet(INTERFACE,parentRegion,err,error,*)

    !Argument variables
    TYPE(InterfaceType), POINTER :: interface !<A pointer to the interface to get the parent region for
    TYPE(RegionType), POINTER :: parentRegion !<On exit, a pointer to the parent region for the interface. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Interface_ParentRegionGet",err,error,*998)

    IF(ASSOCIATED(parentRegion)) CALL FlagError("Parent region is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interface)) CALL FlagError("Interface is not associated.",err,error,*999)
 
    parentRegion=>interface%parentRegion
    IF(.NOT.ASSOCIATED(parentRegion)) THEN
      localError="The parent region for interface number "//TRIM(NumberToVString(interface%userNumber,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("Interface_ParentRegionGet")
    RETURN
999 NULLIFY(parentRegion)
998 ERRORSEXITS("Interface_ParentRegionGet",err,error)
    RETURN 1
    
  END SUBROUTINE Interface_ParentRegionGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the points connectivity for a given interface.
  SUBROUTINE Interface_PointsConnectivityGet(INTERFACE,pointsConnectivity,err,error,*)

    !Argument variables
    TYPE(InterfaceType), POINTER :: interface !<A pointer to the interface to get the points connectivity for
    TYPE(InterfacePointsConnectivityType), POINTER :: pointsConnectivity !<On exit, a pointer to the points connectivity for the interface. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Interface_PointConnectivityGet",err,error,*998)

    IF(ASSOCIATED(pointsConnectivity)) CALL FlagError("Point connectivity is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interface)) CALL FlagError("Interface is not associated.",err,error,*999)
    IF(.NOT.interface%interfaceFinished) CALL FlagError("Interface has not been finished.",err,error,*999)
    
    pointsConnectivity=>interface%pointsConnectivity
    IF(.NOT.ASSOCIATED(pointsConnectivity)) THEN
      localError="Points connectivity is not associated on intereface number "// &
        & TRIM(NumberToVString(interface%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Interface_PointsConnectivityGet")
    RETURN
999 NULLIFY(pointsConnectivity)
998 ERRORSEXITS("Interface_PointsConnectivityGet",err,error)
    RETURN 1
    
  END SUBROUTINE Interface_PointsConnectivityGet

  !
  !================================================================================================================================
  !

  !>Finds and returns a pointer to the interface identified by user number in the given parent region. If no interface with that user number exists null is returned.
  SUBROUTINE Interface_UserNumberFind(userNumber,parentRegion,interface,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number to find.
    TYPE(RegionType), POINTER :: parentRegion !<The parent region to find the interface in    
    TYPE(InterfaceType), POINTER :: interface !<On return a pointer to the interface with the given user number. If no interface with that user number exists then the pointer is returned as NULL. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: interfaceIdx
    TYPE(VARYING_STRING) :: localError

    ENTERS("Interface_UserNumberFind",err,error,*999)

    !Check input arguments
    IF(.NOT.ASSOCIATED(parentRegion)) CALL FlagError("Parent region is not associated.",err,error,*999)
    IF(ASSOCIATED(interface)) CALL FlagError("Interface is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(parentRegion%interfaces)) THEN
      localError="The interfaces on parent region number "// &
        & TRIM(NumberToVString(parentRegion%userNumber,"*",err,error))//" are not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    !Get the interface from the user number
    NULLIFY(INTERFACE)
    IF(ALLOCATED(parentRegion%interfaces%interfaces)) THEN
      DO interfaceIdx=1,parentRegion%interfaces%numberOfInterfaces
        IF(ASSOCIATED(parentRegion%interfaces%interfaces(interfaceIdx)%ptr)) THEN
          IF(parentRegion%interfaces%interfaces(interfaceIdx)%ptr%userNumber==userNumber) THEN
            INTERFACE=>parentRegion%interfaces%interfaces(interfaceIdx)%ptr
            EXIT
          ENDIF
        ELSE
          localError="The interface pointer in interfaces is not associated for interface index "// &
            & TRIM(NumberToVString(interfaceIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDDO !interfaceIdx
    ENDIF
   
    EXITS("Interface_UserNumberFind")
    RETURN
999 ERRORSEXITS("Interface_UserNumberFind",err,error)
    RETURN 1
    
  END SUBROUTINE Interface_UserNumberFind

  !
  !=================================================================================================================================
  !

  !>Assert that a interface mesh connectivity has been finished
  SUBROUTINE InterfaceMeshConnectivity_AssertIsFinished(interfaceMeshConnectivity,err,error,*)

    !Argument Variables
    TYPE(InterfaceMeshConnectivityType), POINTER, INTENT(INOUT) :: interfaceMeshConnectivity !<The interface mesh connectivity to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("InterfaceMeshConnectivity_AssertIsFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceMeshConnectivity)) CALL FlagError("Interface mesh connectivity is not associated.",err,error,*999)

    IF(.NOT.interfaceMeshConnectivity%meshConnectivityFinished) THEN
      localError="Interface mesh connectivity "
      IF(ASSOCIATED(interfaceMeshConnectivity%INTERFACE)) THEN
        localError=localError//" for interface number "// &
          & TRIM(NumberToVString(interfaceMeshConnectivity%interface%userNumber,"*",err,error))
        IF(ASSOCIATED(interfaceMeshConnectivity%interface%parentRegion)) localError=localError//" of parent region number "// &
          & TRIM(NumberToVString(interfaceMeshConnectivity%INTERFACE%parentRegion%userNumber,"*",err,error))
      ENDIF
      localError=localError//" has not been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("InterfaceMeshConnectivity_AssertIsFinished")
    RETURN
999 ERRORSEXITS("InterfaceMeshConnectivity_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMeshConnectivity_AssertIsFinished

  !
  !=================================================================================================================================
  !

  !>Assert that a interface mesh connectivity has not been finished
  SUBROUTINE InterfaceMeshConnectivity_AssertNotFinished(interfaceMeshConnectivity,err,error,*)

    !Argument Variables
    TYPE(InterfaceMeshConnectivityType), POINTER, INTENT(INOUT) :: interfaceMeshConnectivity !<The interface mesh connectivity to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("InterfaceMeshConnectivity_AssertNotFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceMeshConnectivity)) CALL FlagError("Interface mesh connectivity is not associated.",err,error,*999)

    IF(interfaceMeshConnectivity%meshConnectivityFinished) THEN
      localError="Interface mesh connectivity "
      IF(ASSOCIATED(interfaceMeshConnectivity%INTERFACE)) THEN
        localError=localError//" for interface number "// &
          & TRIM(NumberToVString(interfaceMeshConnectivity%interface%userNumber,"*",err,error))
        IF(ASSOCIATED(interfaceMeshConnectivity%interface%parentRegion)) localError=localError//" of parent region number "// &
          & TRIM(NumberToVString(interfaceMeshConnectivity%INTERFACE%parentRegion%userNumber,"*",err,error))
      ENDIF
      localError=localError//" has already been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("InterfaceMeshConnectivity_AssertNotFinished")
    RETURN
999 ERRORSEXITS("InterfaceMeshConnectivity_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMeshConnectivity_AssertNotFinished

  !
  !=================================================================================================================================
  !

  !>Returns the basis for an interface mesh connectivity
  SUBROUTINE InterfaceMeshConnectivity_BasisGet(interfaceMeshConnectivity,basis,err,error,*)

    !Argument Variables
    TYPE(InterfaceMeshConnectivityType), POINTER, INTENT(INOUT) :: interfaceMeshConnectivity !<The interface mesh connectivity to get the basis for
    TYPE(BasisType), POINTER, INTENT(OUT) :: basis !<On return, the basis for the interface mesh connectivity. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMeshConnectivity_BasisGet",err,error,*998)

    IF(ASSOCIATED(basis)) CALL FlagError("Basis is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMeshConnectivity)) &
      & CALL FlagError("Interface mesh connectivity is not associated.",err,error,*999)

    basis=>interfaceMeshConnectivity%basis
    IF(.NOT.ASSOCIATED(basis)) &
      & CALL FlagError("The interface mesh connectivity basis is not associated.",err,error,*999)
    
    EXITS("InterfaceMeshConnectivity_BasisGet")
    RETURN
999 NULLIFY(basis)
998 ERRORSEXITS("InterfaceMeshConnectivity_BasisGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMeshConnectivity_BasisGet

  !
  !=================================================================================================================================
  !

  !>Returns the interface for an interface mesh connectivity
  SUBROUTINE InterfaceMeshConnectivity_InterfaceGet(interfaceMeshConnectivity,interface,err,error,*)

    !Argument Variables
    TYPE(InterfaceMeshConnectivityType), POINTER, INTENT(INOUT) :: interfaceMeshConnectivity !<The interface mesh connectivity to get the interface for
    TYPE(InterfaceType), POINTER, INTENT(OUT) :: INTERFACE !<On return, the interface for the interface mesh connectivity. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMeshConnectivity_InterfaceGet",err,error,*998)

    IF(ASSOCIATED(INTERFACE)) CALL FlagError("Interface is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMeshConnectivity)) &
      & CALL FlagError("Interface mesh connectivity is not associated.",err,error,*999)

    INTERFACE=>interfaceMeshConnectivity%INTERFACE
    IF(.NOT.ASSOCIATED(INTERFACE)) &
      & CALL FlagError("The interface mesh connectivity interface is not associated.",err,error,*999)
    
    EXITS("InterfaceMeshConnectivity_InterfaceGet")
    RETURN
999 NULLIFY(interface)
998 ERRORSEXITS("InterfaceMeshConnectivity_InterfaceGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMeshConnectivity_InterfaceGet

  !
  !=================================================================================================================================
  !

  !>Returns the interface mesh for an interface mesh connectivity
  SUBROUTINE InterfaceMeshConnectivity_InterfaceMeshGet(interfaceMeshConnectivity,interfaceMesh,err,error,*)

    !Argument Variables
    TYPE(InterfaceMeshConnectivityType), POINTER, INTENT(INOUT) :: interfaceMeshConnectivity !<The interface mesh connectivity to get the interface mesh for
    TYPE(MeshType), POINTER, INTENT(OUT) :: interfaceMesh !<On return, the interface mesh for the interface mesh connectivity. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMeshConnectivity_InterfaceMeshGet",err,error,*998)

    IF(ASSOCIATED(interfaceMesh)) CALL FlagError("Interface mesh is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMeshConnectivity)) &
      & CALL FlagError("Interface mesh connectivity is not associated.",err,error,*999)

    interfaceMesh=>interfaceMeshConnectivity%interfaceMesh
    IF(.NOT.ASSOCIATED(interfaceMesh)) &
      & CALL FlagError("The interface mesh connectivity interface mesh is not associated.",err,error,*999)
    
    EXITS("InterfaceMeshConnectivity_InterfaceMeshGet")
    RETURN
999 NULLIFY(interfaceMesh)
998 ERRORSEXITS("InterfaceMeshConnectivity_InterfaceMeshGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMeshConnectivity_InterfaceMeshGet

  !
  !=================================================================================================================================
  !

  !>Assert that a interface points connectivity has been finished
  SUBROUTINE InterfacePointsConnectivity_AssertIsFinished(interfacePointsConnectivity,err,error,*)

    !Argument Variables
    TYPE(InterfacePointsConnectivityType), POINTER, INTENT(INOUT) :: interfacePointsConnectivity !<The interface points connectivity to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("InterfacePointsConnectivity_AssertIsFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(interfacePointsConnectivity)) &
      & CALL FlagError("Interface points connectivity is not associated.",err,error,*999)

    IF(.NOT.interfacePointsConnectivity%pointsConnectivityFinished) THEN
      localError="Interface points connectivity "
      IF(ASSOCIATED(interfacePointsConnectivity%INTERFACE)) THEN
        localError=localError//" for interface number "// &
          & TRIM(NumberToVString(interfacePointsConnectivity%interface%userNumber,"*",err,error))
        IF(ASSOCIATED(interfacePointsConnectivity%interface%parentRegion)) localError=localError//" of parent region number "// &
          & TRIM(NumberToVString(interfacePointsConnectivity%INTERFACE%parentRegion%userNumber,"*",err,error))
      ENDIF
      localError=localError//" has not been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("InterfacePointsConnectivity_AssertIsFinished")
    RETURN
999 ERRORS("InterfacePointsConnectivity_AssertIsFinished",err,error)
    EXITS("InterfacePointsConnectivity_AssertIsFinished")
    RETURN 1
    
  END SUBROUTINE InterfacePointsConnectivity_AssertIsFinished

  !
  !=================================================================================================================================
  !

  !>Assert that a interface points connectivity has not been finished
  SUBROUTINE InterfacePointsConnectivity_AssertNotFinished(interfacePointsConnectivity,err,error,*)

    !Argument Variables
    TYPE(InterfacePointsConnectivityType), POINTER, INTENT(INOUT) :: interfacePointsConnectivity !<The interface points connectivity to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("InterfacePointsConnectivity_AssertNotFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(interfacePointsConnectivity)) &
      & CALL FlagError("Interface points connectivity is not associated.",err,error,*999)

    IF(interfacePointsConnectivity%pointsConnectivityFinished) THEN
      localError="Interface points connectivity "
      IF(ASSOCIATED(interfacePointsConnectivity%INTERFACE)) THEN
        localError=localError//" for interface number "// &
          & TRIM(NumberToVString(interfacePointsConnectivity%interface%userNumber,"*",err,error))
        IF(ASSOCIATED(interfacePointsConnectivity%interface%parentRegion)) localError=localError//" of parent region number "// &
          & TRIM(NumberToVString(interfacePointsConnectivity%INTERFACE%parentRegion%userNumber,"*",err,error))
      ENDIF
      localError=localError//" has already been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("InterfacePointsConnectivity_AssertNotFinished")
    RETURN
999 ERRORS("InterfacePointsConnectivity_AssertNotFinished",err,error)
    EXITS("InterfacePointsConnectivity_AssertNotFinished")
    RETURN 1   
    
  END SUBROUTINE InterfacePointsConnectivity_AssertNotFinished

  !
  !=================================================================================================================================
  !

  !>Returns the data points for an interface points connectivity
  SUBROUTINE InterfacePointsConnectivity_DataPointsGet(interfacePointsConnectivity,dataPoints,err,error,*)

    !Argument Variables
    TYPE(InterfacePointsConnectivityType), POINTER, INTENT(INOUT) :: interfacePointsConnectivity !<The interface points connectivity to get the data points for
    TYPE(DataPointsType), POINTER, INTENT(OUT) :: dataPoints !<On return, the data points for the interface points connectivity. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfacePointsConnectivity_DataPointsGet",err,error,*998)

    IF(ASSOCIATED(dataPoints)) CALL FlagError("Data points is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfacePointsConnectivity)) &
      & CALL FlagError("Interface points connectivity is not associated.",err,error,*999)

    dataPoints=>interfacePointsConnectivity%dataPoints
    IF(.NOT.ASSOCIATED(dataPoints)) &
      & CALL FlagError("The interface points connectivity data points is not associated.",err,error,*999)
    
    EXITS("InterfacePointsConnectivity_DataPointsGet")
    RETURN
999 NULLIFY(dataPoints)
998 ERRORSEXITS("InterfacePointsConnectivity_DataPointsGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfacePointsConnectivity_DataPointsGet

  !
  !=================================================================================================================================
  !

  !>Returns the interface for an interface points connectivity
  SUBROUTINE InterfacePointsConnectivity_InterfaceGet(interfacePointsConnectivity,interface,err,error,*)

    !Argument Variables
    TYPE(InterfacePointsConnectivityType), POINTER, INTENT(INOUT) :: interfacePointsConnectivity !<The interface points connectivity to get the interface for
    TYPE(InterfaceType), POINTER, INTENT(OUT) :: INTERFACE !<On return, the interface for the interface points connectivity. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfacePointsConnectivity_InterfaceGet",err,error,*998)

    IF(ASSOCIATED(INTERFACE)) CALL FlagError("Interface is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfacePointsConnectivity)) &
      & CALL FlagError("Interface points connectivity is not associated.",err,error,*999)

    INTERFACE=>interfacePointsConnectivity%INTERFACE
    IF(.NOT.ASSOCIATED(INTERFACE)) &
      & CALL FlagError("The interface points connectivity interface is not associated.",err,error,*999)
    
    EXITS("InterfacePointsConnectivity_InterfaceGet")
    RETURN
999 NULLIFY(interface)
998 ERRORSEXITS("InterfacePointsConnectivity_InterfaceGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfacePointsConnectivity_InterfaceGet

  !
  !=================================================================================================================================
  !

  !>Returns the interface mesh for an interface points connectivity
  SUBROUTINE InterfacePointsConnectivity_InterfaceMeshGet(interfacePointsConnectivity,interfaceMesh,err,error,*)

    !Argument Variables
    TYPE(InterfacePointsConnectivityType), POINTER, INTENT(INOUT) :: interfacePointsConnectivity !<The interface points connectivity to get the interface mesh for
    TYPE(MeshType), POINTER, INTENT(OUT) :: interfaceMesh !<On return, the interface mesh for the interface points connectivity. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfacePointsConnectivity_InterfaceMeshGet",err,error,*998)

    IF(ASSOCIATED(interfaceMesh)) CALL FlagError("Interface mesh is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfacePointsConnectivity)) &
      & CALL FlagError("Interface points connectivity is not associated.",err,error,*999)

    interfaceMesh=>interfacePointsConnectivity%interfaceMesh
    IF(.NOT.ASSOCIATED(interfaceMesh)) &
      & CALL FlagError("The interface points connectivity interface mesh is not associated.",err,error,*999)
    
    EXITS("InterfacePointsConnectivity_InterfaceMeshGet")
    RETURN
999 NULLIFY(interfaceMesh)
998 ERRORS("InterfacePointsConnectivity_InterfaceMeshGet",err,error)
    EXITS("InterfacePointsConnectivity_InterfaceMeshGet")
    RETURN 1
    
  END SUBROUTINE InterfacePointsConnectivity_InterfaceMeshGet

  !
  !================================================================================================================================
  !

END MODULE InterfaceAccessRoutines
