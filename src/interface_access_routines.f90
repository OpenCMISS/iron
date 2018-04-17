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

  PUBLIC Interface_CoordinateSystemGet

  PUBLIC INTERFACE_COORDINATE_SYSTEM_GET

  PUBLIC Interface_DataPointsGet

  PUBLIC INTERFACE_DATA_POINTS_GET

  PUBLIC Interface_FieldGet

  PUBLIC Interface_InterfaceConditionGet

  PUBLIC Interface_MeshGet

  PUBLIC Interface_MeshConnectivityGet

  PUBLIC Interface_NodesGet

  PUBLIC INTERFACE_NODES_GET

  PUBLIC Interface_PointsConnectivityGet

  PUBLIC Interface_UserNumberFind

  PUBLIC INTERFACE_USER_NUMBER_FIND

CONTAINS

  !
  !================================================================================================================================
  !
  
  !>Returns the coordinate system of an interface. \see OPENCMISS::Iron::cmfe_Interface_CoordinateSystemGet
  SUBROUTINE Interface_CoordinateSystemGet(interface,coordinateSystem,err,error,*)

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: interface !<A pointer to the interface to get the coordinate system for
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: coordinateSystem !<On exit, the coordinate system for the specified interface. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Interface_CoordinateSystemGet",err,error,*999)

    IF(.NOT.ASSOCIATED(INTERFACE)) CALL FlagError("Interface is not associated.",err,error,*999)
    IF(.NOT.INTERFACE%INTERFACE_FINISHED) CALL FlagError("Interface has not been finished.",err,error,*999)
    IF(ASSOCIATED(coordinateSystem)) CALL FlagError("Coordinate system is already associated.",ERR,ERROR,*999)

    coordinateSystem=>interface%COORDINATE_SYSTEM
    IF(.NOT.ASSOCIATED(coordinateSystem)) THEN
      localError="The coordinate system for interface number "//TRIM(NumberToVString(interface%USER_NUMBER,"*",err,error))// &
        & " is not associated."
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
  
  !>Returns a pointer to the data points for a given user number in an interface. \see OPENCMISS::Iron::cmfe_Interface_DataPointsGet
  SUBROUTINE Interface_DataPointsGet(interface,userNumber,dataPoints,err,error,*)

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: interface !<A pointer to the interface to get the data points for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the data points to get.
    TYPE(DataPointsType), POINTER :: dataPoints !<On exit, a pointer to the data points for the interface. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Interface_DataPointsGet",err,error,*998)

    IF(.NOT.ASSOCIATED(interface)) CALL FlagError("Interface is not associated.",err,error,*998)
    IF(.NOT.interface%INTERFACE_FINISHED) CALL FlagError("Interface has not been finished.",err,error,*998)
    IF(ASSOCIATED(dataPoints)) CALL FlagError("Data points is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interface%dataPointSets)) CALL FlagError("Interface data point sets is not associated.",err,error,*998)

    NULLIFY(dataPoints)
    CALL DataPointSets_UserNumberFind(INTERFACE%dataPointSets,userNumber,dataPoints,err,error,*999)
    IF(.NOT.ASSOCIATED(dataPoints)) THEN
      localError="Data points with a user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " do not exist on interface number "//TRIM(NumberToVString(interface%USER_NUMBER,"*",err,error))//"."
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
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to get the field for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the field to get.
    TYPE(FIELD_TYPE), POINTER :: field !<On exit, a pointer to the field for the interface. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Interface_FieldGet",err,error,*998)

    IF(.NOT.ASSOCIATED(interface)) CALL FlagError("Interface is not associated.",err,error,*998)
    IF(.NOT.interface%INTERFACE_FINISHED) CALL FlagError("Interface has not been finished.",err,error,*998)
    IF(ASSOCIATED(field)) CALL FlagError("Field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interface%fields)) CALL FlagError("Interface fields is not associated.",err,error,*998)

    NULLIFY(field)
    CALL Field_UserNumberFind(userNumber,interface,field,err,error,*999)
    IF(.NOT.ASSOCIATED(field)) THEN
      localError="A field with a user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " do not exist on interface number "//TRIM(NumberToVString(interface%USER_NUMBER,"*",err,error))//"."
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
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to get the interface condition for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the interface condition to get.
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: interfaceCondition !<On exit, a pointer to the interface condition for the interface. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Interface_InterfaceConditionGet",err,error,*998)

    IF(.NOT.ASSOCIATED(interface)) CALL FlagError("Interface is not associated.",err,error,*998)
    IF(.NOT.interface%INTERFACE_FINISHED) CALL FlagError("Interface has not been finished.",err,error,*998)
    IF(ASSOCIATED(interfaceCondition)) CALL FlagError("InterfaceCondition is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(INTERFACE%INTERFACE_CONDITIONS)) &
      & CALL FlagError("Interface interface conditions is not associated.",err,error,*998)

    NULLIFY(interfaceCondition)
    CALL InterfaceCondition_UserNumberFind(userNumber,interface,interfaceCondition,err,error,*999)
    IF(.NOT.ASSOCIATED(interfaceCondition)) THEN
      localError="An interface Condition with a user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " do not exist on interface number "//TRIM(NumberToVString(interface%USER_NUMBER,"*",err,error))//"."
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

  !>Returns a pointer to the mesh for a given user number in a interface. \see OPENCMISS::Iron::cmfe_Interface_MeshGet
  SUBROUTINE Interface_MeshGet(interface,userNumber,mesh,err,error,*)

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: interface !<A pointer to the interface to get the mesh for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the mesh to get.
    TYPE(MESH_TYPE), POINTER :: mesh !<On exit, a pointer to the mesh for the interface. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Interface_MeshGet",err,error,*998)

    IF(.NOT.ASSOCIATED(interface)) CALL FlagError("Interface is not associated.",err,error,*998)
    IF(.NOT.interface%INTERFACE_FINISHED) CALL FlagError("Interface has not been finished.",err,error,*998)
    IF(ASSOCIATED(mesh)) CALL FlagError("Mesh is already associated.",err,error,*998)
    
    NULLIFY(mesh)
    CALL Mesh_UserNumberFind(userNumber,interface,mesh,err,error,*999)
    IF(.NOT.ASSOCIATED(mesh)) THEN
      localError="A mesh with a user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " does not exist on interface number "//TRIM(NumberToVString(interface%USER_NUMBER,"*",err,error))//"."
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
    TYPE(INTERFACE_TYPE), POINTER :: interface !<A pointer to the interface to get the mesh connectivity for
    TYPE(INTERFACE_MESH_CONNECTIVITY_TYPE), POINTER :: meshConnectivity !<On exit, a pointer to the mesh connectivity for the interface. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Interface_MeshConnectivityGet",err,error,*998)

    IF(ASSOCIATED(meshConnectivity)) CALL FlagError("Mesh connectivity is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interface)) CALL FlagError("Interface is not associated.",err,error,*998)
    IF(.NOT.interface%INTERFACE_FINISHED) CALL FlagError("Interface has not been finished.",err,error,*998)
    
    meshConnectivity=>interface%MESH_CONNECTIVITY
    IF(.NOT.ASSOCIATED(meshConnectivity)) THEN
      localError="Mesh connectivity is not associated on intereface number "// &
        & TRIM(NumberToVString(interface%USER_NUMBER,"*",err,error))//"."
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

  !>Returns a pointer to the nodes for a interface. \see OPENCMISS::Iron::cmfe_InterfaceNodesGet
  SUBROUTINE Interface_NodesGet(interface,nodes,err,error,*)

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: interface !<A pointer to the interface to get the nodes for
    TYPE(NodesType), POINTER :: nodes !<On exit, a pointer to the nodes for the interface. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Interface_NodesGet",err,error,*998)

    IF(.NOT.ASSOCIATED(interface)) CALL FlagError("Interface is not associated.",err,error,*998)
    IF(.NOT.interface%INTERFACE_FINISHED) CALL FlagError("Interface has not been finished.",err,error,*998)
    IF(ASSOCIATED(nodes)) CALL FlagError("Nodes is already associated.",err,error,*998)

    nodes=>interface%nodes
    IF(.NOT.ASSOCIATED(nodes)) CALL FlagError("Interface nodes is not associated.",err,error,*999)
       
    EXITS("Interface_NodesGet")
    RETURN
999 NULLIFY(nodes)
998 ERRORSEXITS("Interface_NodesGet",err,error)
    RETURN 1
    
  END SUBROUTINE Interface_NodesGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the points connectivity for a given interface.
  SUBROUTINE Interface_PointsConnectivityGet(INTERFACE,pointsConnectivity,err,error,*)

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: interface !<A pointer to the interface to get the points connectivity for
    TYPE(InterfacePointsConnectivityType), POINTER :: pointsConnectivity !<On exit, a pointer to the points connectivity for the interface. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Interface_PointConnectivityGet",err,error,*998)

    IF(ASSOCIATED(pointsConnectivity)) CALL FlagError("Point connectivity is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interface)) CALL FlagError("Interface is not associated.",err,error,*999)
    IF(.NOT.interface%INTERFACE_FINISHED) CALL FlagError("Interface has not been finished.",err,error,*999)
    
    pointsConnectivity=>interface%pointsConnectivity
    IF(.NOT.ASSOCIATED(pointsConnectivity)) THEN
      localError="Points connectivity is not associated on intereface number "// &
        & TRIM(NumberToVString(interface%USER_NUMBER,"*",err,error))//"."
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
    TYPE(REGION_TYPE), POINTER :: parentRegion !<The parent region to find the interface in    
    TYPE(INTERFACE_TYPE), POINTER :: interface !<On return a pointer to the interface with the given user number. If no interface with that user number exists then the pointer is returned as NULL. Must not be associated on entry.
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
        & TRIM(NumberToVString(parentRegion%USER_NUMBER,"*",err,error))//" are not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    !Get the interface from the user number
    NULLIFY(INTERFACE)
    IF(ASSOCIATED(parentRegion%interfaces%interfaces)) THEN
      DO interfaceIdx=1,parentRegion%interfaces%NUMBER_OF_INTERFACES
        IF(ASSOCIATED(parentRegion%interfaces%interfaces(interfaceIdx)%ptr)) THEN
          IF(parentRegion%interfaces%interfaces(interfaceIdx)%ptr%USER_NUMBER==userNumber) THEN
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
  !================================================================================================================================
  !

END MODULE InterfaceAccessRoutines
