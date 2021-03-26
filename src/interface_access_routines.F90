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

  PUBLIC Interface_AssertIsFinished,Interface_AssertNotFinished

  PUBLIC Interface_CoordinateSystemGet

  PUBLIC Interface_CoupledMeshGet

  PUBLIC Interface_DataPointsGet

  PUBLIC Interface_FieldGet

  PUBLIC Interface_FieldsGet

  PUBLIC Interface_InterfaceConditionGet

  PUBLIC Interface_InterfacesGet

  PUBLIC Interface_MeshGet

  PUBLIC Interface_MeshConnectivityGet

  PUBLIC Interface_MeshesGet

  PUBLIC Interface_NodesGet

  PUBLIC Interface_NumberOfCoupledMeshesGet

  PUBLIC Interface_ParentRegionGet

  PUBLIC Interface_PointsConnectivityGet

  PUBLIC Interface_UserNumberFind

  PUBLIC Interface_UserNumberGet

  PUBLIC InterfaceMeshConnectivity_AssertIsFinished,InterfaceMeshConnectivity_AssertNotFinished

  PUBLIC InterfaceMeshConnectivity_BasisGet

  PUBLIC InterfaceMeshConnectivity_CoupledElementNumberGet

  PUBLIC InterfaceMeshConnectivity_CoupledNodeNumberGet

  PUBLIC InterfaceMeshConnectivity_InterfaceGet

  PUBLIC InterfaceMeshConnectivity_InterfaceMeshGet

  PUBLIC InterfaceMeshConnectivity_NumberOfCoupledMeshesGet
  
  PUBLIC InterfaceMeshConnectivity_NumberOfInterfaceElementsGet

  PUBLIC InterfaceMeshConnectivity_NumberOfInterfaceNodesGet

  PUBLIC InterfacePointsConnectivity_AssertIsFinished,InterfacePointsConnectivity_AssertNotFinished

  PUBLIC InterfacePointsConnectivity_DataPointsGet

  PUBLIC InterfacePointsConnectivity_InterfaceGet

  PUBLIC InterfacePointsConnectivity_InterfaceMeshGet

  PUBLIC InterfacePointsConnectivity_MaximumCoupledElementsGet

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

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(INTERFACE)) CALL FlagError("Interface is not associated.",err,error,*999)
#endif    

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

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(INTERFACE)) CALL FlagError("Interface is not associated.",err,error,*999)
#endif    

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
  
  !>Returns the coordinate system of an interface. \see OpenCMISS::Iron::cmfe_Interface_CoordinateSystemGet
  SUBROUTINE Interface_CoordinateSystemGet(interface,coordinateSystem,err,error,*)

    !Argument variables
    TYPE(InterfaceType), POINTER :: interface !<A pointer to the interface to get the coordinate system for
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem !<On exit, the coordinate system for the specified interface. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("Interface_CoordinateSystemGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(INTERFACE)) CALL FlagError("Interface is not associated.",err,error,*999)
    IF(.NOT.INTERFACE%interfaceFinished) CALL FlagError("Interface has not been finished.",err,error,*999)
    IF(ASSOCIATED(coordinateSystem)) CALL FlagError("Coordinate system is already associated.",err,error,*999)
#endif    

    coordinateSystem=>interface%coordinateSystem

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(coordinateSystem)) THEN
      localError="Coordinate system is not associated for interface number "// &
        & TRIM(NumberToVString(INTERFACE%userNumber,"*",err,error))
      IF(ASSOCIATED(INTERFACE%parentRegion)) localError=localError// &
        & " of parent region number "//TRIM(NumberToVString(INTERFACE%parentRegion%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif   
    
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
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("Interface_CoupledMeshGet",err,error,*998)

#ifdef WITH_PRECHECKS    
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
#endif
    
    coupledMesh=>interface%coupledMeshes(coupledMeshIndex)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(coupledMesh)) THEN
      localError="The coupled mesh for coupled mesh index "//TRIM(NumberToVString(coupledMeshIndex,"*",err,error))// &
        & " is not associated for interface number "//TRIM(NumberToVString(INTERFACE%userNumber,"*",err,error))
      IF(ASSOCIATED(INTERFACE%parentRegion)) localError=localError// &
        & " of parent region number "//TRIM(NumberToVString(INTERFACE%parentRegion%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("Interface_CoupledMeshGet")
    RETURN
999 NULLIFY(coupledMesh)
998 ERRORSEXITS("Interface_CoupledMeshGet",err,error)
    RETURN 1
    
  END SUBROUTINE Interface_CoupledMeshGet

  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the data points for a given user number in an interface. \see OpenCMISS::Iron::cmfe_Interface_DataPointsGet
  SUBROUTINE Interface_DataPointsGet(interface,userNumber,dataPoints,err,error,*)

    !Argument variables
    TYPE(InterfaceType), POINTER :: interface !<A pointer to the interface to get the data points for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the data points to get.
    TYPE(DataPointsType), POINTER :: dataPoints !<On exit, a pointer to the data points for the interface. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("Interface_DataPointsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dataPoints)) CALL FlagError("Data points is already associated.",err,error,*998)
    CALL Interface_AssertIsFinished(INTERFACE,err,error,*999)
    IF(.NOT.ASSOCIATED(INTERFACE%dataPointSets)) CALL FlagError("Interface data point sets is not associated.",err,error,*999)
#endif    

    NULLIFY(dataPoints)
    CALL DataPointSets_UserNumberFind(INTERFACE%dataPointSets,userNumber,dataPoints,err,error,*999)
#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(dataPoints)) THEN
      localError="Data points with a user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " do not exist on interface number "//TRIM(NumberToVString(interface%userNumber,"*",err,error))
      IF(ASSOCIATED(INTERFACE%parentRegion)) localError=localError// &
        & " of parent region number "//TRIM(NumberToVString(INTERFACE%parentRegion%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("Interface_DataPointsGet")
    RETURN
999 NULLIFY(dataPoints)
998 ERRORSEXITS("Interface_DataPointsGet",err,error)
    RETURN 1
    
  END SUBROUTINE Interface_DataPointsGet    

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the fields for an interface.
  SUBROUTINE Interface_FieldsGet(interface,fields,err,error,*)

    !Argument variables
    TYPE(InterfaceType), POINTER :: INTERFACE !<A pointer to the interface to get the fields for
    TYPE(FieldsType), POINTER :: fields !<On exit, a pointer to the fields for the interface. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("Interface_FieldsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(fields)) CALL FlagError("Fields is already associated.",err,error,*998)
    CALL Interface_AssertIsFinished(INTERFACE,err,error,*999)
#endif    

    fields=>interface%fields
    
#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(fields)) THEN
      localError="The interface fields are not associated on interface number "// &
        & TRIM(NumberToVString(interface%userNumber,"*",err,error))
      IF(ASSOCIATED(INTERFACE%parentRegion)) localError=localError// &
        & " of parent region number "//TRIM(NumberToVString(INTERFACE%parentRegion%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
    EXITS("Interface_FieldsGet")
    RETURN
999 NULLIFY(fields)
998 ERRORSEXITS("Interface_FieldsGet",err,error)
    RETURN 1
    
  END SUBROUTINE Interface_FieldsGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to a field for a given user number in an interface. \see OpenCMISS::Iron::cmfe_Interface_FieldGet
  SUBROUTINE Interface_FieldGet(interface,userNumber,field,err,error,*)

    !Argument variables
    TYPE(InterfaceType), POINTER :: INTERFACE !<A pointer to the interface to get the field for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the field to get.
    TYPE(FieldType), POINTER :: field !<On exit, a pointer to the field for the interface. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("Interface_FieldGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(field)) CALL FlagError("Field is already associated.",err,error,*998)
    CALL Interface_AssertIsFinished(INTERFACE,err,error,*999)
    IF(.NOT.ASSOCIATED(INTERFACE%fields)) CALL FlagError("Interface fields is not associated.",err,error,*998)
#endif    

    NULLIFY(field)
    CALL Field_UserNumberFind(userNumber,INTERFACE,field,err,error,*999)
#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(field)) THEN
      localError="A field with a user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " do not exist on interface number "//TRIM(NumberToVString(interface%userNumber,"*",err,error))
      IF(ASSOCIATED(INTERFACE%parentRegion)) localError=localError// &
        & " of parent region number "//TRIM(NumberToVString(INTERFACE%parentRegion%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
    EXITS("Interface_FieldGet")
    RETURN
999 NULLIFY(field)
998 ERRORSEXITS("Interface_FieldGet",err,error)
    RETURN 1
    
  END SUBROUTINE Interface_FieldGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to a interface condition for a given user number in an interface. \see OpenCMISS::Iron::cmfe_Interface_InterfaceConditionGet
  SUBROUTINE Interface_InterfaceConditionGet(interface,userNumber,interfaceCondition,err,error,*)

    !Argument variables
    TYPE(InterfaceType), POINTER :: INTERFACE !<A pointer to the interface to get the interface condition for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the interface condition to get.
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<On exit, a pointer to the interface condition for the interface. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("Interface_InterfaceConditionGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is already associated.",err,error,*998)
    CALL Interface_AssertIsFinished(INTERFACE,err,error,*999)
    IF(.NOT.ASSOCIATED(interface%interfaceConditions)) &
      & CALL FlagError("Interface interface conditions is not associated.",err,error,*998)
#endif    

    NULLIFY(interfaceCondition)
    CALL InterfaceCondition_UserNumberFind(userNumber,INTERFACE,interfaceCondition,err,error,*999)
#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(interfaceCondition)) THEN
      localError="An interface Condition with a user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " do not exist on interface number "//TRIM(NumberToVString(interface%userNumber,"*",err,error))
      IF(ASSOCIATED(INTERFACE%parentRegion)) localError=localError// &
        & " of parent region number "//TRIM(NumberToVString(INTERFACE%parentRegion%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
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
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("Interface_InterfacesGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(interfaces)) CALL FlagError("Interfaces is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(INTERFACE)) CALL FlagError("Interface is not associated.",err,error,*999)
#endif    

    interfaces=>interface%interfaces

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(interfaces)) THEN
      localError="Interfaces is not associated for interface number "// &
        & TRIM(NumberToVString(interface%userNumber,"*",err,error))
      IF(ASSOCIATED(INTERFACE%parentRegion)) localError=localError// &
        & " of parent region number "//TRIM(NumberToVString(INTERFACE%parentRegion%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    EXITS("Interface_InterfacesGet")
    RETURN
999 NULLIFY(interfaces)
998 ERRORSEXITS("Interface_InterfacesGet",err,error)
    RETURN 1
    
  END SUBROUTINE Interface_InterfacesGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the mesh for a given user number in a interface. \see OpenCMISS::Iron::cmfe_Interface_MeshGet
  SUBROUTINE Interface_MeshGet(interface,userNumber,mesh,err,error,*)

    !Argument variables
    TYPE(InterfaceType), POINTER :: interface !<A pointer to the interface to get the mesh for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the mesh to get.
    TYPE(MeshType), POINTER :: mesh !<On exit, a pointer to the mesh for the interface. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("Interface_MeshGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(mesh)) CALL FlagError("Mesh is already associated.",err,error,*998)
    CALL Interface_AssertIsFinished(INTERFACE,err,error,*999)
#endif    
    
    NULLIFY(mesh)
    CALL Mesh_UserNumberFind(userNumber,INTERFACE,mesh,err,error,*999)
#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(mesh)) THEN
      localError="A mesh with a user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " does not exist on interface number "//TRIM(NumberToVString(interface%userNumber,"*",err,error))
      IF(ASSOCIATED(INTERFACE%parentRegion)) localError=localError// &
        & " of parent region number "//TRIM(NumberToVString(INTERFACE%parentRegion%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
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
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("Interface_MeshConnectivityGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(meshConnectivity)) CALL FlagError("Mesh connectivity is already associated.",err,error,*998)
    CALL Interface_AssertIsFinished(INTERFACE,err,error,*999)
#endif    
    
    meshConnectivity=>interface%meshConnectivity

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(meshConnectivity)) THEN
      localError="Mesh connectivity is not associated on intereface number "// &
        & TRIM(NumberToVString(interface%userNumber,"*",err,error))
      IF(ASSOCIATED(INTERFACE%parentRegion)) localError=localError// &
        & " of parent region number "//TRIM(NumberToVString(INTERFACE%parentRegion%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
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
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("Interface_MeshesGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(meshes)) CALL FlagError("Meshes is already associated.",err,error,*998)
    CALL Interface_AssertIsFinished(INTERFACE,err,error,*999)
#endif    
 
    meshes=>interface%meshes

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(meshes)) THEN
      localError="Meshes is not associated for interface number "//TRIM(NumberToVString(INTERFACE%userNumber,"*",err,error))
      IF(ASSOCIATED(INTERFACE%parentRegion)) localError=localError// &
        & " of parent region number "//TRIM(NumberToVString(INTERFACE%parentRegion%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
    EXITS("Interface_MeshesGet")
    RETURN
999 NULLIFY(meshes)
998 ERRORSEXITS("Interface_MeshesGet",err,error)
    RETURN 1
    
  END SUBROUTINE Interface_MeshesGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the nodes for a interface. \see OpenCMISS::Iron::cmfe_Interface_NodesGet
  SUBROUTINE Interface_NodesGet(interface,nodes,err,error,*)

    !Argument variables
    TYPE(InterfaceType), POINTER :: interface !<A pointer to the interface to get the nodes for
    TYPE(NodesType), POINTER :: nodes !<On exit, a pointer to the nodes for the interface. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("Interface_NodesGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(nodes)) CALL FlagError("Nodes is already associated.",err,error,*998)
    CALL Interface_AssertIsFinished(INTERFACE,err,error,*999)
#endif    

    nodes=>interface%nodes

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(nodes)) THEN
      localError="Nodes is not associated for interface number "//TRIM(NumberToVString(INTERFACE%userNumber,"*",err,error))
      IF(ASSOCIATED(INTERFACE%parentRegion)) localError=localError// &
        & " of parent region number "//TRIM(NumberToVString(INTERFACE%parentRegion%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
    EXITS("Interface_NodesGet")
    RETURN
999 NULLIFY(nodes)
998 ERRORSEXITS("Interface_NodesGet",err,error)
    RETURN 1
    
  END SUBROUTINE Interface_NodesGet

  !
  !================================================================================================================================
  !

  !Gets the number of coupled meshes in a interface.
  SUBROUTINE Interface_NumberOfCoupledMeshesGet(INTERFACE,numberOfCoupledMeshes,err,error,*)

    !Argument variables
    TYPE(InterfaceType), POINTER :: interface !<A pointer to the interface to get the number of coupled meshes for
    INTEGER(INTG), INTENT(OUT) :: numberOfCoupledMeshes !<On exit, the number of coupled meshes for the interface.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("Interface_NumberOfCoupledMeshesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    CALL Interface_AssertIsFinished(INTERFACE,err,error,*999)
#endif    

    numberOfCoupledMeshes=INTERFACE%numberOfCoupledMeshes
      
    EXITS("Interface_NumberOfCoupledMeshesGet")
    RETURN
999 ERRORSEXITS("Interface_NumberOfCoupledMeshesGet",err,error)
    RETURN 1
    
  END SUBROUTINE Interface_NumberOfCoupledMeshesGet

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
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("Interface_ParentRegionGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(parentRegion)) CALL FlagError("Parent region is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(INTERFACE)) CALL FlagError("Interface is not associated.",err,error,*999)
#endif    
 
    parentRegion=>interface%parentRegion

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(parentRegion)) THEN
      localError="The parent region for interface number "//TRIM(NumberToVString(interface%userNumber,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
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
#ifdef WITH_POSTCHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("Interface_PointConnectivityGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(pointsConnectivity)) CALL FlagError("Point connectivity is already associated.",err,error,*998)
#endif    
    CALL Interface_AssertIsFinished(INTERFACE,err,error,*999)
     
    pointsConnectivity=>interface%pointsConnectivity

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(pointsConnectivity)) THEN
      localError="Points connectivity is not associated on intereface number "// &
        & TRIM(NumberToVString(interface%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
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
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("Interface_UserNumberFind",err,error,*999)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(ASSOCIATED(INTERFACE)) CALL FlagError("Interface is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(parentRegion)) CALL FlagError("Parent region is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(parentRegion%interfaces)) THEN
      localError="The interfaces on parent region number "// &
        & TRIM(NumberToVString(parentRegion%userNumber,"*",err,error))//" are not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    !Get the interface from the user number
    NULLIFY(INTERFACE)
    IF(ALLOCATED(parentRegion%interfaces%interfaces)) THEN
      DO interfaceIdx=1,parentRegion%interfaces%numberOfInterfaces
#ifdef WITH_PRECHECKS        
        IF(.NOT.ASSOCIATED(parentRegion%interfaces%interfaces(interfaceIdx)%ptr)) THEN
          localError="The interface pointer in interfaces is not associated for interface index "// &
            & TRIM(NumberToVString(interfaceIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
#endif        
        IF(parentRegion%interfaces%interfaces(interfaceIdx)%ptr%userNumber==userNumber) THEN
          INTERFACE=>parentRegion%interfaces%interfaces(interfaceIdx)%ptr
          EXIT
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

  !>Returns the user number for a given interface.
  SUBROUTINE Interface_UserNumberGet(INTERFACE,userNumber,err,error,*)

    !Argument variables
    TYPE(InterfaceType), POINTER :: interface !<A pointer to the interface to get the user number for
    INTEGER(INTG), INTENT(OUT) :: userNumber !<On exit, the user number for the interface.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Interface_UserNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(interface)) CALL FlagError("Interface is not associated.",err,error,*999)
#endif    
     
    userNumber=interface%userNumber
    
    EXITS("Interface_UserNumberGet")
    RETURN
999 ERRORSEXITS("Interface_UserNumberGet",err,error)
    RETURN 1
    
  END SUBROUTINE Interface_UserNumberGet

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

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(interfaceMeshConnectivity)) CALL FlagError("Interface mesh connectivity is not associated.",err,error,*999)
#endif    

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

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(interfaceMeshConnectivity)) CALL FlagError("Interface mesh connectivity is not associated.",err,error,*999)
#endif    

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

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(basis)) CALL FlagError("Basis is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMeshConnectivity)) &
      & CALL FlagError("Interface mesh connectivity is not associated.",err,error,*999)
#endif    

    basis=>interfaceMeshConnectivity%basis

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(basis)) &
      & CALL FlagError("The interface mesh connectivity basis is not associated.",err,error,*999)
#endif    
    
    EXITS("InterfaceMeshConnectivity_BasisGet")
    RETURN
999 NULLIFY(basis)
998 ERRORSEXITS("InterfaceMeshConnectivity_BasisGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMeshConnectivity_BasisGet

  !
  !=================================================================================================================================
  !

  !>Returns the coupled element number for an interface element number and a couple mesh index for an interface mesh connectivity
  SUBROUTINE InterfaceMeshConnectivity_CoupledElementNumberGet(interfaceMeshConnectivity,interfaceElementNumber, &
    & meshIndex,coupledElementNumber,err,error,*)

    !Argument Variables
    TYPE(InterfaceMeshConnectivityType), POINTER, INTENT(INOUT) :: interfaceMeshConnectivity !<The interface mesh connectivity to get the coupled element for
    INTEGER(INTG), INTENT(IN) :: interfaceElementNumber !<The interface element number to get the coupled element number for
    INTEGER(INTG), INTENT(IN) :: meshIndex !<The mesh index to get the coupled element number for.
    INTEGER(INTG), INTENT(OUT) :: coupledElementNumber !<On exit, the coupled element number in the mesh.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("InterfaceMeshConnectivity_CoupledElementNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(interfaceMeshConnectivity)) &
      & CALL FlagError("Interface mesh connectivity is not associated.",err,error,*999)
    IF(interfaceElementNumber<1.OR.interfaceElementNumber>interfaceMeshConnectivity%numberOfInterfaceElements) THEN
      localError="The specified interface element number of "//TRIM(NumberToVString(interfaceElementNumber,"*",err,error))// &
        & " is invalid. The interface element number should be >= 1 and <= "// &
        & TRIM(NumberToVString(interfaceMeshConnectivity%numberOfInterfaceElements,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(meshIndex<1.OR.meshIndex>interfaceMeshConnectivity%numberOfCoupledMeshes) THEN
      localError="The specified mesh index of "//TRIM(NumberToVString(meshIndex,"*",err,error))// &
        & " is invalid. The mesh index should be >= 1 and <= "// &
        & TRIM(NumberToVString(interfaceMeshConnectivity%numberOfCoupledMeshes,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(interfaceMeshConnectivity%elementConnectivity)) &
      & CALL FlagError("The element connectivity array is not allocated for the interface mesh connectivity.",err,error,*999)
#endif    

    coupledElementNumber=interfaceMeshConnectivity%elementConnectivity(interfaceElementNumber,meshIndex)%coupledElementNumber
    
    EXITS("InterfaceMeshConnectivity_CoupledElementNumberGet")
    RETURN
999 ERRORSEXITS("InterfaceMeshConnectivity_CoupledElementNumberGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMeshConnectivity_CoupledElementNumberGet

  !
  !=================================================================================================================================
  !

  !>Returns the coupled node number for an interface node number and a couple mesh index for an interface mesh connectivity
  SUBROUTINE InterfaceMeshConnectivity_CoupledNodeNumberGet(interfaceMeshConnectivity,interfaceNodeNumber, &
    & meshIndex,coupledNodeNumber,err,error,*)

    !Argument Variables
    TYPE(InterfaceMeshConnectivityType), POINTER, INTENT(INOUT) :: interfaceMeshConnectivity !<The interface mesh connectivity to get the coupled node for
    INTEGER(INTG), INTENT(IN) :: interfaceNodeNumber !<The interface node number to get the coupled node number for
    INTEGER(INTG), INTENT(IN) :: meshIndex !<The mesh index to get the coupled node number for.
    INTEGER(INTG), INTENT(OUT) :: coupledNodeNumber !<On exit, the coupled node number in the mesh.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("InterfaceMeshConnectivity_CoupledNodeNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(interfaceMeshConnectivity)) &
      & CALL FlagError("Interface mesh connectivity is not associated.",err,error,*999)
    IF(interfaceNodeNumber<1.OR.interfaceNodeNumber>interfaceMeshConnectivity%numberOfInterfaceNodes) THEN
      localError="The specified interface node number of "//TRIM(NumberToVString(interfaceNodeNumber,"*",err,error))// &
        & " is invalid. The interface node number should be >= 1 and <= "// &
        & TRIM(NumberToVString(interfaceMeshConnectivity%numberOfInterfaceNodes,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(meshIndex<1.OR.meshIndex>interfaceMeshConnectivity%numberOfCoupledMeshes) THEN
      localError="The specified mesh index of "//TRIM(NumberToVString(meshIndex,"*",err,error))// &
        & " is invalid. The mesh index should be >= 1 and <= "// &
        & TRIM(NumberToVString(interfaceMeshConnectivity%numberOfCoupledMeshes,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(interfaceMeshConnectivity%coupledNodes)) &
      & CALL FlagError("The coupled nodes array is not allocated for the interface mesh connectivity.",err,error,*999)
#endif    

    coupledNodeNumber=interfaceMeshConnectivity%coupledNodes(meshIndex,interfaceNodeNumber)
    
    EXITS("InterfaceMeshConnectivity_CoupledNodeNumberGet")
    RETURN
999 ERRORS("InterfaceMeshConnectivity_CoupledNodeNumberGet",err,error)
    EXITS("InterfaceMeshConnectivity_CoupledNodeNumberGet")
    RETURN 1
    
  END SUBROUTINE InterfaceMeshConnectivity_CoupledNodeNumberGet

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

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(INTERFACE)) CALL FlagError("Interface is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMeshConnectivity)) &
      & CALL FlagError("Interface mesh connectivity is not associated.",err,error,*999)
#endif    

    INTERFACE=>interfaceMeshConnectivity%INTERFACE

#ifdef WITH_POSTCHECKS      
    IF(.NOT.ASSOCIATED(INTERFACE)) &
      & CALL FlagError("The interface mesh connectivity interface is not associated.",err,error,*999)
#endif    
    
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

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(interfaceMesh)) CALL FlagError("Interface mesh is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMeshConnectivity)) &
      & CALL FlagError("Interface mesh connectivity is not associated.",err,error,*999)
#endif    

    interfaceMesh=>interfaceMeshConnectivity%interfaceMesh

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(interfaceMesh)) &
      & CALL FlagError("The interface mesh connectivity interface mesh is not associated.",err,error,*999)
#endif    
    
    EXITS("InterfaceMeshConnectivity_InterfaceMeshGet")
    RETURN
999 NULLIFY(interfaceMesh)
998 ERRORSEXITS("InterfaceMeshConnectivity_InterfaceMeshGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMeshConnectivity_InterfaceMeshGet

  !
  !=================================================================================================================================
  !

  !>Returns the number of coupled meshes for an interface mesh connectivity
  SUBROUTINE InterfaceMeshConnectivity_NumberOfCoupledMeshesGet(interfaceMeshConnectivity,numberOfCoupledMeshes,err,error,*)

    !Argument Variables
    TYPE(InterfaceMeshConnectivityType), POINTER, INTENT(INOUT) :: interfaceMeshConnectivity !<The interface mesh connectivity to get the number of coupled meshes for
    INTEGER(INTG), INTENT(OUT) :: numberOfCoupledMeshes !<On return, the number of coupled meshes for the interface mesh connectivity.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMeshConnectivity_NumberOfCoupledMeshesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(interfaceMeshConnectivity)) &
      & CALL FlagError("Interface mesh connectivity is not associated.",err,error,*999)
#endif    

    numberOfCoupledMeshes=interfaceMeshConnectivity%numberOfCoupledMeshes

    EXITS("InterfaceMeshConnectivity_NumberOfCoupledMeshesGet")
    RETURN
999 ERRORS("InterfaceMeshConnectivity_NumberOfCoupledMeshesGet",err,error)
    EXITS("InterfaceMeshConnectivity_NumberOfCoupledMeshesGet")
    RETURN 1
    
  END SUBROUTINE InterfaceMeshConnectivity_NumberOfCoupledMeshesGet

  !
  !=================================================================================================================================
  !

  !>Returns the number of interface elements for an interface mesh connectivity
  SUBROUTINE InterfaceMeshConnectivity_NumberOfInterfaceElementsGet(interfaceMeshConnectivity,numberOfInterfaceElements,err,error,*)

    !Argument Variables
    TYPE(InterfaceMeshConnectivityType), POINTER, INTENT(INOUT) :: interfaceMeshConnectivity !<The interface mesh connectivity to get the number of interface elements for
    INTEGER(INTG), INTENT(OUT) :: numberOfInterfaceElements !<On return, the number of interface elements for the interface mesh connectivity.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMeshConnectivity_NumberOfInterfaceElementsGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(interfaceMeshConnectivity)) &
      & CALL FlagError("Interface mesh connectivity is not associated.",err,error,*999)
#endif    

    numberOfInterfaceElements=interfaceMeshConnectivity%numberOfInterfaceElements

    EXITS("InterfaceMeshConnectivity_NumberOfInterfaceElementsGet")
    RETURN
999 ERRORS("InterfaceMeshConnectivity_NumberOfInterfaceElementsGet",err,error)
    EXITS("InterfaceMeshConnectivity_NumberOfInterfaceElementsGet")
    RETURN 1
    
  END SUBROUTINE InterfaceMeshConnectivity_NumberOfInterfaceElementsGet

  !
  !=================================================================================================================================
  !

  !>Returns the number of interface nodes for an interface mesh connectivity
  SUBROUTINE InterfaceMeshConnectivity_NumberOfInterfaceNodesGet(interfaceMeshConnectivity,numberOfInterfaceNodes,err,error,*)

    !Argument Variables
    TYPE(InterfaceMeshConnectivityType), POINTER, INTENT(INOUT) :: interfaceMeshConnectivity !<The interface mesh connectivity to get the number of interface nodes for
    INTEGER(INTG), INTENT(OUT) :: numberOfInterfaceNodes !<On return, the number of interface nodes for the interface mesh connectivity.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMeshConnectivity_NumberOfInterfaceNodesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(interfaceMeshConnectivity)) &
      & CALL FlagError("Interface mesh connectivity is not associated.",err,error,*999)
#endif    

    numberOfInterfaceNodes=interfaceMeshConnectivity%numberOfInterfaceNodes

    EXITS("InterfaceMeshConnectivity_NumberOfInterfaceNodesGet")
    RETURN
999 ERRORS("InterfaceMeshConnectivity_NumberOfInterfaceNodesGet",err,error)
    EXITS("InterfaceMeshConnectivity_NumberOfInterfaceNodesGet")
    RETURN 1
    
  END SUBROUTINE InterfaceMeshConnectivity_NumberOfInterfaceNodesGet

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

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(interfacePointsConnectivity)) &
      & CALL FlagError("Interface points connectivity is not associated.",err,error,*999)
#endif    

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

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(interfacePointsConnectivity)) &
      & CALL FlagError("Interface points connectivity is not associated.",err,error,*999)
#endif    

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

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dataPoints)) CALL FlagError("Data points is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfacePointsConnectivity)) &
      & CALL FlagError("Interface points connectivity is not associated.",err,error,*999)
#endif    

    dataPoints=>interfacePointsConnectivity%dataPoints

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(dataPoints)) &
      & CALL FlagError("The interface points connectivity data points is not associated.",err,error,*999)
#endif    
    
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

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(INTERFACE)) CALL FlagError("Interface is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfacePointsConnectivity)) &
      & CALL FlagError("Interface points connectivity is not associated.",err,error,*999)
#endif    

    INTERFACE=>interfacePointsConnectivity%INTERFACE

#ifdef WITH_POSTCHECKS      
    IF(.NOT.ASSOCIATED(INTERFACE)) &
      & CALL FlagError("The interface points connectivity interface is not associated.",err,error,*999)
#endif    
    
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

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(interfaceMesh)) CALL FlagError("Interface mesh is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfacePointsConnectivity)) &
      & CALL FlagError("Interface points connectivity is not associated.",err,error,*999)
#endif    

    interfaceMesh=>interfacePointsConnectivity%interfaceMesh

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(interfaceMesh)) &
      & CALL FlagError("The interface points connectivity interface mesh is not associated.",err,error,*999)
#endif    
    
    EXITS("InterfacePointsConnectivity_InterfaceMeshGet")
    RETURN
999 NULLIFY(interfaceMesh)
998 ERRORS("InterfacePointsConnectivity_InterfaceMeshGet",err,error)
    EXITS("InterfacePointsConnectivity_InterfaceMeshGet")
    RETURN 1
    
  END SUBROUTINE InterfacePointsConnectivity_InterfaceMeshGet

  !
  !=================================================================================================================================
  !

  !>Returns the maximum number of coupled elements for an interface points connectivity
  SUBROUTINE InterfacePointsConnectivity_MaximumCoupledElementsGet(interfacePointsConnectivity,coupledMeshIdx, &
    & maximumCoupledElements,err,error,*)

    !Argument Variables
    TYPE(InterfacePointsConnectivityType), POINTER, INTENT(INOUT) :: interfacePointsConnectivity !<The interface points connectivity to get the maximum number of coupled elements for
    INTEGER(INTG), INTENT(IN) :: coupledMeshIdx !<The coupled mesh index to get the maximum number of coupled elements for
    INTEGER(INTG), INTENT(OUT) :: maximumCoupledElements !<On return, the maximum number of coupled elements for the specified coupled mesh index in the interface points connectivity.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif
 
    ENTERS("InterfacePointsConnectivity_MaximumCoupledElementsGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(interfacePointsConnectivity)) &
      & CALL FlagError("Interface points connectivity is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(interfacePointsConnectivity%maxNumberOfCoupledElements)) &
      & CALL FlagError("The maximum number of coupled elements is not allocated for the interface points connectivity", &
      & err,error,*999)
    IF(coupledMeshIdx<1.OR.coupledMeshIdx>SIZE(interfacePointsConnectivity%maxNumberOfCoupledElements,1)) THEN
      localError="The specified coupled mesh index of "//TRIM(NumberToVString(coupledMeshIdx,"*",err,error))// &
        & " is invalid. The coupled mesh index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(interfacePointsConnectivity%maxNumberOfCoupledElements,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    maximumCoupledElements=interfacePointsConnectivity%maxNumberOfCoupledElements(coupledMeshIdx)
    
    EXITS("InterfacePointsConnectivity_MaximumCoupledElementsGet")
    RETURN
999 ERRORS("InterfacePointsConnectivity_MaximumCoupledElementsGet",err,error)
    EXITS("InterfacePointsConnectivity_MaximumCoupledElementsGet")
    RETURN 1
    
  END SUBROUTINE InterfacePointsConnectivity_MaximumCoupledElementsGet

  !
  !================================================================================================================================
  !

END MODULE InterfaceAccessRoutines
