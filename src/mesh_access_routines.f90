!> \file
!> \author Chris Bradley
!> \brief This module contains all mesh access method routines.
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

!> This module contains all mesh access method routines.
MODULE MeshAccessRoutines
  
  USE BaseRoutines
  USE DecompositionAccessRoutines
  USE Kinds
  USE ISO_VARYING_STRING
  USE Strings
  USE Types
  
#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  INTERFACE MESH_TOPOLOGY_ELEMENTS_GET
    MODULE PROCEDURE Mesh_MeshElementsGet
  END INTERFACE MESH_TOPOLOGY_ELEMENTS_GET

  INTERFACE MeshTopology_NodesGet
    MODULE PROCEDURE Mesh_MeshNodesGet
  END INTERFACE MeshTopology_NodesGet

  INTERFACE Mesh_UserNumberFind
    MODULE PROCEDURE Mesh_UserNumberFindInterface
    MODULE PROCEDURE Mesh_UserNumberFindRegion
  END INTERFACE Mesh_UserNumberFind

   INTERFACE MESH_USER_NUMBER_FIND
    MODULE PROCEDURE Mesh_UserNumberFindInterface
    MODULE PROCEDURE Mesh_UserNumberFindRegion
  END INTERFACE MESH_USER_NUMBER_FIND

  PUBLIC Mesh_DecompositionGet

  PUBLIC Mesh_MeshElementsGet

  PUBLIC Mesh_MeshNodesGet

  PUBLIC Mesh_RegionGet

  PUBLIC Mesh_UserNumberFind

  PUBLIC Mesh_UserNumberFindGeneric

  PUBLIC MESH_USER_NUMBER_FIND

  PUBLIC Mesh_UserNumberGet

  PUBLIC MESH_TOPOLOGY_ELEMENTS_GET

  PUBLIC MeshTopology_NodesGet


CONTAINS

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the decomposition for a given user number in a mesh. \see OPENCMISS::Iron::cmfe_Mesh_DecompositionGet
  SUBROUTINE Mesh_DecompositionGet(mesh,userNumber,decomposition,err,error,*)

    !Argument variables
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh to get the decomposition for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the decomposition to get.
    TYPE(DecompositionType), POINTER :: decomposition !<On exit, a pointer to the decomposition for the mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Mesh_DecompositionGet",err,error,*998)

    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Mesh is not associated.",err,error,*998)
    IF(.NOT.mesh%meshFinished) CALL FlagError("Mesh has not been finished.",err,error,*998)
    IF(ASSOCIATED(decomposition)) CALL FlagError("Decomposition is already associated.",err,error,*998)
    
    NULLIFY(decomposition)
    CALL Decomposition_UserNumberFind(userNumber,mesh,decomposition,err,error,*999)
    IF(.NOT.ASSOCIATED(decomposition)) THEN
      localError="A decomposition with a user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " do not exist on mesh number "//TRIM(NumberToVString(mesh%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Mesh_DecompositionGet")
    RETURN
999 NULLIFY(decomposition)
998 ERRORSEXITS("Mesh_DecompositionGet",err,error)
    RETURN 1
    
  END SUBROUTINE Mesh_DecompositionGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the mesh elements for a given mesh component.
  SUBROUTINE Mesh_MeshElementsGet(mesh,meshComponentNumber,meshElements,err,error,*)

    !Argument variables
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh to get the elements for
    INTEGER(INTG), INTENT(IN) :: meshComponentNumber !<The mesh component number to get the elements for
    TYPE(MeshElementsType), POINTER :: meshElements !<On return, a pointer to the mesh elements
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Mesh_MeshElementsGet",err,error,*998)
    
    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Mesh is not associated",err,error,*998)
    IF(ASSOCIATED(meshElements)) CALL FlagError("Elements is already associated.",err,error,*998)
    IF(meshComponentNumber<=0.OR.meshComponentNumber>mesh%NUMBER_OF_COMPONENTS) THEN
      localError="The specified mesh component number of "//TRIM(NumberToVString(meshComponentNumber,"*",err,error))// &
        & " is invalid for mesh number "//TRIM(NumberToVString(mesh%userNumber,"*",err,error))// &
        & ". The component number must be between 1 and "//TRIM(NumberToVString(mesh%NUMBER_OF_COMPONENTS,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ASSOCIATED(mesh%topology(meshComponentNumber)%ptr)) THEN
      localError="The mesh topology is not associated for mesh component number "// &
        & TRIM(NumberToVString(meshComponentNumber,"*",err,error))//" of mesh number "// &
        & TRIM(NumberToVString(mesh%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    meshElements=>mesh%topology(meshComponentNumber)%ptr%elements
    IF(.NOT.ASSOCIATED(meshElements)) THEN
      localError="The mesh topology elements is not associated for mesh component number "// &
        & TRIM(NumberToVString(meshComponentNumber,"*",err,error))//" of mesh number "// &
        & TRIM(NumberToVString(mesh%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    EXITS("Mesh_MeshElementsGet")
    RETURN
999 NULLIFY(meshElements)
998 ERRORSEXITS("Mesh_MeshElementsGet",err,error)
    RETURN 1
    
  END SUBROUTINE Mesh_MeshElementsGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the mesh nodes for a given mesh component.
  SUBROUTINE Mesh_MeshNodesGet(mesh,meshComponentNumber,meshNodes,err,error,*)

    !Argument variables
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh to get the nodes for
    INTEGER(INTG), INTENT(IN) :: meshComponentNumber !<The mesh component number to get the nodes for
    TYPE(MeshNodesType), POINTER :: meshNodes !<On return, a pointer to the mesh nodes. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Mesh_MeshNodesGet",err,error,*998)
    
    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Mesh is not associated",err,error,*998)
    IF(ASSOCIATED(meshNodes)) CALL FlagError("Nodes is already associated.",err,error,*998)    
    IF(meshComponentNumber<=0.OR.meshComponentNumber>mesh%NUMBER_OF_COMPONENTS) THEN
      localError="The specified mesh component number of "//TRIM(NumberToVString(meshComponentNumber,"*",err,error))// &
        & " is invalid for mesh number "//TRIM(NumberToVString(mesh%userNumber,"*",err,error))// &
        & ". The component number must be between 1 and "//TRIM(NumberToVString(mesh%NUMBER_OF_COMPONENTS,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ASSOCIATED(mesh%topology(meshComponentNumber)%ptr)) THEN
      localError="The mesh topology is not associated for mesh component number "// &
        & TRIM(NumberToVString(meshComponentNumber,"*",err,error))//" of mesh number "// &
        & TRIM(NumberToVString(mesh%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    meshNodes=>mesh%topology(meshComponentNumber)%ptr%nodes
    IF(.NOT.ASSOCIATED(meshNodes)) THEN
      localError="The mesh topology nodes is not associated for mesh component number "// &
        & TRIM(NumberToVString(meshComponentNumber,"*",err,error))//" of mesh number "// &
        & TRIM(NumberToVString(mesh%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    EXITS("Mesh_MeshNodesGet")
    RETURN
999 NULLIFY(meshNodes)
998 ERRORSEXITS("Mesh_MeshNodesGet",err,error)
    RETURN 1
    
  END SUBROUTINE Mesh_MeshNodesGet

  !
  !================================================================================================================================
  !

  !>Returns the region for a mesh accounting for regions and interfaces
  SUBROUTINE Mesh_RegionGet(mesh,region,err,error,*)

    !Argument variables
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh to get the region for
    TYPE(RegionType), POINTER :: region !<On return, the meshes region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(InterfaceType), POINTER :: interface
    TYPE(RegionType), POINTER :: parentRegion
    TYPE(VARYING_STRING) :: localError

    ENTERS("Mesh_RegionGet",err,error,*999)

    !Check input arguments
    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Mesh is not associated.",err,error,*999)
    IF(ASSOCIATED(region)) CALL FlagError("Region is already associated.",err,error,*999)
    
    NULLIFY(region)
    NULLIFY(interface)
    region=>mesh%region
    IF(.NOT.ASSOCIATED(region)) THEN
      interface=>mesh%interface
      IF(ASSOCIATED(interface)) THEN
        parentRegion=>interface%parentRegion
        IF(ASSOCIATED(parentRegion)) THEN
          region=>parentRegion
        ELSE
          localError="The parent region is not associated for mesh number "// &
            & TRIM(NumberToVString(mesh%userNumber,"*",err,error))//" of interface number "// &
            & TRIM(NumberToVString(interface%userNumber,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        localError="The region or interface is not associated for mesh number "// &
          & TRIM(NumberToVString(mesh%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDIF

    EXITS("Mesh_RegionGet")
    RETURN
999 ERRORSEXITS("Mesh_RegionGet",err,error)
    RETURN 1
    
  END SUBROUTINE Mesh_RegionGet

  !
  !================================================================================================================================
  !

  !>Finds and returns a pointer to the mesh identified by user number in the given list of meshes. If no mesh with that number exits mesh is left nullified.
  SUBROUTINE Mesh_UserNumberFindGeneric(userNumber,meshes,mesh,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the mesh to find
    TYPE(MeshesType), POINTER :: meshes !<The list of meshes containing the mesh.
    TYPE(MeshType), POINTER :: mesh !<On return, a pointer to the mesh of the specified user number. In no mesh with the specified user number exists the pointer is returned NULL. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: meshIdx
    TYPE(VARYING_STRING) :: localError

    ENTERS("Mesh_UserNumberFindGeneric",err,error,*999)

    !Check input arguments
    IF(.NOT.ASSOCIATED(meshes)) CALL FlagError("Meshes is not associated",err,error,*999)
    IF(ASSOCIATED(mesh)) CALL FlagError("Mesh is already associated.",err,error,*999)
    
    !Get the mesh from the user number
    NULLIFY(mesh)
    IF(ASSOCIATED(meshes%meshes)) THEN
      DO meshIdx=1,meshes%numberOfMeshes
        IF(ASSOCIATED(meshes%meshes(meshIdx)%ptr)) THEN
          IF(meshes%meshes(meshIdx)%ptr%userNumber==userNumber) THEN
            mesh=>meshes%meshes(meshIdx)%ptr
            EXIT
          ENDIF
        ELSE
          localError="The mesh pointer in meshes is not associated for mesh index "// &
            & TRIM(NumberToVString(meshIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)          
        ENDIF
      ENDDO !meshIdx      
    ENDIF
    
    EXITS("Mesh_UserNumberFindGeneric")
    RETURN
999 ERRORSEXITS("Mesh_UserNumberFindGeneric",err,error)
    RETURN 1
    
  END SUBROUTINE Mesh_UserNumberFindGeneric

  !
  !================================================================================================================================
  !

  !>Finds and returns a pointer to the mesh identified by user number in the given interface. If no mesh with that number exits mesh is left nullified.
  SUBROUTINE Mesh_UserNumberFindInterface(userNumber,interface,mesh,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the mesh to find
    TYPE(InterfaceType), POINTER :: interface !<A pointer to the interface containing the mesh
    TYPE(MeshType), POINTER :: mesh !<On return, a pointer to the mesh of the specified user number. In no mesh with the specified user number exists the pointer is returned NULL. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
  
    ENTERS("Mesh_UserNumberFindInterface",err,error,*999)

    !Check input arguments
    IF(.NOT.ASSOCIATED(interface)) CALL FlagError("Interface is not associated",err,error,*999)
    
    CALL Mesh_UserNumberFindGeneric(userNumber,interface%meshes,mesh,err,error,*999)
     
    EXITS("Mesh_UserNumberFindInterface")
    RETURN
999 ERRORSEXITS("Mesh_UserNumberFindInterface",err,error)
    RETURN 1
    
  END SUBROUTINE Mesh_UserNumberFindInterface

  !
  !================================================================================================================================
  !

  !>Finds and returns a pointer to the mesh identified by user number in the given. If no mesh with that number exits mesh is left nullified.
  SUBROUTINE Mesh_UserNumberFindRegion(userNumber,region,mesh,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the mesh to find
    TYPE(RegionType), POINTER :: REGION !<The region containing the mesh
    TYPE(MeshType), POINTER :: MESH !<On return, a pointer to the mesh of the specified user number. In no mesh with the specified user number exists the pointer is returned NULL. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("Mesh_UserNumberFindRegion",err,error,*999)

    !Check input arguments
    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated",err,error,*999)
    
    CALL Mesh_UserNumberFindGeneric(userNumber,region%meshes,mesh,err,error,*999)
     
    EXITS("Mesh_UserNumberFindRegion")
    RETURN
999 ERRORSEXITS("Mesh_UserNumberFindRegion",err,error)
    RETURN 1
    
  END SUBROUTINE Mesh_UserNumberFindRegion

  !
  !================================================================================================================================
  !

  !>Returns the user number for a mesh.
  SUBROUTINE Mesh_UserNumberGet(mesh,userNumber,err,error,*)

    !Argument variables
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh to get the user number for
    INTEGER(INTG), INTENT(OUT) :: userNumber !<On exit, the user number of the mesh.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Mesh_UserNumberGet",err,error,*999)

    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Mesh is not associated.",err,error,*999)

    userNumber=mesh%userNumber
  
    EXITS("Mesh_UserNumberGet")
    RETURN
999 ERRORSEXITS("Mesh_UserNumberGet",err,error)
    RETURN 1
    
  END SUBROUTINE Mesh_UserNumberGet

  !
  !================================================================================================================================
  !
  
END MODULE MeshAccessRoutines
