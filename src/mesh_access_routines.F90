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
  USE NodeAccessRoutines
  USE Strings
  USE Trees
  USE Types
  
#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup MESH_ROUTINES_MeshBoundaryTypes MESH_ROUTINES::MeshBoundaryTypes
  !> \brief The types of whether or not a node/element is on a mesh domain boundary.
  !> \see MESH_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: MESH_OFF_DOMAIN_BOUNDARY=0 !<The node/element is not on the mesh domain boundary. \see MESH_ROUTINES_MeshBoundaryTypes,MESH_ROUTINES
  INTEGER(INTG), PARAMETER :: MESH_ON_DOMAIN_BOUNDARY=1 !<The node/element is on the mesh domain boundary. \see MESH_ROUTINES_MeshBoundaryTypes,MESH_ROUTINES
  !>@}
  
  !Module types

  !Module variables

  !Interfaces

  INTERFACE Mesh_UserNumberFind
    MODULE PROCEDURE Mesh_UserNumberFindInterface
    MODULE PROCEDURE Mesh_UserNumberFindRegion
  END INTERFACE Mesh_UserNumberFind

  PUBLIC MESH_ON_DOMAIN_BOUNDARY,MESH_OFF_DOMAIN_BOUNDARY

  PUBLIC Mesh_AssertIsFinished,Mesh_AssertNotFinished

  PUBLIC Mesh_DecompositionGet

  PUBLIC Mesh_DecompositionsGet

  PUBLIC Mesh_GeneratedMeshExists

  PUBLIC Mesh_GeneratedMeshGet

  PUBLIC Mesh_InterfaceGet

  PUBLIC Mesh_IsInterfaceMesh

  PUBLIC Mesh_IsRegionMesh

  PUBLIC Mesh_MeshesGet

  PUBLIC Mesh_MeshElementsGet

  PUBLIC Mesh_MeshNodesGet

  PUBLIC Mesh_NodesGet

  PUBLIC Mesh_NumberOfComponentsGet

  PUBLIC Mesh_NumberOfElementsGet

  PUBLIC Mesh_RegionGet

  PUBLIC Mesh_MeshTopologyGet

  PUBLIC Mesh_UserNumberFind

  PUBLIC Mesh_UserNumberGet

  PUBLIC MeshElements_AssertIsFinished,MeshElements_AssertNotFinished

  PUBLIC MeshElements_BasisGet

  PUBLIC MeshElements_ElementCheckExists

  PUBLIC MeshElements_GlobalElementNumberGet
    
  PUBLIC MeshElements_MeshElementGet

  PUBLIC MeshElements_MeshTopologyGet
  
  PUBLIC MeshElements_UserElementNumberGet

  PUBLIC MeshNodes_GlobalNodeNumberGet
  
  PUBLIC MeshNodes_MeshNodeGet

  PUBLIC MeshNodes_MeshNodeNumberGet

  PUBLIC MeshNodes_MeshTopologyGet
  
  PUBLIC MeshNodes_NodeOnBoundaryGet

  PUBLIC MeshNodes_NodeDerivativesGet

  PUBLIC MeshNodes_NodeNumberOfDerivativesGet

  PUBLIC MeshNodes_NodeNumberOfVersionsGet

  PUBLIC MeshNodes_NodeCheckExists    

  PUBLIC MeshNodes_NumberOfNodesGet

  PUBLIC MeshNodes_UserNodeNumberGet

  PUBLIC MeshTopology_MeshGet

  PUBLIC MeshTopology_MeshDataPointsGet
  
  PUBLIC MeshTopology_MeshDofsGet

  PUBLIC MeshTopology_MeshElementsGet

  PUBLIC MeshTopology_MeshNodesGet

CONTAINS

  !
  !=================================================================================================================================
  !

  !>Assert that a mesh has been finished
  SUBROUTINE Mesh_AssertIsFinished(mesh,err,error,*)

    !Argument Variables
    TYPE(MeshType), POINTER, INTENT(INOUT) :: mesh !<The mesh to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(InterfaceType), POINTER :: interface
    TYPE(RegionType), POINTER :: parentRegion,region
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Mesh_AssertIsFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Mesh is not associated.",err,error,*999)
#endif    

    IF(.NOT.mesh%meshFinished) THEN
      region=>mesh%region
      IF(ASSOCIATED(region)) THEN
        localError="Mesh number "//TRIM(NumberToVString(mesh%userNumber,"*",err,error))// &
          & " on region number "//TRIM(NumberToVString(region%userNumber,"*",err,error))// &
          & " has not been finished."        
      ELSE
        interface=>mesh%interface
        IF(ASSOCIATED(interface)) THEN
          parentRegion=>interface%parentRegion
          IF(ASSOCIATED(parentRegion)) THEN
            localError="Mesh number "//TRIM(NumberToVString(mesh%userNumber,"*",err,error))// &
              & " on interface number "//TRIM(NumberToVString(interface%userNumber,"*",err,error))// &
              & " in parent region number "//TRIM(NumberToVString(parentRegion%userNumber,"*",err,error))// &
              & " has not been finished."        
          ELSE
            localError="Mesh number "//TRIM(NumberToVString(mesh%userNumber,"*",err,error))// &
              & " on interface number "//TRIM(NumberToVString(interface%userNumber,"*",err,error))// &
              & " has not been finished."        
          ENDIF
        ELSE
          localError="Mesh number "//TRIM(NumberToVString(mesh%userNumber,"*",err,error))// &
            & " has not been finished."
        ENDIF
      ENDIF
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Mesh_AssertIsFinished")
    RETURN
999 ERRORSEXITS("Mesh_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE Mesh_AssertIsFinished

  !
  !=================================================================================================================================
  !

  !>Assert that a mesh has not been finished
  SUBROUTINE Mesh_AssertNotFinished(mesh,err,error,*)

    !Argument Variables
    TYPE(MeshType), POINTER, INTENT(INOUT) :: mesh !<The mesh to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(InterfaceType), POINTER :: interface
    TYPE(RegionType), POINTER :: parentRegion,region
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Mesh_AssertNotFinished",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Mesh is not associated.",err,error,*999)
#endif    

    IF(mesh%meshFinished) THEN
      region=>mesh%region
      IF(ASSOCIATED(region)) THEN
        localError="Mesh number "//TRIM(NumberToVString(mesh%userNumber,"*",err,error))// &
          & " on region number "//TRIM(NumberToVString(region%userNumber,"*",err,error))// &
          & " has been finished."        
      ELSE
        interface=>mesh%interface
        IF(ASSOCIATED(interface)) THEN
          parentRegion=>interface%parentRegion
          IF(ASSOCIATED(parentRegion)) THEN
            localError="Mesh number "//TRIM(NumberToVString(mesh%userNumber,"*",err,error))// &
              & " on interface number "//TRIM(NumberToVString(interface%userNumber,"*",err,error))// &
              & " in parent region number "//TRIM(NumberToVString(parentRegion%userNumber,"*",err,error))// &
              & " has been finished."        
          ELSE
            localError="Mesh number "//TRIM(NumberToVString(mesh%userNumber,"*",err,error))// &
              & " on interface number "//TRIM(NumberToVString(interface%userNumber,"*",err,error))// &
              & " has been finished."        
          ENDIF
        ELSE
          localError="Mesh number "//TRIM(NumberToVString(mesh%userNumber,"*",err,error))// &
            & " has been finished."
        ENDIF
      ENDIF
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("Mesh_AssertNotFinished")
    RETURN
999 ERRORSEXITS("Mesh_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE Mesh_AssertNotFinished

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the decomposition for a given user number in a mesh. \see OpenCMISS::Iron::cmfe_Mesh_DecompositionGet
  SUBROUTINE Mesh_DecompositionGet(mesh,userNumber,decomposition,err,error,*)

    !Argument variables
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh to get the decomposition for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the decomposition to get.
    TYPE(DecompositionType), POINTER :: decomposition !<On exit, a pointer to the decomposition for the mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("Mesh_DecompositionGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(decomposition)) CALL FlagError("Decomposition is already associated.",err,error,*998)
    CALL Mesh_AssertIsFinished(mesh,err,error,*999)
#endif    
    
    NULLIFY(decomposition)
    CALL Decomposition_UserNumberFind(userNumber,mesh,decomposition,err,error,*999)

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(decomposition)) THEN
      localError="A decomposition with a user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " do not exist on mesh number "//TRIM(NumberToVString(mesh%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    EXITS("Mesh_DecompositionGet")
    RETURN
999 NULLIFY(decomposition)
998 ERRORSEXITS("Mesh_DecompositionGet",err,error)
    RETURN 1
    
  END SUBROUTINE Mesh_DecompositionGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the decompositions for a given mesh. 
  SUBROUTINE Mesh_DecompositionsGet(mesh,decompositions,err,error,*)

    !Argument variables
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh to get the decompositions for
    TYPE(DecompositionsType), POINTER :: decompositions !<On exit, a pointer to the decompositions for the mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("Mesh_DecompositionsGet",err,error,*998)

#ifdef WITH_PRECHECKS
    IF(ASSOCIATED(decompositions)) CALL FlagError("Decompositions is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Mesh is not associated.",err,error,*999)
#endif    

    decompositions=>mesh%decompositions

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(decompositions)) THEN      
      localError="The decompositions is not associated for mesh number "//TRIM(NumberToVString(mesh%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    EXITS("Mesh_DecompositionsGet")
    RETURN
999 NULLIFY(decompositions)
998 ERRORSEXITS("Mesh_DecompositionsGet",err,error)
    RETURN 1
    
  END SUBROUTINE Mesh_DecompositionsGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the generated mesh for a given mesh to check if the generated mesh eists 
  SUBROUTINE Mesh_GeneratedMeshExists(mesh,generatedMesh,err,error,*)

    !Argument variables
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh to get the generated mesh for
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<On exit, a pointer to the generated mesh for the mesh if it exists. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Mesh_GeneratedMeshExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(generatedMesh)) CALL FlagError("Generated mesh is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Mesh is not associated.",err,error,*999)
#endif    

    generatedMesh=>mesh%generatedMesh
    
    EXITS("Mesh_GeneratedMeshExists")
    RETURN
999 NULLIFY(generatedMesh)
998 ERRORSEXITS("Mesh_GeneratedMeshExists",err,error)
    RETURN 1
    
  END SUBROUTINE Mesh_GeneratedMeshExists
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the generated mesh for a given mesh. 
  SUBROUTINE Mesh_GeneratedMeshGet(mesh,generatedMesh,err,error,*)

    !Argument variables
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh to get the generated mesh for
    TYPE(GeneratedMeshType), POINTER :: generatedMesh !<On exit, a pointer to the generated mesh for the mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("Mesh_GeneratedMeshGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(generatedMesh)) CALL FlagError("Generated mesh is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Mesh is not associated.",err,error,*999)
#endif    

    generatedMesh=>mesh%generatedMesh

#ifdef WITH_POSTCHECKS      
    IF(.NOT.ASSOCIATED(generatedMesh)) THEN      
      localError="The generated mesh is not associated for mesh number "//TRIM(NumberToVString(mesh%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    EXITS("Mesh_GeneratedMeshGet")
    RETURN
999 NULLIFY(generatedMesh)
998 ERRORSEXITS("Mesh_GeneratedMeshGet",err,error)
    RETURN 1
    
  END SUBROUTINE Mesh_GeneratedMeshGet
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the interface for a given mesh. 
  SUBROUTINE Mesh_InterfaceGet(mesh,interface,err,error,*)

    !Argument variables
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh to get the interface for
    TYPE(InterfaceType), POINTER :: interface !<On exit, a pointer to the interface for the mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("Mesh_InterfaceGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(interface)) CALL FlagError("Interface is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Mesh is not associated.",err,error,*999)
#endif    

    INTERFACE=>mesh%INTERFACE

#ifdef WITH_POSTCHECKS      
    IF(.NOT.ASSOCIATED(INTERFACE)) THEN      
      localError="The interface is not associated for mesh number "//TRIM(NumberToVString(mesh%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    EXITS("Mesh_InterfaceGet")
    RETURN
999 NULLIFY(interface)
998 ERRORSEXITS("Mesh_InterfaceGet",err,error)
    RETURN 1
    
  END SUBROUTINE Mesh_InterfaceGet
  
  !
  !================================================================================================================================
  !

  !>Determines if the given mesh is an interface mesh or not. 
  SUBROUTINE Mesh_IsInterfaceMesh(mesh,interfaceMesh,err,error,*)

    !Argument variables
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh to determine if it is an interface mesh or not.
    LOGICAL :: interfaceMesh !<On exit, .TRUE. if the given mesh is in an interface region, .FALSE. if not. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Mesh_IsInterfaceMesh",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Mesh is not associated.",err,error,*999)
#endif    

    interfaceMesh = ASSOCIATED(mesh%interface)
    
    EXITS("Mesh_IsInterfaceMesh")
    RETURN
999 ERRORSEXITS("Mesh_IsInterfaceMesh",err,error)
    RETURN 1
    
  END SUBROUTINE Mesh_IsInterfaceMesh

  !
  !================================================================================================================================
  !

  !>Determines if the given mesh is a region mesh or not. 
  SUBROUTINE Mesh_IsRegionMesh(mesh,regionMesh,err,error,*)

    !Argument variables
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh to determine if it is an region mesh or not.
    LOGICAL :: regionMesh !<On exit, .TRUE. if the given mesh is in a region, .FALSE. if not. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Mesh_IsRegionMesh",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Mesh is not associated.",err,error,*999)
#endif
    
    regionMesh = ASSOCIATED(mesh%region)
    
    EXITS("Mesh_IsRegionMesh")
    RETURN
999 ERRORSEXITS("Mesh_IsRegionMesh",err,error)
    RETURN 1
    
  END SUBROUTINE Mesh_IsRegionMesh

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the meshes for a given mesh. 
  SUBROUTINE Mesh_MeshesGet(mesh,meshes,err,error,*)

    !Argument variables
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh to get the meshes for
    TYPE(MeshesType), POINTER :: meshes !<On exit, a pointer to the meshes for the mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif
    
    ENTERS("Mesh_MeshesGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(meshes)) CALL FlagError("Meshes is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Mesh is not associated.",err,error,*999)
#endif
    
    meshes=>mesh%meshes

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(meshes)) THEN      
      localError="The meshes is not associated for mesh number "//TRIM(NumberToVString(mesh%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    EXITS("Mesh_MeshesGet")
    RETURN
999 NULLIFY(meshes)
998 ERRORSEXITS("Mesh_MeshesGet",err,error)
    RETURN 1
    
  END SUBROUTINE Mesh_MeshesGet

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
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif
    
    ENTERS("Mesh_MeshElementsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Mesh is not associated",err,error,*998)
    IF(ASSOCIATED(meshElements)) CALL FlagError("Elements is already associated.",err,error,*998)
    IF(meshComponentNumber<=0.OR.meshComponentNumber>mesh%numberOfComponents) THEN
      localError="The specified mesh component number of "//TRIM(NumberToVString(meshComponentNumber,"*",err,error))// &
        & " is invalid for mesh number "//TRIM(NumberToVString(mesh%userNumber,"*",err,error))// &
        & ". The component number must be between 1 and "//TRIM(NumberToVString(mesh%numberOfComponents,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ASSOCIATED(mesh%topology(meshComponentNumber)%ptr)) THEN
      localError="The mesh topology is not associated for mesh component number "// &
        & TRIM(NumberToVString(meshComponentNumber,"*",err,error))//" of mesh number "// &
        & TRIM(NumberToVString(mesh%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    meshElements=>mesh%topology(meshComponentNumber)%ptr%elements

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(meshElements)) THEN
      localError="The mesh topology elements is not associated for mesh component number "// &
        & TRIM(NumberToVString(meshComponentNumber,"*",err,error))//" of mesh number "// &
        & TRIM(NumberToVString(mesh%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

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
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("Mesh_MeshNodesGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Mesh is not associated",err,error,*998)
    IF(ASSOCIATED(meshNodes)) CALL FlagError("Nodes is already associated.",err,error,*998)    
    IF(meshComponentNumber<=0.OR.meshComponentNumber>mesh%numberOfComponents) THEN
      localError="The specified mesh component number of "//TRIM(NumberToVString(meshComponentNumber,"*",err,error))// &
        & " is invalid for mesh number "//TRIM(NumberToVString(mesh%userNumber,"*",err,error))// &
        & ". The component number must be between 1 and "//TRIM(NumberToVString(mesh%numberOfComponents,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ASSOCIATED(mesh%topology(meshComponentNumber)%ptr)) THEN
      localError="The mesh topology is not associated for mesh component number "// &
        & TRIM(NumberToVString(meshComponentNumber,"*",err,error))//" of mesh number "// &
        & TRIM(NumberToVString(mesh%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    meshNodes=>mesh%topology(meshComponentNumber)%ptr%nodes

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(meshNodes)) THEN
      localError="The mesh topology nodes is not associated for mesh component number "// &
        & TRIM(NumberToVString(meshComponentNumber,"*",err,error))//" of mesh number "// &
        & TRIM(NumberToVString(mesh%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    EXITS("Mesh_MeshNodesGet")
    RETURN
999 NULLIFY(meshNodes)
998 ERRORSEXITS("Mesh_MeshNodesGet",err,error)
    RETURN 1
    
  END SUBROUTINE Mesh_MeshNodesGet

  !
  !================================================================================================================================
  !

  !>Returns a region nodes pointer corresponding to the mesh global nodes accounting for interfaces.
  SUBROUTINE Mesh_NodesGet(mesh,nodes,err,error,*)

    !Argument variables
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh to get the global nodes for
    TYPE(NodesType), POINTER :: nodes !<On return, the nodes pointer corresponding to the global nodes for the mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(InterfaceType), POINTER :: interface
    TYPE(RegionType), POINTER :: region
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Mesh_NodesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(nodes)) CALL FlagError("Nodes is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Mesh is not associated.",err,error,*999)
#endif
    
    region=>mesh%region
    IF(ASSOCIATED(region)) THEN
      nodes=>region%nodes
    ELSE
      INTERFACE=>mesh%INTERFACE
      IF(ASSOCIATED(interface)) THEN
        nodes=>interface%nodes
      ELSE
        localError="Mesh number "//TRIM(NumberToVString(mesh%userNumber,"*",err,error))// &
          & " does not have an associated region or interface."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDIF

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(nodes)) THEN
      IF(ASSOCIATED(region)) THEN
        localError="Mesh number "//TRIM(NumberToVString(mesh%userNumber,"*",err,error))// &
          & " does not have any nodes associated with the mesh region."
      ELSE
        localError="Mesh number "//TRIM(NumberToVString(mesh%userNumber,"*",err,error))// &
          & " does not have any nodes associated with the mesh interface."
      ENDIF
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
    EXITS("Mesh_NodesGet")
    RETURN
999 ERRORSEXITS("Mesh_NodesGet",err,error)
    RETURN 1
    
  END SUBROUTINE Mesh_NodesGet
  
  !
  !================================================================================================================================
  !
  
  !>Gets the number of mesh components for a mesh identified by a pointer. \see OpenCMISS::Iron::cmfe_Mesh_NumberOfComponentsGet
  SUBROUTINE Mesh_NumberOfComponentsGet(mesh,numberOfComponents,err,error,*)

    !Argument variables
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh to get the number of components for
    INTEGER(INTG), INTENT(OUT) :: numberOfComponents !<On return, the number of components in the specified mesh.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("Mesh_NumberOfComponentsGet",err,error,*999)

    CALL Mesh_AssertIsFinished(mesh,err,error,*999)
    
    numberOfComponents=mesh%numberOfComponents
     
    EXITS("Mesh_NumberOfComponentsGet")
    RETURN
999 ERRORSEXITS("Mesh_NumberOfComponentsGet",err,error)    
    RETURN 1
    
  END SUBROUTINE Mesh_NumberOfComponentsGet

  !
  !================================================================================================================================
  !
  
  !>Gets the number of elements for a mesh identified by a pointer. \see OpenCMISS::Iron::cmfe_Mesh_NumberOfElementsGet
  SUBROUTINE Mesh_NumberOfElementsGet(mesh,numberOfElements,err,error,*)

    !Argument variables
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh to get the number of elements for
    INTEGER(INTG), INTENT(OUT) :: numberOfElements !<On return, the number of elements in the specified mesh
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Mesh_NumberOfElementsGet",err,error,*999)

    CALL Mesh_AssertIsFinished(mesh,err,error,*999)

    numberOfElements=mesh%numberOfElements
    
    EXITS("Mesh_NumberOfElementsGet")
    RETURN
999 ERRORSEXITS("Mesh_NumberOfElementsGet",err,error)    
    RETURN 1
    
  END SUBROUTINE Mesh_NumberOfElementsGet

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

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Mesh is not associated.",err,error,*999)
    IF(ASSOCIATED(region)) CALL FlagError("Region is already associated.",err,error,*999)
#endif    
    
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

  !>Returns the mesh topology for a mesh and mesh component number.
  SUBROUTINE Mesh_MeshTopologyGet(mesh,meshComponentNumber,meshTopology,err,error,*)

    !Argument variables
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh to get the mesh topology for
    INTEGER(INTG), INTENT(IN) :: meshComponentNumber !<The mesh component number of the topology to get
    TYPE(MeshTopologyType), POINTER :: meshTopology !<On return, the mesh topology for the mesh and mesh component number. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("Mesh_MeshTopologyGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(ASSOCIATED(meshTopology)) CALL FlagError("Mesh topology is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Mesh is not associated.",err,error,*998)
    IF(meshComponentNumber<0.OR.meshComponentNumber>mesh%numberOfComponents) THEN
      localError="The specified mesh component number of "//TRIM(NumberToVString(meshComponentNumber,"*",err,error))// &
        & " is invalid. The mesh component number must be >= 1 and <= "// &
        & TRIM(NumberToVString(mesh%numberOfComponents,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(mesh%topology)) THEN
      localError="Topology is not allocated for mesh number "//TRIM(NumberToVString(mesh%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    meshTopology=>mesh%topology(meshComponentNumber)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(meshTopology)) THEN    
      localError="The mesh topology is not associated for mesh component number "// &
        & TRIM(NumberToVString(meshComponentNumber,"*",err,error))//" of mesh number "// &
        & TRIM(NumberToVString(mesh%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    EXITS("Mesh_MeshTopologyGet")
    RETURN
999 NULLIFY(meshTopology)
998 ERRORSEXITS("Mesh_MeshTopologyGet",err,error)
    RETURN 1
    
  END SUBROUTINE Mesh_MeshTopologyGet

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

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(.NOT.ASSOCIATED(meshes)) CALL FlagError("Meshes is not associated.",err,error,*999)
    IF(ASSOCIATED(mesh)) CALL FlagError("Mesh is already associated.",err,error,*999)
#endif    
    
    !Get the mesh from the user number
    NULLIFY(mesh)
    IF(ALLOCATED(meshes%meshes)) THEN
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

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(.NOT.ASSOCIATED(INTERFACE)) CALL FlagError("Interface is not associated",err,error,*999)
#endif    
    
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

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated",err,error,*999)
#endif
    
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

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Mesh is not associated.",err,error,*999)
#endif
    
    userNumber=mesh%userNumber
  
    EXITS("Mesh_UserNumberGet")
    RETURN
999 ERRORSEXITS("Mesh_UserNumberGet",err,error)
    RETURN 1
    
  END SUBROUTINE Mesh_UserNumberGet

  !
  !=================================================================================================================================
  !

  !>Assert that a mesh elements has been finished
  SUBROUTINE MeshElements_AssertIsFinished(meshElements,err,error,*)

    !Argument Variables
    TYPE(MeshElementsType), POINTER, INTENT(INOUT) :: meshElements !<The mesh elements to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: meshComponentNumber
    TYPE(MeshType), POINTER :: mesh
    TYPE(MeshTopologyType), POINTER :: meshTopology
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Mesh_AssertIsFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(meshElements)) CALL FlagError("Mesh elements is not associated.",err,error,*999)
#endif    

    IF(.NOT.meshElements%elementsFinished) THEN
      meshTopology=>meshElements%meshTopology
      IF(ASSOCIATED(meshTopology)) THEN
        meshComponentNumber=meshTopology%meshComponentNumber
        mesh=>meshTopology%mesh
        IF(ASSOCIATED(mesh)) THEN
          localError="Mesh elements of mesh component number "//TRIM(NumberToVString(meshComponentNumber,"*",err,error))// &
            & " of mesh number "//TRIM(NumberToVString(mesh%userNumber,"*",err,error))// &
            & " have not been finished."        
        ELSE
          localError="Mesh elements of mesh component number "//TRIM(NumberToVString(meshComponentNumber,"*",err,error))// &
            & " have not been finished."        
        ENDIF
      ELSE
        localError="Mesh elements have not been finished."        
      ENDIF
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("MeshElements_AssertIsFinished")
    RETURN
999 ERRORSEXITS("MeshElements_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE MeshElements_AssertIsFinished

  !
  !=================================================================================================================================
  !

  !>Assert that a mesh elements has not been finished
  SUBROUTINE MeshElements_AssertNotFinished(meshElements,err,error,*)

    !Argument Variables
    TYPE(MeshElementsType), POINTER, INTENT(INOUT) :: meshElements !<The mesh elements to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: meshComponentNumber
    TYPE(MeshType), POINTER :: mesh
    TYPE(MeshTopologyType), POINTER :: meshTopology
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("MeshElements_AssertNotFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(meshElements)) CALL FlagError("Mesh elements is not associated.",err,error,*999)
#endif
    
    IF(meshElements%elementsFinished) THEN
      meshTopology=>meshElements%meshTopology
      IF(ASSOCIATED(meshTopology)) THEN
        meshComponentNumber=meshTopology%meshComponentNumber
        mesh=>meshTopology%mesh
        IF(ASSOCIATED(mesh)) THEN
          localError="Mesh elements of mesh component number "//TRIM(NumberToVString(meshComponentNumber,"*",err,error))// &
            & " of mesh number "//TRIM(NumberToVString(mesh%userNumber,"*",err,error))// &
            & " have been finished."        
        ELSE
          localError="Mesh elements of mesh component number "//TRIM(NumberToVString(meshComponentNumber,"*",err,error))// &
            & " have been finished."        
        ENDIF
      ELSE
        localError="Mesh elements have been finished."        
      ENDIF
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("MeshElements_AssertNotFinished")
    RETURN
999 ERRORSEXITS("MeshElements_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE MeshElements_AssertNotFinished

  !  
  !================================================================================================================================
  !

  !>Get the basis for an element in the mesh elements identified by its global number
  SUBROUTINE MeshElements_BasisGet(meshElements,globalElementNumber,basis,err,error,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: meshElements !<A pointer to the mesh elements to get the element basis for
    INTEGER(INTG), INTENT(IN) :: globalElementNumber !<The element global number to get the basis for
    TYPE(BasisType), POINTER, INTENT(OUT) :: basis !<On return, a pointer to the basis for the element. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif
    
    ENTERS("MeshElements_BasisGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(basis)) CALL FlagError("Basis is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(meshElements)) CALL FlagError("Mesh elements is not associated.",err,error,*999)
    IF(globalElementNumber<=0.OR.globalElementNumber>meshElements%numberOfElements) THEN
      localError="The global element number of "//TRIM(NumberToVString(globalElementNumber,"*",err,error))// &
        & " is invalid. The global number must be >= 1 and <= "// &
        & TRIM(NumberToVString(meshElements%numberOfElements,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(meshElements%elements)) CALL FlagError("Mesh elements elements has not been allocated.",err,error,*999)
#endif
    
    basis=>meshElements%elements(globalElementNumber)%basis

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(basis)) THEN
      localError="The basis for global element number "//TRIM(NumberToVString(globalElementNumber,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    EXITS("MeshElements_BasisGet")
    RETURN
999 NULLIFY(basis)
998 ERRORSEXITS("MeshElements_BasisGet",err,error)
    RETURN 1
    
  END SUBROUTINE MeshElements_BasisGet

  !
  !================================================================================================================================
  !

  !>Checks that a user element number exists in a mesh elements. 
  SUBROUTINE MeshElements_ElementCheckExists(meshElements,userElementNumber,elementExists,globalElementNumber,err,error,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: meshElements !<A pointer to the mesh elements to check the element exists on
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number to check if it exists
    LOGICAL, INTENT(OUT) :: elementExists !<On exit, is .TRUE. if the element user number exists in the mesh elements, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: globalElementNumber !<On exit, if the element exists the global number corresponding to the user element number. If the element does not exist then global number will be 0.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(TreeNodeType), POINTER :: treeNode
    
    ENTERS("MeshElements_ElementCheckExists",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(meshElements)) CALL FlagError("Mesh elements is not associated.",err,error,*999)
#endif
    
    elementExists=.FALSE.
    globalElementNumber=0
    NULLIFY(treeNode)
    CALL Tree_Search(meshElements%elementsTree,userElementNumber,treeNode,err,error,*999)
    IF(ASSOCIATED(treeNode)) THEN
      CALL Tree_NodeValueGet(meshElements%elementsTree,treeNode,globalElementNumber,err,error,*999)
      elementExists=.TRUE.
    ENDIF
    
    EXITS("MeshElements_ElementCheckExists")
    RETURN
999 ERRORSEXITS("MeshElements_ElementCheckExists",err,error)
    RETURN 1
    
  END SUBROUTINE MeshElements_ElementCheckExists
  
  !  
  !================================================================================================================================
  !

  !>Get the global element number in the mesh elements from a user element number
  SUBROUTINE MeshElements_GlobalElementNumberGet(meshElements,userElementNumber,globalElementNumber,err,error,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: meshElements !<A pointer to the mesh elements to get the element for
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number to get the global element number for
    INTEGER(INTG), INTENT(OUT) :: globalElementNumber !<On return, the global element number corresponding to the user element number
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    LOGICAL :: elementExists
    TYPE(MeshType), POINTER :: mesh
    TYPE(MeshTopologyType), POINTER :: meshTopology
    TYPE(VARYING_STRING) :: localError

    ENTERS("MeshElements_GlobalElementNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(meshElements)) CALL FlagError("Mesh elements is not associated.",err,error,*999)
#endif    
    
    CALL MeshElements_ElementCheckExists(meshElements,userElementNumber,elementExists,globalElementNumber,err,error,*999)
    IF(.NOT.elementExists) THEN
      meshTopology=>meshElements%meshTopology
      IF(ASSOCIATED(meshTopology)) THEN
        mesh=>meshTopology%mesh
        IF(ASSOCIATED(mesh)) THEN
          localError="The user element number "//TRIM(NumberToVString(userElementNumber,"*",err,error))// &
            & " does not exist in mesh number "//TRIM(NumberToVString(mesh%userNumber,"*",err,error))//"."
        ELSE
          localError="The user element number "//TRIM(NumberToVString(userElementNumber,"*",err,error))//" does not exist."
        ENDIF
      ELSE
        localError="The user element number "//TRIM(NumberToVString(userElementNumber,"*",err,error))//" does not exist."
      ENDIF
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("MeshElements_GlobalElementNumberGet")
    RETURN
999 ERRORSEXITS("MeshElements_GlobalElementNumberGet",err,error)
    RETURN 1
    
  END SUBROUTINE MeshElements_GlobalElementNumberGet

  !  
  !================================================================================================================================
  !

  !>Get the mesh element in the mesh elements identified by its global number
  SUBROUTINE MeshElements_MeshElementGet(meshElements,globalElementNumber,meshElement,err,error,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: meshElements !<A pointer to the mesh elements to get the element for
    INTEGER(INTG), INTENT(IN) :: globalElementNumber !<The element global number to get the element for
    TYPE(MeshElementType), POINTER :: meshElement !<On return, a pointer to the mesh element. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif
    
    ENTERS("MeshElements_MeshElementGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(meshElement)) CALL FlagError("Mesh element is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(meshElements)) CALL FlagError("Mesh elements is not associated.",err,error,*999)
    IF(globalElementNumber<=0.OR.globalElementNumber>meshElements%numberOfElements) THEN
      localError="The global element number of "//TRIM(NumberToVString(globalElementNumber,"*",err,error))// &
        & " is invalid. The global number must be >= 1 and <= "// &
        & TRIM(NumberToVString(meshElements%numberOfElements,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(meshElements%elements)) CALL FlagError("Mesh elements elements has not been allocated.",err,error,*999)
#endif
    
    meshElement=>meshElements%elements(globalElementNumber)

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(meshElement)) THEN
      localError="The mesh element for global element number "//TRIM(NumberToVString(globalElementNumber,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    EXITS("MeshElements_MeshElementGet")
    RETURN
999 NULLIFY(meshElement)
998 ERRORSEXITS("MeshElements_MeshElementGet",err,error)
    RETURN 1
    
  END SUBROUTINE MeshElements_MeshElementGet

  !  
  !================================================================================================================================
  !

  !>Get the mesh topology for mesh elements 
  SUBROUTINE MeshElements_MeshTopologyGet(meshElements,meshTopology,err,error,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER, INTENT(IN) :: meshElements !<A pointer to the mesh elements to get the mesh topology for
    TYPE(MeshTopologyType), POINTER, INTENT(INOUT) :: meshTopology !<On return, a pointer to the mesh topology for the mesh elements. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshElements_MeshTopologyGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(meshTopology)) CALL FlagError("Mesh topology is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(meshElements)) CALL FlagError("Mesh elements is not associated.",err,error,*999)
#endif    
    
    meshTopology=>meshElements%meshTopology

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(meshTopology)) CALL FlagError("Mesh topology is not associated for the mesh elements.",err,error,*999)
#endif
    
    EXITS("MeshElements_MeshTopologyGet")
    RETURN
999 NULLIFY(meshTopology)
998 ERRORSEXITS("MeshElements_MeshTopologyGet",err,error)
    RETURN 1
    
  END SUBROUTINE MeshElements_MeshTopologyGet

  !
  !================================================================================================================================
  !
  
  !>Returns the user element number for an element in mesh elements
  SUBROUTINE MeshElements_UserElementNumberGet(meshElements,elementNumber,userElementNumber,err,error,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: meshElements !<A pointer to the mesh elements containing the elemnet to get the user number for
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number of the element to get the user number for
    INTEGER(INTG), INTENT(OUT) :: userElementNumber !<On return, the user element number of the specified element
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("MeshElements_UserElementNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(meshElements)) CALL FlagError("Mesh elements is not associated.",err,error,*999)
    IF(elementNumber<1.OR.elementNumber>meshElements%numberOfElements) THEN
      localError="The specified mesh element number of "//TRIM(NumberToVString(elementNumber,"*",err,error))// &
        & " is invalid. The element number must be >= 1 and <= "// &
        & TRIM(NumberToVString(meshElements%numberOfElements,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(meshElements%elements)) CALL FlagError("Mesh elements elements is not allocated.",err,error,*999)
#endif
    
    userElementNumber=meshElements%elements(elementNumber)%userNumber
        
    EXITS("MeshElemnets_UserElementNumberGet")
    RETURN
999 ERRORSEXITS("MeshElements_UserElementNumberGet",err,error)    
    RETURN 1
   
  END SUBROUTINE MeshElements_UserElementNumberGet
  
  !  
  !================================================================================================================================
  !

  !>Get the mesh node in the mesh nodes identified by its mesh number
  SUBROUTINE MeshNodes_MeshNodeGet(meshNodes,meshNodeNumber,meshNode,err,error,*)

    !Argument variables
    TYPE(MeshNodesType), POINTER :: meshNodes !<A pointer to the mesh nodes to get the node for
    INTEGER(INTG), INTENT(IN) :: meshNodeNumber !<The mesh node number to get the node for
    TYPE(MeshNodeType), POINTER :: meshNode !<On return, a pointer to the mesh node. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("MeshNodes_MeshNodeGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(meshNode)) CALL FlagError("Mesh node is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(meshNodes)) CALL FlagError("Mesh nodes is not associated.",err,error,*999)
    IF(meshNodeNumber<=0.OR.meshNodeNumber>meshNodes%numberOfNodes) THEN
      localError="The mesh node number of "//TRIM(NumberToVString(meshNodeNumber,"*",err,error))// &
        & " is invalid. The mesh number must be >= 1 and <= "// &
        & TRIM(NumberToVString(meshNodes%numberOfNodes,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(meshNodes%nodes)) CALL FlagError("Mesh nodes nodes has not been allocated.",err,error,*999)
#endif
    
    meshNode=>meshNodes%nodes(meshNodeNumber)

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(meshNode)) THEN
      localError="The mesh node for mesh node number "//TRIM(NumberToVString(meshNodeNumber,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    EXITS("MeshNodes_MeshNodeGet")
    RETURN
999 NULLIFY(meshNode)
998 ERRORSEXITS("MeshNodes_MeshNodeGet",err,error)
    RETURN 1
    
  END SUBROUTINE MeshNodes_MeshNodeGet

  !  
  !================================================================================================================================
  !

  !>Get the global node number in the mesh nodes from a user node number
  SUBROUTINE MeshNodes_GlobalNodeNumberGet(meshNodes,userNodeNumber,globalNodeNumber,err,error,*)

    !Argument variables
    TYPE(MeshNodesType), POINTER :: meshNodes !<A pointer to the mesh nodes to get the node for
    INTEGER(INTG), INTENT(IN) :: userNodeNumber !<The user node number to get the global node number for
    INTEGER(INTG), INTENT(OUT) :: globalNodeNumber !<On return, the global node number corresponding to the user node number
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: meshComponentNumber,meshNodeNumber
    LOGICAL :: nodeExists
    TYPE(MeshType), POINTER :: mesh
    TYPE(MeshTopologyType), POINTER :: meshTopology
    TYPE(VARYING_STRING) :: localError

    ENTERS("MeshNodes_GlobalNodeNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(meshNodes)) CALL FlagError("Mesh nodes is not associated.",err,error,*999)
#endif
    
    CALL MeshNodes_NodeCheckExists(meshNodes,userNodeNumber,nodeExists,globalNodeNumber,meshNodeNumber,err,error,*999)
    IF(.NOT.nodeExists) THEN
      meshTopology=>meshNodes%meshTopology
      IF(ASSOCIATED(meshTopology)) THEN
        meshComponentNumber=meshTopology%meshComponentNumber
        mesh=>meshTopology%mesh
        IF(ASSOCIATED(mesh)) THEN
          localError="The user node number "//TRIM(NumberToVString(userNodeNumber,"*",err,error))// &
            & " does not exist in mesh component number "//TRIM(NumberToVString(meshComponentNumber,"*",err,error))// &
            & " of mesh number "//TRIM(NumberToVString(mesh%userNumber,"*",err,error))//"."
        ELSE
          localError="The user node number "//TRIM(NumberToVString(userNodeNumber,"*",err,error))// &
            & " does not exist in mesh component number "//TRIM(NumberToVString(meshComponentNumber,"*",err,error))//"."
        ENDIF
      ELSE
        localError="The user node number "//TRIM(NumberToVString(userNodeNumber,"*",err,error))//" does not exist."
      ENDIF
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("MeshNodes_GlobalNodeNumberGet")
    RETURN
999 ERRORSEXITS("MeshNodes_GlobalNodeNumberGet",err,error)
    RETURN 1
    
  END SUBROUTINE MeshNodes_GlobalNodeNumberGet

  !  
  !================================================================================================================================
  !

  !>Get the mesh node number in the mesh nodes from a user node number
  SUBROUTINE MeshNodes_MeshNodeNumberGet(meshNodes,userNodeNumber,meshNodeNumber,err,error,*)

    !Argument variables
    TYPE(MeshNodesType), POINTER :: meshNodes !<A pointer to the mesh nodes to get the node for
    INTEGER(INTG), INTENT(IN) :: userNodeNumber !<The user node number to get the mesh node number for
    INTEGER(INTG), INTENT(OUT) :: meshNodeNumber !<On return, the mesh node number corresponding to the user node number
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: globalNodeNumber,meshComponentNumber
    LOGICAL :: nodeExists
    TYPE(MeshType), POINTER :: mesh
    TYPE(MeshTopologyType), POINTER :: meshTopology
    TYPE(VARYING_STRING) :: localError

    ENTERS("MeshNodes_MeshNodeNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(meshNodes)) CALL FlagError("Mesh nodes is not associated.",err,error,*999)
#endif
    
    CALL MeshNodes_NodeCheckExists(meshNodes,userNodeNumber,nodeExists,globalNodeNumber,meshNodeNumber,err,error,*999)
    IF(.NOT.nodeExists) THEN
      meshTopology=>meshNodes%meshTopology
      IF(ASSOCIATED(meshTopology)) THEN
        meshComponentNumber=meshTopology%meshComponentNumber
        mesh=>meshTopology%mesh
        IF(ASSOCIATED(mesh)) THEN
          localError="The user node number "//TRIM(NumberToVString(userNodeNumber,"*",err,error))// &
            & " does not exist in mesh component number "//TRIM(NumberToVString(meshComponentNumber,"*",err,error))// &
            & " of mesh number "//TRIM(NumberToVString(mesh%userNumber,"*",err,error))//"."
        ELSE
          localError="The user node number "//TRIM(NumberToVString(userNodeNumber,"*",err,error))// &
            & " does not exist in mesh component number "//TRIM(NumberToVString(meshComponentNumber,"*",err,error))//"."
        ENDIF
      ELSE
        localError="The user node number "//TRIM(NumberToVString(userNodeNumber,"*",err,error))//" does not exist."
      ENDIF
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("MeshNodes_MeshNodeNumberGet")
    RETURN
999 ERRORSEXITS("MeshNodes_MeshNodeNumberGet",err,error)
    RETURN 1
    
  END SUBROUTINE MeshNodes_MeshNodeNumberGet

  !  
  !================================================================================================================================
  !

  !>Get the mesh topology for mesh nodes 
  SUBROUTINE MeshNodes_MeshTopologyGet(meshNodes,meshTopology,err,error,*)

    !Argument variables
    TYPE(MeshNodesType), POINTER, INTENT(IN) :: meshNodes !<A pointer to the mesh nodes to get the mesh topology for
    TYPE(MeshTopologyType), POINTER, INTENT(INOUT) :: meshTopology !<On return, a pointer to the mesh topology for the mesh nodes. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshNodes_MeshTopologyGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(meshTopology)) CALL FlagError("Mesh topology is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(meshNodes)) CALL FlagError("Mesh nodes is not associated.",err,error,*999)
#endif
    
    meshTopology=>meshNodes%meshTopology

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(meshTopology)) CALL FlagError("Mesh topology is not associated for the mesh nodes.",err,error,*999)
#endif
    
    EXITS("MeshNodes_MeshTopologyGet")
    RETURN
999 NULLIFY(meshTopology)
998 ERRORSEXITS("MeshNodes_MeshTopologyGet",err,error)
    RETURN 1
    
  END SUBROUTINE MeshNodes_MeshTopologyGet

  !
  !================================================================================================================================
  !
  
  !>Returns if the node in a mesh is on the boundary or not
  SUBROUTINE MeshNodes_NodeOnBoundaryGet(meshNodes,userNodeNumber,onBoundary,err,error,*)

    !Argument variables
    TYPE(MeshNodesType), POINTER :: meshNodes !<A pointer to the mesh nodes containing the nodes to get the boundary type for
    INTEGER(INTG), INTENT(IN) :: userNodeNumber !<The user node number to get the boundary type for
    INTEGER(INTG), INTENT(OUT) :: onBoundary !<On return, the boundary type of the specified user node number.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: meshNodeNumber
 
    ENTERS("MeshNodes_NodeOnBoundaryGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(meshNodes)) CALL FlagError("Mesh nodes is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(meshNodes%nodes)) CALL FlagError("Mesh nodes nodes is not allocated.",err,error,*999)
#endif
    
    CALL MeshNodes_MeshNodeNumberGet(meshNodes,userNodeNumber,meshNodeNumber,err,error,*999)
    IF(meshNodes%nodes(meshNodeNumber)%boundaryNode) THEN
      onBoundary=MESH_ON_DOMAIN_BOUNDARY
    ELSE
      onBoundary=MESH_OFF_DOMAIN_BOUNDARY
    ENDIF
    
    EXITS("MeshNodes_NodeOnBoundaryGet")
    RETURN
999 ERRORSEXITS("MeshNodes_NodeOnBoundaryGet",err,error)    
    RETURN 1
   
  END SUBROUTINE MeshNodes_NodeOnBoundaryGet
  
  !
  !================================================================================================================================
  !

  !>Checks that a user node number exists in a mesh component. 
  SUBROUTINE MeshNodes_NodeCheckExists(meshNodes,userNodeNumber,nodeExists,globalNodeNumber,meshNodeNumber,err,error,*)

    !Argument variables
    TYPE(MeshNodesType), POINTER :: meshNodes !<A pointer to the mesh nodes to check the node exists on
    INTEGER(INTG), INTENT(IN) :: userNodeNumber !<The user node number to check if it exists
    LOGICAL, INTENT(OUT) :: nodeExists !<On exit, is .TRUE. if the node user number exists in the mesh component, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: globalNodeNumber !<On exit, if the node exists the global number corresponding to the user node number. If the node does not exist then global number will be 0.
    INTEGER(INTG), INTENT(OUT) :: meshNodeNumber !<On exit, if the node exists the mesh number corresponding to the user node number. If the node does not exist then mesh number will be 0.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(MeshType), POINTER :: mesh
    TYPE(MeshTopologyType), POINTER :: meshTopology
    TYPE(NodesType), POINTER :: nodes
    TYPE(TreeNodeType), POINTER :: treeNode
    
    ENTERS("MeshNodes_NodeCheckExists",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(meshNodes)) CALL FlagError("Mesh nodes is not associated.",err,error,*999)
#endif
    
    nodeExists=.FALSE.
    meshNodeNumber=0

    NULLIFY(meshTopology)
    CALL MeshNodes_MeshTopologyGet(meshNodes,meshTopology,err,error,*999)
    NULLIFY(mesh)
    CALL MeshTopology_MeshGet(meshTopology,mesh,err,error,*999)
    NULLIFY(nodes)
    CALL Mesh_NodesGet(mesh,nodes,err,error,*999)
    CALL Nodes_NodeCheckExists(nodes,userNodeNumber,nodeExists,globalNodeNumber,err,error,*999)
    IF(nodeExists) THEN
      NULLIFY(treeNode)
      CALL Tree_Search(meshNodes%nodesTree,globalNodeNumber,treeNode,err,error,*999)
      IF(ASSOCIATED(treeNode)) THEN
        CALL Tree_NodeValueGet(meshNodes%nodesTree,treeNode,meshNodeNumber,err,error,*999)
        nodeExists=.TRUE.
      ELSE
        nodeExists=.FALSE.
      ENDIF
    ENDIF
    
    EXITS("MeshNodes_NodeCheckExists")
    RETURN
999 ERRORSEXITS("MeshNodes_NodeCheckExists",err,error)
    RETURN 1
    
  END SUBROUTINE MeshNodes_NodeCheckExists
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of derivatives for a node in a mesh
  SUBROUTINE MeshNodes_NodeNumberOfDerivativesGet(meshNodes,userNodeNumber,numberOfDerivatives,err,error,*)

    !Argument variables
    TYPE(MeshNodesType), POINTER :: meshNodes !<A pointer to the mesh nodes containing the nodes to get the number of derivatives for
    INTEGER(INTG), INTENT(IN) :: userNodeNumber !<The user node number to get the number of derivatives for
    INTEGER(INTG), INTENT(OUT) :: numberOfDerivatives !<On return, the number of global derivatives at the specified user node number.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: meshNodeNumber

    ENTERS("MeshNodes_NodeNumberOfDerivativesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(meshNodes)) CALL FlagError("Mesh nodes is not associated.",err,error,*999)
#endif
    
    CALL MeshNodes_MeshNodeNumberGet(meshNodes,userNodeNumber,meshNodeNumber,err,error,*999)
    
    numberOfDerivatives=meshNodes%nodes(meshNodeNumber)%numberOfDerivatives
    
    EXITS("MeshNodes_NodeNumberOfDerivativesGet")
    RETURN
999 ERRORSEXITS("MeshNodes_NodeNumberOfDerivativesGet",err,error)    
    RETURN 1
   
  END SUBROUTINE MeshNodes_NodeNumberOfDerivativesGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the global derivative numbers for a node in mesh nodes
  SUBROUTINE MeshNodes_NodeDerivativesGet(meshNodes,userNodeNumber,derivatives,err,error,*)

    !Argument variables
    TYPE(MeshNodesType), POINTER :: meshNodes !<A pointer to the mesh nodes containing the node to get the derivatives for
    INTEGER(INTG), INTENT(IN) :: userNodeNumber !<The user number of the node to get the derivatives for
    INTEGER(INTG), INTENT(OUT) :: derivatives(:) !<On return, the global derivatives at the specified node
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: derivativeIdx,meshNodeNumber,numberOfDerivatives
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("MeshNodes_NodeDerivativesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(meshNodes)) CALL FlagError("Mesh nodes is not associated.",err,error,*999)
#endif
    
    CALL MeshNodes_MeshNodeNumberGet(meshNodes,userNodeNumber,meshNodeNumber,err,error,*999)
    numberOfDerivatives=meshNodes%nodes(meshNodeNumber)%numberOfDerivatives
#ifdef WITH_PRECHECKS    
    IF(SIZE(derivatives,1)<numberOfDerivatives) THEN
      localError="The size of the supplied derivatives array of "// &
        & TRIM(NumberToVString(SIZE(derivatives,1),"*",err,error))// &
        & " is too small. The size should be >= "// &
        & TRIM(NumberToVString(numberOfDerivatives,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    DO derivativeIdx=1,numberOfDerivatives
      derivatives(derivativeIdx)=meshNodes%nodes(meshNodeNumber)%derivatives(derivativeIdx)%globalDerivativeIndex
    ENDDO !derivativeIdx
    
    EXITS("MeshNodes_NodeDerivativesGet")
    RETURN
999 ERRORSEXITS("MeshNodes_NodeDerivativesGet",err,error)    
    RETURN 1
   
  END SUBROUTINE MeshNodes_NodeDerivativesGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of versions for a derivative of a node in mesh nodes
  SUBROUTINE MeshNodes_NodeNumberOfVersionsGet(meshNodes,derivativeNumber,userNodeNumber,numberOfVersions,err,error,*)

    !Argument variables
    TYPE(MeshNodesType), POINTER :: meshNodes !<A pointer to the mesh nodes containing the node to get the number of versions for
    INTEGER(INTG), INTENT(IN) :: derivativeNumber !<The derivative number of the node to get the number of versions for
    INTEGER(INTG), INTENT(IN) :: userNodeNumber !<The user number of the node to get the number of versions for
    INTEGER(INTG), INTENT(OUT) :: numberOfVersions !<On return, the number of versions for the specified derivative of the specified node
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: meshNodeNumber
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("MeshNodes_NodeNumberOfVersionsGet",err,error,*999)

#ifdef WITH_PRECECKS    
    IF(.NOT.ASSOCIATED(meshNodes)) CALL FlagError("Mesh nodes is not associated.",err,error,*999)
#endif    
    CALL MeshNodes_MeshNodeNumberGet(meshNodes,userNodeNumber,meshNodeNumber,err,error,*999)
#ifdef WITH_PRECHECKS    
     IF(derivativeNumber<1.OR.derivativeNumber>meshNodes%nodes(meshNodeNumber)%numberOfDerivatives) THEN
      localError="The specified derivative index of "// &
        & TRIM(NumberToVString(derivativeNumber,"*",err,error))// &
        & " is invalid. The derivative index must be >= 1 and <= "// &
        & TRIM(NumberToVString(meshNodes%nodes(meshNodeNumber)%numberOfDerivatives,"*",err,error))// &
        & " for user node number "//TRIM(NumberToVString(userNodeNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    numberOfVersions=meshNodes%nodes(meshNodeNumber)%derivatives(derivativeNumber)%numberOfVersions
    
    EXITS("MeshNodes_NodeNumberOfVersionsGet")
    RETURN
999 ERRORSEXITS("MeshNodes_NodeNumberOfVersionsGet",err,error)    
    RETURN 1
   
  END SUBROUTINE MeshNodes_NodeNumberOfVersionsGet
  
  !
  !================================================================================================================================
  !

  !>Returns the number of nodes for a node in a mesh
  SUBROUTINE MeshNodes_NumberOfNodesGet(meshNodes,numberOfNodes,err,error,*)

    !Argument variables
    TYPE(MeshNodesType), POINTER :: meshNodes !<A pointer to the mesh nodes containing the nodes to get the number of nodes for
    INTEGER(INTG), INTENT(OUT) :: numberOfNodes !<On return, the number of nodes in the mesh.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
  
    ENTERS("MeshNodes_NumberOfNodesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(meshNodes)) CALL FlagError("Mesh nodes is not associated.",err,error,*999)
#endif
    
    numberOfNodes=meshNodes%numberOfNodes
    
    EXITS("MeshNodes_NumberOfNodesGet")
    RETURN
999 ERRORSEXITS("MeshNodes_NumberOfNodesGet",err,error)    
    RETURN 1
   
  END SUBROUTINE MeshNodes_NumberOfNodesGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the user node number for a node in mesh nodes
  SUBROUTINE MeshNodes_UserNodeNumberGet(meshNodes,meshNodeNumber,userNodeNumber,err,error,*)

    !Argument variables
    TYPE(MeshNodesType), POINTER :: meshNodes !<A pointer to the mesh nodes containing the node to get the user number for
    INTEGER(INTG), INTENT(IN) :: meshNodeNumber !<The mesh number of the node to get the user number for
    INTEGER(INTG), INTENT(OUT) :: userNodeNumber !<On return, the user node number of the specified node
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif
    ENTERS("MeshNodes_UserNodeNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(meshNodes)) CALL FlagError("Mesh nodes is not associated.",err,error,*999)
    IF(meshNodeNumber<1.OR.meshNodeNumber>meshNodes%numberOfNodes) THEN
      localError="The specified mesh node number of "//TRIM(NumberToVString(meshNodeNumber,"*",err,error))// &
        & " is invalid. The mesh node number must be >= 1 and <= "// &
        & TRIM(NumberToVString(meshNodes%numberOfNodes,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(meshNodes%nodes)) CALL FlagError("Mesh nodes nodes is not allocated.",err,error,*999)
#endif    
    
    userNodeNumber=meshNodes%nodes(meshNodeNumber)%userNumber
        
    EXITS("MeshNodes_UserNodeNumberGet")
    RETURN
999 ERRORSEXITS("MeshNodes_UserNodeNumberGet",err,error)    
    RETURN 1
   
  END SUBROUTINE MeshNodes_UserNodeNumberGet
  
  !
  !================================================================================================================================
  !

  !>Returns the mesh data points for a mesh topology
  SUBROUTINE MeshTopology_MeshDataPointsGet(meshTopology,meshDataPoints,err,error,*)

    !Argument variables
    TYPE(MeshTopologyType), POINTER :: meshTopology !<A pointer to the mesh topology to get the mesh data points for
    TYPE(MeshDataPointsType), POINTER :: meshDataPoints !<On return, the mesh data points for the mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshTopology_MeshDataPointsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(ASSOCIATED(meshDataPoints)) CALL FlagError("Mesh data points is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(meshTopology)) CALL FlagError("Mesh topology is not associated.",err,error,*999)
#endif
    
    meshDataPoints=>meshTopology%dataPoints

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(meshDataPoints)) CALL FlagError("Mesh data points is not associated for the mesh topology.",err,error,*999)
#endif    

    EXITS("MeshTopology_MeshDataPointsGet")
    RETURN
999 NULLIFY(meshDataPoints)
998 ERRORSEXITS("MeshTopology_MeshDataPointsGet",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_MeshDataPointsGet

  !
  !================================================================================================================================
  !

  !>Returns the mesh dofs for a mesh topology
  SUBROUTINE MeshTopology_MeshDofsGet(meshTopology,meshDofs,err,error,*)

    !Argument variables
    TYPE(MeshTopologyType), POINTER :: meshTopology !<A pointer to the mesh topology to get the mesh dofs for
    TYPE(MeshDofsType), POINTER :: meshDofs !<On return, the mesh dofs for the mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshTopology_MeshDofsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(ASSOCIATED(meshDofs)) CALL FlagError("Mesh dofs is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(meshTopology)) CALL FlagError("Mesh topology is not associated.",err,error,*999)
#endif
    
    meshDofs=>meshTopology%dofs

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(meshDofs)) CALL FlagError("Mesh dofs is not associated for the mesh topology.",err,error,*999)
#endif    

    EXITS("MeshTopology_MeshDofsGet")
    RETURN
999 NULLIFY(meshDofs)
998 ERRORSEXITS("MeshTopology_MeshDofsGet",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_MeshDofsGet

  !
  !================================================================================================================================
  !

  !>Returns the mesh elements for a mesh topology
  SUBROUTINE MeshTopology_MeshElementsGet(meshTopology,meshElements,err,error,*)

    !Argument variables
    TYPE(MeshTopologyType), POINTER :: meshTopology !<A pointer to the mesh topology to get the mesh elements for
    TYPE(MeshElementsType), POINTER :: meshElements !<On return, the mesh elements for the mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshTopology_MeshElementsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(ASSOCIATED(meshElements)) CALL FlagError("Mesh elements is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(meshTopology)) CALL FlagError("Mesh topology is not associated.",err,error,*998)
#endif
    
    meshElements=>meshTopology%elements

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(meshElements)) CALL FlagError("Mesh elments is not associated for the mesh topology.",err,error,*999)
#endif
    
    EXITS("MeshTopology_MeshElementsGet")
    RETURN
999 NULLIFY(meshElements)
998 ERRORSEXITS("MeshTopology_MeshElementsGet",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_MeshElementsGet

  !
  !================================================================================================================================
  !

  !>Returns the mesh for a mesh topology
  SUBROUTINE MeshTopology_MeshGet(meshTopology,mesh,err,error,*)

    !Argument variables
    TYPE(MeshTopologyType), POINTER :: meshTopology !<A pointer to the mesh topology to get the mesh for
    TYPE(MeshType), POINTER :: mesh !<On return, the mesh for the mesh topology. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshTopology_MeshGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(ASSOCIATED(mesh)) CALL FlagError("Mesh is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(meshTopology)) CALL FlagError("Mesh topology is not associated.",err,error,*998)
#endif
    
    mesh=>meshTopology%mesh

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Mesh is not associated for the mesh topology.",err,error,*999)
#endif
    
    EXITS("MeshTopology_MeshGet")
    RETURN
999 NULLIFY(mesh)
998 ERRORSEXITS("MeshTopology_MeshGet",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_MeshGet

  !
  !================================================================================================================================
  !

  !>Returns the mesh nodes for a mesh topology
  SUBROUTINE MeshTopology_MeshNodesGet(meshTopology,meshNodes,err,error,*)

    !Argument variables
    TYPE(MeshTopologyType), POINTER :: meshTopology !<A pointer to the mesh topology to get the mesh nodes
    TYPE(MeshNodesType), POINTER :: meshNodes !<On return, the mesh nodes for the mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshTopology_MeshNodesGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    !Check input arguments
    IF(ASSOCIATED(meshNodes)) CALL FlagError("Mesh nodes is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(meshTopology)) CALL FlagError("Mesh topology is not associated.",err,error,*998)
#endif
    
    meshNodes=>meshTopology%nodes

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(meshNodes)) CALL FlagError("Mesh nodes is not associated for the mesh topology.",err,error,*999)
#endif
    
    EXITS("MeshTopology_MeshNodesGet")
    RETURN
999 NULLIFY(meshNodes)
998 ERRORSEXITS("MeshTopology_MeshNodesGet",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_MeshNodesGet

  !
  !================================================================================================================================
  !
  
END MODULE MeshAccessRoutines
