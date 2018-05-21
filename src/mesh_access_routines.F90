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
  USE Trees
  USE Types
  
#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  INTERFACE Mesh_UserNumberFind
    MODULE PROCEDURE Mesh_UserNumberFindInterface
    MODULE PROCEDURE Mesh_UserNumberFindRegion
  END INTERFACE Mesh_UserNumberFind

  INTERFACE MESH_USER_NUMBER_FIND
    MODULE PROCEDURE Mesh_UserNumberFindInterface
    MODULE PROCEDURE Mesh_UserNumberFindRegion
  END INTERFACE MESH_USER_NUMBER_FIND

  PUBLIC Mesh_AssertIsFinished,Mesh_AssertNotFinished

  PUBLIC Mesh_DecompositionGet

  PUBLIC Mesh_MeshElementsGet

  PUBLIC Mesh_MeshNodesGet

  PUBLIC Mesh_NodesGet

  PUBLIC Mesh_NumberOfComponentsGet

  PUBLIC Mesh_NumberOfElementsGet

  PUBLIC Mesh_RegionGet

  PUBLIC Mesh_MeshTopologyGet

  PUBLIC Mesh_UserNumberFind

  PUBLIC MESH_USER_NUMBER_FIND

  PUBLIC Mesh_UserNumberGet

  PUBLIC MeshElements_AssertIsFinished,MeshElements_AssertNotFinished

  PUBLIC MeshElements_BasisGet

  PUBLIC MeshElements_ElementCheckExists

  PUBLIC MeshElements_GlobalElementNumberGet
    
  PUBLIC MeshElements_MeshElementGet
  
  PUBLIC MeshNodes_GlobalNodeNumberGet
  
  PUBLIC MeshNodes_MeshNodeGet

  PUBLIC MeshNodes_MeshNodeNumberGet
  
  PUBLIC MeshNodes_NodeCheckExists    

  PUBLIC MeshTopology_ElementCheckExists
  
  PUBLIC MeshTopology_MeshElementsGet

  PUBLIC MeshTopology_NodeCheckExists

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

    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Mesh is not associated.",err,error,*999)

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

    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Mesh is not associated.",err,error,*999)

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
    
    IF(ASSOCIATED(nodes)) CALL FlagError("Nodes is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Mesh is not associated.",err,error,*999)
    
    region=>mesh%region
    IF(ASSOCIATED(region)) THEN
      nodes=>region%nodes
    ELSE
      interface=>mesh%interface
      IF(ASSOCIATED(interface)) THEN
        nodes=>interface%nodes
      ELSE
        localError="Mesh number "//TRIM(NumberToVString(mesh%userNumber,"*",err,error))// &
          & " does not have an associated region or interface."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDIF
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

  !>Returns the mesh topology for a mesh and mesh component number.
  SUBROUTINE Mesh_MeshTopologyGet(mesh,meshComponentNumber,meshTopology,err,error,*)

    !Argument variables
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh to get the mesh topology for
    INTEGER(INTG), INTENT(IN) :: meshComponentNumber !<The mesh component number of the topology to get
    TYPE(MeshTopologyType), POINTER :: meshTopology !<On return, the mesh topology for the mesh and mesh component number. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Mesh_MeshTopologyGet",err,error,*998)

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

    meshTopology=>mesh%topology(meshComponentNumber)%ptr
    IF(.NOT.ASSOCIATED(meshTopology)) THEN    
      localError="The mesh topology is not associated for mesh component number "// &
        & TRIM(NumberToVString(meshComponentNumber,"*",err,error))//" of mesh number "// &
        & TRIM(NumberToVString(mesh%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

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

    !Check input arguments
    IF(.NOT.ASSOCIATED(meshes)) CALL FlagError("Meshes is not associated",err,error,*999)
    IF(ASSOCIATED(mesh)) CALL FlagError("Mesh is already associated.",err,error,*999)
    
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

    IF(.NOT.ASSOCIATED(meshElements)) CALL FlagError("Mesh elements is not associated.",err,error,*999)

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

    IF(.NOT.ASSOCIATED(meshElements)) CALL FlagError("Mesh elements is not associated.",err,error,*999)

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
    TYPE(VARYING_STRING) :: localError

    ENTERS("MeshElements_BasisGet",err,error,*998)

    IF(ASSOCIATED(basis)) CALL FlagError("Basis is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(meshElements)) CALL FlagError("Mesh elements is not associated.",err,error,*999)
    IF(globalElementNumber<=0.OR.globalElementNumber>meshElements%numberOfElements) THEN
      localError="The global element number of "//TRIM(NumberToVString(globalElementNumber,"*",err,error))// &
        & " is invalid. The global number must be >= 1 and <= "// &
        & TRIM(NumberToVString(meshElements%numberOfElements,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(meshElements%elements)) CALL FlagError("Mesh elements elements has not been allocated.",err,error,*999)
    
    basis=>meshElements%elements(globalElementNumber)%basis
    IF(.NOT.ASSOCIATED(basis)) THEN
      localError="The basis for global element number "//TRIM(NumberToVString(globalElementNumber,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
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
    TYPE(TREE_NODE_TYPE), POINTER :: treeNode
    
    ENTERS("MeshElements_ElementCheckExists",err,error,*999)

    elementExists=.FALSE.
    globalElementNumber=0
    IF(.NOT.ASSOCIATED(meshElements)) CALL FlagError("Mesh elements is not associated.",err,error,*999)
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

    IF(.NOT.ASSOCIATED(meshElements)) CALL FlagError("Mesh elements is not associated.",err,error,*999)
    
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
    TYPE(VARYING_STRING) :: localError

    ENTERS("MeshElements_MeshElementGet",err,error,*998)

    IF(ASSOCIATED(meshElement)) CALL FlagError("Mesh element is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(meshElements)) CALL FlagError("Mesh elements is not associated.",err,error,*999)
    IF(globalElementNumber<=0.OR.globalElementNumber>meshElements%numberOfElements) THEN
      localError="The global element number of "//TRIM(NumberToVString(globalElementNumber,"*",err,error))// &
        & " is invalid. The global number must be >= 1 and <= "// &
        & TRIM(NumberToVString(meshElements%numberOfElements,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(meshElements%elements)) CALL FlagError("Mesh elements elements has not been allocated.",err,error,*999)
    
    meshElement=>meshElements%elements(globalElementNumber)
    IF(.NOT.ASSOCIATED(meshElement)) THEN
      localError="The mesh element for global element number "//TRIM(NumberToVString(globalElementNumber,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("MeshElements_MeshElementGet")
    RETURN
999 NULLIFY(meshElement)
998 ERRORSEXITS("MeshElements_MeshElementGet",err,error)
    RETURN 1
    
  END SUBROUTINE MeshElements_MeshElementGet

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
    TYPE(VARYING_STRING) :: localError

    ENTERS("MeshNodes_MeshNodeGet",err,error,*998)

    IF(ASSOCIATED(meshNode)) CALL FlagError("Mesh node is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(meshNodes)) CALL FlagError("Mesh nodes is not associated.",err,error,*999)
    IF(meshNodeNumber<=0.OR.meshNodeNumber>meshNodes%numberOfNodes) THEN
      localError="The mesh node number of "//TRIM(NumberToVString(meshNodeNumber,"*",err,error))// &
        & " is invalid. The mesh number must be >= 1 and <= "// &
        & TRIM(NumberToVString(meshNodes%numberOfNodes,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(meshNodes%nodes)) CALL FlagError("Mesh nodes nodes has not been allocated.",err,error,*999)
    
    meshNode=>meshNodes%nodes(meshNodeNumber)
    IF(.NOT.ASSOCIATED(meshNode)) THEN
      localError="The mesh node for mesh node number "//TRIM(NumberToVString(meshNodeNumber,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
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

    IF(.NOT.ASSOCIATED(meshNodes)) CALL FlagError("Mesh nodes is not associated.",err,error,*999)
    
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

    IF(.NOT.ASSOCIATED(meshNodes)) CALL FlagError("Mesh nodes is not associated.",err,error,*999)
    
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
    TYPE(TREE_NODE_TYPE), POINTER :: treeNode
    
    ENTERS("MeshNodes_NodeCheckExists",err,error,*999)

    nodeExists=.FALSE.
    meshNodeNumber=0
    IF(.NOT.ASSOCIATED(meshNodes)) CALL FlagError("Mesh nodes is not associated.",err,error,*999)

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
        CALL Tree_NodeValuesGet(meshNodes%nodesTree,treeNode,meshNodeNumber,err,error,*999)
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

  !>Returns the mesh for a mesh topology
  SUBROUTINE MeshTopology_MeshGet(meshTopology,mesh,err,error,*)

    !Argument variables
    TYPE(MeshTopologyType), POINTER :: meshTopology !<A pointer to the mesh topology to get the mesh for
    TYPE(MeshType), POINTER :: mesh !<On return, the mesh for the mesh topology. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshTopology_MeshGet",err,error,*998)

    !Check input arguments
    IF(ASSOCIATED(mesh)) CALL FlagError("Mesh is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(meshTopology)) CALL FlagError("Mesh topology is not associated.",err,error,*998)

    mesh=>meshTopology%mesh
    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Mesh is not associated for the mesh topology.",err,error,*999)

    EXITS("MeshTopology_MeshGet")
    RETURN
999 NULLIFY(mesh)
998 ERRORSEXITS("MeshTopology_MeshGet",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_MeshGet

  !
  !================================================================================================================================
  !

  !>Checks that a user element number exists in a mesh elements. 
  SUBROUTINE MeshTopology_ElementCheckExists(meshTopology,userElementNumber,elementExists,globalElementNumber,err,error,*)

    !Argument variables
    TYPE(MeshTopologyType), POINTER :: meshTopology !<A pointer to the mesh topology to check the element exists on
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number to check if it exists
    LOGICAL, INTENT(OUT) :: elementExists !<On exit, is .TRUE. if the element user number exists in the mesh elements, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: globalElementNumber !<On exit, if the element exists the global number corresponding to the user element number. If the element does not exist then global number will be 0.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(MeshElementsType), POINTER :: meshElements
    TYPE(TREE_NODE_TYPE), POINTER :: treeNode
    
    ENTERS("MeshTopology_ElementCheckExists",err,error,*999)

    elementExists=.FALSE.
    globalElementNumber=0
    IF(.NOT.ASSOCIATED(meshTopology)) CALL FlagError("Mesh topology is not associated.",err,error,*999)
    
    NULLIFY(meshElements)
    CALL MeshTopology_MeshElementsGet(meshTopology,meshElements,err,error,*999)
    NULLIFY(treeNode)
    CALL Tree_Search(meshElements%elementsTree,userElementNumber,treeNode,err,error,*999)
    IF(ASSOCIATED(treeNode)) THEN
      CALL Tree_NodeValueGet(meshElements%elementsTree,treeNode,globalElementNumber,err,error,*999)
      elementExists=.TRUE.      
    ENDIF
    
    EXITS("MeshTopology_ElementCheckExists")
    RETURN
999 ERRORSEXITS("MeshTopology_ElementCheckExists",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_ElementCheckExists
  
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

    !Check input arguments
    IF(ASSOCIATED(meshElements)) CALL FlagError("Mesh elements is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(meshTopology)) CALL FlagError("Mesh topology is not associated.",err,error,*998)

    meshElements=>meshTopology%elements
    IF(.NOT.ASSOCIATED(meshElements)) CALL FlagError("Mesh elments is not associated for the mesh topology.",err,error,*999)

    EXITS("MeshTopology_MeshElementsGet")
    RETURN
999 NULLIFY(meshElements)
998 ERRORSEXITS("MeshTopology_MeshElementsGet",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_MeshElementsGet

  !
  !================================================================================================================================
  !

  !>Checks that a user node number exists in a mesh component. 
  SUBROUTINE MeshTopology_NodeCheckExists(meshTopology,userNodeNumber,nodeExists,globalNodeNumber,meshNodeNumber,err,error,*)

    !Argument variables
    TYPE(MeshTopologyType), POINTER :: meshTopology !<A pointer to the mesh topology to check the node exists on
    INTEGER(INTG), INTENT(IN) :: userNodeNumber !<The user node number to check if it exists
    LOGICAL, INTENT(OUT) :: nodeExists !<On exit, is .TRUE. if the node user number exists in the mesh component, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: globalNodeNumber !<On exit, if the node exists the global number corresponding to the user node number. If the node does not exist then global number will be 0.
    INTEGER(INTG), INTENT(OUT) :: meshNodeNumber !<On exit, if the node exists the mesh number corresponding to the user node number. If the node does not exist then mesh number will be 0.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(MeshNodesType), POINTER :: meshNodes
    
    ENTERS("MeshTopology_NodeCheckExists",err,error,*999)

    nodeExists=.FALSE.
    meshNodeNumber=0
    IF(.NOT.ASSOCIATED(meshTopology)) CALL FlagError("Mesh topology is not associated.",err,error,*999)
    NULLIFY(meshNodes)
    CALL MeshTopology_MeshNodesGet(meshTopology,meshNodes,err,error,*999)
    CALL MeshNodes_NodeCheckExists(meshNodes,userNodeNumber,nodeExists,globalNodeNumber,meshNodeNumber,err,error,*999)
   
    EXITS("MeshTopology_NodeCheckExists")
    RETURN
999 ERRORSEXITS("MeshTopology_NodeCheckExists",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_NodeCheckExists
  
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

    !Check input arguments
    IF(ASSOCIATED(meshNodes)) CALL FlagError("Mesh nodes is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(meshTopology)) CALL FlagError("Mesh topology is not associated.",err,error,*998)

    meshNodes=>meshTopology%nodes
    IF(.NOT.ASSOCIATED(meshNodes)) CALL FlagError("Mesh nodes is not associated for the mesh topology.",err,error,*999)

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
