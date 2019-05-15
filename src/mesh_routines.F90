!> \file
!> \author Chris Bradley
!> \brief This module handles all mesh (node and element) routines.
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

!> This module handles all mesh (node and element) routines.
MODULE MeshRoutines

  USE BaseRoutines
  USE BasisRoutines
  USE BasisAccessRoutines
  USE Constants
  USE ContextAccessRoutines
  USE COORDINATE_ROUTINES
  USE CoordinateSystemAccessRoutines
  USE DataProjectionAccessRoutines
  USE DecompositionRoutines
  USE Kinds
  USE INPUT_OUTPUT
  USE InterfaceAccessRoutines
  USE ISO_VARYING_STRING
  USE Lists
  USE MeshAccessRoutines
#ifndef NOMPIMOD
  USE MPI
#endif
  USE NodeAccessRoutines
  USE RegionAccessRoutines
  USE Strings
  USE Trees
  USE Types

#include "macros.h"  

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces
  
  !>Starts the process of creating a mesh
  INTERFACE Mesh_CreateStart
    MODULE PROCEDURE Mesh_CreateStartInterface
    MODULE PROCEDURE Mesh_CreateStartRegion
  END INTERFACE Mesh_CreateStart

  !>Initialises the meshes for a region or interface.
  INTERFACE Meshes_Initialise
    MODULE PROCEDURE Meshes_InitialiseInterface
    MODULE PROCEDURE Meshes_InitialiseRegion
  END INTERFACE Meshes_Initialise

  PUBLIC Mesh_CreateStart,Mesh_CreateFinish

  PUBLIC Mesh_Destroy
  
  PUBLIC Mesh_NumberOfComponentsSet

  PUBLIC Mesh_NumberOfElementsGet,Mesh_NumberOfElementsSet

  PUBLIC Mesh_SurroundingElementsCalculateSet

  PUBLIC MeshElements_CreateStart,MeshElements_CreateFinish

  PUBLIC MeshElements_Destroy

  PUBLIC MeshElements_ElementBasisGet,MeshElements_ElementBasisSet

  PUBLIC MeshElements_AdjacentElementGet

  PUBLIC MeshElements_ElementNodesGet

  PUBLIC MeshElements_ElementNodesSet

  PUBLIC MeshElements_ElementNodeVersionSet,MeshElements_UserNodeVersionSet

  PUBLIC MeshElements_ElementOnBoundaryGet

  PUBLIC MeshElements_ElementUserNumberGet,MeshElements_ElementUserNumberSet
  
  PUBLIC MeshElements_ElementsUserNumbersAllSet

  PUBLIC MeshTopology_DataPointsCalculateProjection

  PUBLIC MeshTopology_NodesDestroy
  
  PUBLIC MESH_EMBEDDING_CREATE,MESH_EMBEDDING_SET_CHILD_NODE_POSITION

  PUBLIC MESH_EMBEDDING_SET_GAUSS_POINT_DATA

  PUBLIC Meshes_Initialise,Meshes_Finalise

CONTAINS

  !
  !================================================================================================================================
  !

  !>Finialises the mesh adjacent elements information and deallocates all memory
  SUBROUTINE Mesh_AdjacentElementFinalise(meshAdjacentElement,err,error,*)
    
    !Argument variables
    TYPE(MeshAdjacentElementType) :: meshAdjacentElement !<The mesh adjacent element to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Mesh_AdjacentElementFinalise",err,error,*999)

    meshAdjacentElement%numberOfAdjacentElements=0
    IF(ALLOCATED(meshAdjacentElement%adjacentElements)) DEALLOCATE(meshAdjacentElement%adjacentElements)
       
    EXITS("Mesh_AdjacentElementFinalise")
    RETURN
999 ERRORSEXITS("Mesh_AdjacentElementFinalise",err,error)    
    RETURN 1
   
  END SUBROUTINE Mesh_AdjacentElementFinalise

  !
  !================================================================================================================================
  !
  !>Initalises the mesh adjacent elements information.
  SUBROUTINE Mesh_AdjacentElementInitialise(meshAdjacentElement,err,error,*)
    
    !Argument variables
    TYPE(MeshAdjacentElementType) :: meshAdjacentElement !<The mesh adjacent element to initialise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Mesh_AdjacentElementInitialise",err,error,*999)

    meshAdjacentElement%numberOfAdjacentElements=0
       
    EXITS("Mesh_AdjacentElementInitialise")
    RETURN
999 ERRORSEXITS("Mesh_AdjacentElementInitialise",err,error)    
    RETURN 1
   
  END SUBROUTINE Mesh_AdjacentElementInitialise

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a mesh. \see OpenCMISS::Iron::cmfe_Mesh_CreateFinish
  SUBROUTINE Mesh_CreateFinish(mesh,err,error,*)

    !Argument variables
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh to finish creating
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx
    TYPE(MeshTopologyType), POINTER :: meshTopology
    TYPE(MeshElementsType), POINTER :: meshElements

    ENTERS("Mesh_CreateFinish",err,error,*999)

    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Mesh is not associated.",err,error,*999)
    
    !Check that the mesh component elements have been finished
    DO componentIdx=1,mesh%numberOfComponents
      NULLIFY(meshTopology)
      CALL Mesh_MeshTopologyGet(mesh,componentIdx,meshTopology,err,error,*999)
      NULLIFY(meshElements)
      CALL MeshTopology_MeshElementsGet(meshTopology,meshElements,err,error,*999)
      CALL MeshElements_AssertIsFinished(meshElements,err,error,*999)
    ENDDO !componentIdx
    mesh%meshFinished=.TRUE.
    !Calulcate the mesh topology
    DO componentIdx=1,mesh%numberOfComponents
      NULLIFY(meshTopology)
      CALL Mesh_MeshTopologyGet(mesh,componentIdx,meshTopology,err,error,*999)
      CALL MeshTopology_Calculate(meshTopology,err,error,*999)
    ENDDO !componentIdx
    
    IF(diagnostics1) THEN
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Mesh user number       = ",mesh%userNumber,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Global number        = ",mesh%globalNumber,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of dimensions = ",mesh%numberOfDimensions,err,error,*999)
    ENDIF
    
    EXITS("Mesh_CreateFinish")
    RETURN
999 ERRORSEXITS("Mesh_CreateFinish",err,error)    
    RETURN 1
   
  END SUBROUTINE Mesh_CreateFinish
        
  !
  !================================================================================================================================
  !

  !>Starts the process of creating a mesh 
  SUBROUTINE Mesh_CreateStartGeneric(meshes,userNumber,numberOfDimensions,mesh,err,error,*)
    
    !Argument variables
    TYPE(MeshesType), POINTER :: meshes !<The pointer to the meshes
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the mesh to create
    INTEGER(INTG), INTENT(IN) :: numberOfDimensions !<The number of dimensions in the mesh.
    TYPE(MeshType), POINTER :: mesh !<On return, a pointer to the mesh. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,meshIdx
    TYPE(MeshType), POINTER :: newMesh
    TYPE(MeshPtrType), ALLOCATABLE :: newMeshes(:)
    TYPE(VARYING_STRING) :: dummyError

    NULLIFY(newMesh)

    ENTERS("Mesh_CreateStartGeneric",err,error,*997)

    IF(ASSOCIATED(mesh)) CALL FlagError("Mesh is already associated.",err,error,*997)    
    IF(.NOT.ASSOCIATED(meshes)) CALL FlagError("Meshes is not associated.",err,error,*997)
    
    CALL Mesh_Initialise(newMesh,err,error,*999)
    !Set default mesh values
    newMesh%userNumber=userNumber
    newMesh%globalNumber=meshes%numberOfMeshes+1
    newMesh%meshes=>meshes
    newMesh%numberOfDimensions=numberOfDimensions
    newMesh%numberOfComponents=1
    newMesh%surroundingElementsCalculate=.TRUE. !default true
    !Initialise mesh topology and decompositions
    CALL Mesh_TopologiesInitialise(newMesh,err,error,*999)
    CALL Mesh_DecompositionsInitialise(newMesh,err,error,*999)
    !Add new mesh into list of meshes 
    ALLOCATE(newMeshes(meshes%numberOfMeshes+1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new meshes",err,error,*999)
    DO meshIdx=1,meshes%numberOfMeshes
      newMeshes(meshIdx)%ptr=>meshes%meshes(meshIdx)%ptr
    ENDDO !meshIdx
    newMeshes(meshes%numberOfMeshes+1)%ptr=>newMesh
    CALL MOVE_ALLOC(newMeshes,meshes%meshes)
    meshes%numberOfMeshes=meshes%numberOfMeshes+1
    !Return the pointer to the new mesh
    mesh=>newMesh
      
    EXITS("Mesh_CreateStartGeneric")
    RETURN
999 CALL Mesh_Finalise(newMesh,dummyErr,dummyError,*998)
998 IF(ALLOCATED(newMeshes)) DEALLOCATE(newMeshes)
    NULLIFY(mesh)    
997 ERRORSEXITS("Mesh_CreateStartGeneric",err,error)    
    RETURN 1
   
  END SUBROUTINE Mesh_CreateStartGeneric

  !
  !================================================================================================================================
  !

  !>Starts the process of creating a mesh defined by a user number with the specified numberOfDimensions in an interface. \see OpenCMISS::Iron::cmfe_Mesh_CreateStart
  !>Default values set for the meshes's attributes are:
  !>- NUMBER_OF_COMPONENTS: 1
  SUBROUTINE Mesh_CreateStartInterface(userNumber,interface,numberOfMeshDimensions,mesh,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the mesh to create
    TYPE(InterfaceType), POINTER :: interface !<A pointer to the interface to create the mesh on
    INTEGER(INTG), INTENT(IN) :: numberOfMeshDimensions !<The number of dimensions in the mesh.
    TYPE(MeshType), POINTER :: mesh !<On exit, a pointer to the created mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfInterfaceDimensions
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem
    TYPE(MeshesType), POINTER :: meshes
    TYPE(VARYING_STRING) :: localError

    ENTERS("Mesh_CreateStartInterface",err,error,*999)

    IF(ASSOCIATED(mesh)) CALL FlagError("Mesh is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(interface)) CALL FlagError("Interface is not associated.",err,error,*999)
    IF(numberOfMeshDimensions<=0) THEN
      localError="The specified number of mesh dimensions of "//TRIM(NumberToVString(numberOfMeshDimensions,"*",err,error))// &
        & " is invalid. The number of dimensions must be > 0."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    NULLIFY(coordinateSystem)
    CALL Interface_CoordinateSystemGet(interface,coordinateSystem,err,error,*999)
    CALL CoordinateSystem_DimensionGet(coordinateSystem,numberOfInterfaceDimensions,err,error,*999)
    IF(numberOfMeshDimensions>numberOfInterfaceDimensions) THEN
      localError="The specified number of mesh dimensions of "//TRIM(NumberToVString(numberOfMeshDimensions,"*",err,error))// &
        & " is invalid for the number of interface dimensions of "// &
        & TRIM(NumberToVString(numberOfInterfaceDimensions,"*",err,error))// &
        & ". The number of mesh dimensions must be <= the number of interface dimensions."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    NULLIFY(meshes)
    CALL Interface_MeshesGet(interface,meshes,err,error,*999)
    CALL Mesh_UserNumberFind(userNumber,interface,mesh,err,error,*999)
    IF(ASSOCIATED(mesh)) THEN
      localError="Mesh number "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " has already been created on interface number "// &
        & TRIM(NumberToVString(INTERFACE%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    CALL Mesh_CreateStartGeneric(meshes,userNumber,numberOfMeshDimensions,mesh,err,error,*999)
    mesh%interface=>interface
    
    EXITS("Mesh_CreateStartInterface")
    RETURN
999 ERRORSEXITS("Mesh_CreateStartInterface",err,error)    
    RETURN 1
   
  END SUBROUTINE Mesh_CreateStartInterface

  !
  !================================================================================================================================
  !

  !>Starts the process of creating a mesh defined by a user number with the specified numberOfDimensions in the region identified by region. \see OpenCMISS::Iron::cmfe_Mesh_CreateStart
  !>Default values set for the mesh's attributes are:
  !>- NUMBER_OF_COMPONENTS: 1
  SUBROUTINE Mesh_CreateStartRegion(userNumber,region,numberOfMeshDimensions,mesh,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the mesh to create
    TYPE(RegionType), POINTER :: region !<A pointer to the region to create the mesh on
    INTEGER(INTG), INTENT(IN) :: numberOfMeshDimensions !<The number of dimensions in the mesh.
    TYPE(MeshType), POINTER :: mesh !<On exit, a pointer to the created mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfRegionDimensions
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem
    TYPE(MeshesType), POINTER :: meshes
    TYPE(VARYING_STRING) :: localError

    ENTERS("Mesh_CreateStartRegion",err,error,*999)

    IF(ASSOCIATED(mesh)) CALL FlagError("Mesh is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated.",err,error,*999)
    IF(numberOfMeshDimensions<=0) THEN
      localError="The specified number of mesh dimensions of "//TRIM(NumberToVString(numberOfMeshDimensions,"*",err,error))// &
        & " is invalid. The number of dimensions must be > 0."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    NULLIFY(coordinateSystem)
    CALL Region_CoordinateSystemGet(region,coordinateSystem,err,error,*999)
    CALL CoordinateSystem_DimensionGet(coordinateSystem,numberOfRegionDimensions,err,error,*999)
    IF(numberOfMeshDimensions>numberOfRegionDimensions) THEN
      localError="The specified number of mesh dimensions of "//TRIM(NumberToVString(numberOfMeshDimensions,"*",err,error))// &
        & " is invalid for the number of region dimensions of "// &
        & TRIM(NumberToVString(numberOfRegionDimensions,"*",err,error))// &
        & ". The number of mesh dimensions must be <= the number of region dimensions."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    NULLIFY(meshes)
    CALL Region_MeshesGet(region,meshes,err,error,*999)
    CALL Mesh_UserNumberFind(userNumber,region,mesh,err,error,*999)
    IF(ASSOCIATED(mesh)) THEN
      localError="Mesh number "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " has already been created on region number "// &
        & TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    CALL Mesh_CreateStartGeneric(meshes,userNumber,numberOfMeshDimensions,mesh,err,error,*999)
    mesh%region=>region
    
    EXITS("Mesh_CreateStartRegion")
    RETURN
999 ERRORSEXITS("Mesh_CreateStartRegion",err,error)    
    RETURN 1
   
  END SUBROUTINE Mesh_CreateStartRegion

  !
  !================================================================================================================================
  !

  !>Finalises the decompositions and deallocate all memory.
  SUBROUTINE Mesh_DecompositionsFinalise(decompositions,err,error,*)

    !Argument variables
    TYPE(DecompositionsType), POINTER :: decompositions !<A pointer to the decompositions to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Mesh_DecompositionsFinalise",err,error,*999)

    IF(ASSOCIATED(decompositions)) THEN
      DO WHILE(decompositions%numberOfDecompositions>0)
        CALL Decomposition_Destroy(decompositions%decompositions(1)%ptr,err,error,*999)
      ENDDO
      DEALLOCATE(decompositions)
    ENDIF
    
    EXITS("Mesh_DecompositionsFinalise")
    RETURN
999 ERRORSEXITS("Mesh_DecompositionsFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE Mesh_DecompositionsFinalise
  
  !
  !================================================================================================================================
  !

  !>Initialises the domain decompositions for a given mesh.
  SUBROUTINE Mesh_DecompositionsInitialise(mesh,err,error,*)

    !Argument variables
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh to initialise the decompositions for

    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("Mesh_DecompositionsInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Mesh is not associated.",err,error,*998)
    IF(ASSOCIATED(mesh%decompositions)) CALL FlagError("Mesh already has decompositions associated.",err,error,*998)
    
    ALLOCATE(mesh%decompositions,STAT=err)
    IF(err/=0) CALL FlagError("Mesh decompositions could not be allocated.",err,error,*999)
    mesh%decompositions%numberOfDecompositions=0
    mesh%decompositions%mesh=>mesh
     
    EXITS("Mesh_DecompositionsInitialise")
    RETURN
999 CALL Mesh_DecompositionsFinalise(mesh%decompositions,dummyErr,dummyError,*999)
998 ERRORSEXITS("Mesh_DecompositionsInitialise",err,error)    
    RETURN 1
    
  END SUBROUTINE Mesh_DecompositionsInitialise
  
  !
  !================================================================================================================================
  !

  !>Destroys the mesh and deallocates all memory. \see OpenCMISS::Iron::cmfe_Mesh_Destroy
  SUBROUTINE Mesh_Destroy(mesh,err,error,*)

    !Argument variables
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh to destroy.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: meshIdx,meshPosition
    TYPE(MeshesType), POINTER :: meshes
    TYPE(MeshPtrType), ALLOCATABLE :: newMeshes(:)

    ENTERS("Mesh_Destroy",err,error,*999)

    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Mesh is not associated.",err,error,*999)
    
    NULLIFY(meshes)
    CALL Mesh_MeshesGet(mesh,meshes,err,error,*999)
    meshPosition=mesh%globalNumber
    
    CALL Mesh_Finalise(mesh,err,error,*999)

    !Remove the mesh from the list of meshes
    IF(meshes%numberOfMeshes>1) THEN
      ALLOCATE(newMeshes(meshes%numberOfMeshes-1),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new meshes.",err,error,*999)
      DO meshIdx=1,meshes%numberOfMeshes
        IF(meshIdx<meshPosition) THEN
          newMeshes(meshIdx)%ptr=>meshes%meshes(meshIdx)%ptr
        ELSE IF(meshIdx>meshPosition) THEN
          meshes%meshes(meshIdx)%ptr%globalNumber=meshes%meshes(meshIdx)%ptr%globalNumber-1
          newMeshes(meshIdx-1)%ptr=>meshes%meshes(meshIdx)%ptr
        ENDIF
      ENDDO !meshIdx
      CALL MOVE_ALLOC(newMeshes,meshes%meshes)
      meshes%numberOfMeshes=meshes%numberOfMeshes-1
    ELSE
      DEALLOCATE(meshes%meshes)
      meshes%numberOfMeshes=0
    ENDIF

    EXITS("Mesh_Destroy")
    RETURN
999 IF(ALLOCATED(newMeshes)) DEALLOCATE(newMeshes)
    ERRORSEXITS("Mesh_Destroy",err,error)
    RETURN 1
    
  END SUBROUTINE Mesh_Destroy

  !
  !================================================================================================================================
  !

  !>Finalises a mesh and deallocates all memory.
  SUBROUTINE Mesh_Finalise(mesh,err,error,*)

    !Argument variables
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Mesh_Finalise",err,error,*999)

    IF(ASSOCIATED(mesh)) THEN
      CALL Mesh_TopologiesFinalise(mesh,err,error,*999)
      CALL Mesh_DecompositionsFinalise(mesh%decompositions,err,error,*999)
      DEALLOCATE(mesh)
    ENDIF
 
    EXITS("Mesh_Finalise")
    RETURN
999 ERRORSEXITS("Mesh_Finalise",err,error)
    RETURN 1
   
  END SUBROUTINE Mesh_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises a mesh.
  SUBROUTINE Mesh_Initialise(mesh,err,error,*)

    !Argument variables
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Mesh_Initialise",err,error,*999)

    IF(ASSOCIATED(mesh)) CALL FlagError("Mesh is already associated.",err,error,*999)
    
    ALLOCATE(mesh,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new mesh.",err,error,*999)
    mesh%userNumber=0
    mesh%globalNumber=0
    mesh%meshFinished=.FALSE.
    NULLIFY(mesh%meshes)
    NULLIFY(mesh%region)
    NULLIFY(mesh%interface)
    NULLIFY(mesh%generatedMesh)
    mesh%numberOfDimensions=0
    mesh%numberOfComponents=0
    mesh%numberOfElements=0
    mesh%MESH_EMBEDDED=.FALSE.
    NULLIFY(mesh%EMBEDDING_MESH)
    mesh%NUMBER_OF_EMBEDDED_MESHES=0
    NULLIFY(mesh%EMBEDDED_MESHES)
    NULLIFY(mesh%decompositions)
  
    EXITS("Mesh_Initialise")
    RETURN
999 ERRORSEXITS("Mesh_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE Mesh_Initialise

  !
  !================================================================================================================================
  !

  !>Changes/sets the number of mesh components for a mesh. \see OpenCMISS::Iron::cmfe_Mesh_NumberOfComponentsSet
  SUBROUTINE Mesh_NumberOfComponentsSet(mesh,numberOfComponents,err,error,*)

    !Argument variables
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh to set the number of components for
    INTEGER(INTG), INTENT(IN) :: numberOfComponents !<The number of components to set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx
    TYPE(MeshTopologyPtrType), ALLOCATABLE :: newTopology(:)
    TYPE(VARYING_STRING) :: localError
   
    ENTERS("Mesh_NumberOfComponentsSet",err,error,*999)

    CALL Mesh_AssertNotFinished(mesh,err,error,*999)
    IF(numberOfComponents<=0) THEN
      localError="The specified number of mesh components of "//TRIM(NumberToVString(numberOfComponents,"*",err,error))// &
        & " is illegal. You must have >0 mesh components."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(numberOfComponents/=mesh%numberOfComponents) THEN
      ALLOCATE(newTopology(numberOfComponents),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new topology.",err,error,*999)
      IF(numberOfComponents<mesh%numberOfComponents) THEN
        DO componentIdx=1,numberOfComponents
          newTopology(componentIdx)%ptr=>mesh%TOPOLOGY(componentIdx)%ptr
        ENDDO !componentIdx
      ELSE !numberOfComponents>mesh%numberOfComponents
        DO componentIdx=1,mesh%numberOfComponents
          newTopology(componentIdx)%ptr=>mesh%TOPOLOGY(componentIdx)%ptr
        ENDDO !componentIdx
        DO componentIdx=mesh%numberOfComponents+1,numberOfComponents
          NULLIFY(newTopology(componentIdx)%ptr)
          CALL MeshTopology_Initialise(newTopology(componentIdx)%ptr,err,error,*999)
          newTopology(componentIdx)%ptr%mesh=>mesh
          newTopology(componentIdx)%ptr%meshComponentNumber=componentIdx
        ENDDO !componentIdx
      ENDIF
      CALL MOVE_ALLOC(newTopology,mesh%topology)
      mesh%numberOfComponents=numberOfComponents
    ENDIF
    
    EXITS("Mesh_NumberOfComponentsSet")
    RETURN
999 IF(ALLOCATED(newTopology)) DEALLOCATE(newTopology)   
    ERRORSEXITS("Mesh_NumberOfComponentsSet",err,error)    
    RETURN 1
   
  END SUBROUTINE Mesh_NumberOfComponentsSet

  !
  !================================================================================================================================
  !

  !>Changes/sets the number of elements for a mesh. \see OpenCMISS::Iron::cmfe_Mesh_NumberOfElementsSet
  SUBROUTINE Mesh_NumberOfElementsSet(mesh,numberOfElements,err,error,*)

    !Argument variables
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh to set the number of elements for
    INTEGER(INTG), INTENT(IN) :: numberOfElements !<The number of elements to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx
    TYPE(MeshElementsType), POINTER :: meshElements
    TYPE(MeshTopologyType), POINTER :: meshTopology
    TYPE(VARYING_STRING) :: localError

    ENTERS("Mesh_NumberOfElementsSet",err,error,*999)

    CALL Mesh_AssertNotFinished(mesh,err,error,*999)
    
    IF(numberOfElements<=0) THEN
      localError="The specified number of elements of "//TRIM(NumberToVString(numberOfElements,"*",err,error))// &
        & " is invalid. You must have > 0 elements."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(numberOfElements/=mesh%numberOfElements) THEN
      DO componentIdx=1,mesh%numberOfComponents
        NULLIFY(meshTopology)
        CALL Mesh_MeshTopologyGet(mesh,componentIdx,meshTopology,err,error,*999)
        NULLIFY(meshElements)
        CALL MeshTopology_MeshElementsGet(meshTopology,meshElements,err,error,*999)
        IF(meshElements%numberOfElements>0) THEN
!!TODO: Reallocate the elements and copy information. 
          CALL FlagError("Not implemented.",err,error,*999)
        ENDIF
      ENDDO !componentIdx
      mesh%numberOfElements=numberOfElements
    ENDIF
    
    EXITS("Mesh_NumberOfElementsSet")
    RETURN
999 ERRORSEXITS("Mesh_NumberOfElementsSet",err,error)    
    RETURN 1
   
  END SUBROUTINE Mesh_NumberOfElementsSet

  !
  !================================================================================================================================
  !

  !>Changes/sets the surrounding elements calculate flag. \see OpenCMISS::Iron::cmfe_Mesh_SurroundingElementsCalculateSet
  SUBROUTINE Mesh_SurroundingElementsCalculateSet(mesh,surroundingElementsCalculateFlag,err,error,*)

    !Argument variables
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh to set the surrounding elements calculate flag for
    LOGICAL, INTENT(IN) :: surroundingElementsCalculateFlag !<The surrounding elements calculate flag
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("Mesh_SurroundingElementsCalculateSet",err,error,*999)

    CALL Mesh_AssertNotFinished(mesh,err,error,*999)
    
    mesh%surroundingElementsCalculate=surroundingElementsCalculateFlag
   
    EXITS("Mesh_SurroundingElementsCalculateSet")
    RETURN
999 ERRORSEXITS("Mesh_SurroundingElementsCalculateSet",err,error)    
    RETURN 1
   
  END SUBROUTINE Mesh_SurroundingElementsCalculateSet

  !
  !================================================================================================================================
  !

  !>Calculates the mesh topology.
  SUBROUTINE MeshTopology_Calculate(meshTopology,err,error,*)    

    !Argument variables
    TYPE(MeshTopologyType), POINTER :: meshTopology !<A pointer to the mesh topology to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshTopology_Calculate",err,error,*999)

    IF(.NOT.ASSOCIATED(meshTopology)) CALL FlagError("Mesh topology is not associated.",err,error,*999)
    
    !Calculate the nodes used in the mesh
    CALL MeshTopology_NodesCalculate(meshTopology,err,error,*999)
    !Calculate the elements surrounding the nodes in a mesh
    CALL MeshTopology_SurroundingElementsCalculate(meshTopology,err,error,*999)
    !Calculate the number of derivatives at each node in a mesh
    CALL MeshTopology_NodesDerivativesCalculate(meshTopology,err,error,*999)
    !Calculate the number of versions for each derivative at each node in a mesh
    CALL MeshTopology_NodesVersionCalculate(meshTopology,err,error,*999)
    !Calculate the elements surrounding the elements in the mesh
    CALL MeshTopology_ElementsAdjacentElementsCalculate(meshTopology,err,error,*999)
    !Calculate the boundary nodes and elements in the mesh
    CALL MeshTopology_BoundaryCalculate(meshTopology,err,error,*999)
    !Calculate the elements surrounding the elements in the mesh
    CALL MeshTopology_DofsCalculate(meshTopology,err,error,*999)
    
    EXITS("MeshTopology_Calculate")
    RETURN
999 ERRORSEXITS("MeshTopology_Calculate",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_Calculate
  
  !
  !===============================================================================================================================
  !

  !>Calculates the boundary nodes and elements for a mesh topology. 
  SUBROUTINE MeshTopology_BoundaryCalculate(meshTopology,err,error,*)

    !Argument variables
    TYPE(MeshTopologyType), POINTER :: meshTopology !<A pointer to the mesh topology to calculate the boundary for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementIdx,localNodeIdx,matchIndex,nodeIdx,xiCoordIdx,xiDirection
    TYPE(BasisType), POINTER :: basis
    TYPE(MeshElementsType), POINTER :: meshElements
    TYPE(MeshNodesType), POINTER :: meshNodes
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("MeshTopology_BoundaryCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(meshTopology)) CALL FlagError("Mesh topology is not associated.",err,error,*999)
    
    NULLIFY(meshNodes)
    CALL MeshTopology_MeshNodesGet(meshTopology,meshNodes,err,error,*999)
    NULLIFY(meshElements)
    CALL MeshTopology_MeshElementsGet(meshTopology,meshElements,err,error,*999)
    DO elementIdx=1,meshElements%numberOfElements
      NULLIFY(basis)
      CALL MeshElements_BasisGet(meshElements,elementIdx,basis,err,error,*999)
      SELECT CASE(basis%type)
      CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
        DO xiCoordIdx=-basis%numberOfXiCoordinates,basis%numberOfXiCoordinates
          IF(xiCoordIdx/=0) THEN
            IF(meshElements%elements(elementIdx)%adjacentElements(xiCoordIdx)%numberOfAdjacentElements==0) THEN
              meshElements%elements(elementIdx)%boundaryElement=.TRUE.
              IF(xiCoordIdx<0) THEN
                xiDirection=-xiCoordIdx
                matchIndex=1
              ELSE
                xiDirection=xiCoordIdx
                matchIndex=basis%numberOfNodesXiC(xiCoordIdx)
              ENDIF
              DO localNodeIdx=1,basis%numberOfNodes
                IF(basis%nodePositionIndex(localNodeIdx,xiDirection)==matchIndex) THEN
                  nodeIdx=meshElements%elements(elementIdx)%meshElementNodes(localNodeIdx)
                  meshNodes%nodes(nodeIdx)%boundaryNode=.TRUE.
                ENDIF
              ENDDO !localNodeIdx
            ENDIF
          ENDIF
        ENDDO !xiCoordIdx            
      CASE(BASIS_SIMPLEX_TYPE)
        meshElements%elements(elementIdx)%boundaryElement=.FALSE.
        DO xiCoordIdx=1,basis%numberOfXiCoordinates
          meshElements%elements(elementIdx)%boundaryElement=meshElements%elements(elementIdx)%boundaryElement.OR. &
            & meshElements%elements(elementIdx)%adjacentElements(xiCoordIdx)%numberOfAdjacentElements==0
          IF(meshElements%elements(elementIdx)%adjacentElements(xiCoordIdx)%numberOfAdjacentElements==0) THEN
            DO localNodeIdx=1,basis%numberOfNodes
              IF(basis%nodePositionIndex(localNodeIdx,xiCoordIdx)==1) THEN
                nodeIdx=meshElements%elements(elementIdx)%meshElementNodes(localNodeIdx)
                meshNodes%nodes(nodeIdx)%boundaryNode=.TRUE.
              ENDIF
            ENDDO !localNodeIdx                  
          ENDIF
        ENDDO !xiCoordIdx
      CASE(BASIS_SERENDIPITY_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(BASIS_AUXILLIARY_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(BASIS_B_SPLINE_TP_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(BASIS_FOURIER_LAGRANGE_HERMITE_TP_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(BASIS_EXTENDED_LAGRANGE_TP_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The basis type of "//TRIM(NumberToVString(basis%type,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDDO !elementIdx

    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Boundary elements:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Number of elements = ",meshElements%numberOfElements,err,error,*999)
      DO elementIdx=1,meshElements%numberOfElements
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Element : ",elementIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Boundary element = ",meshElements%elements(elementIdx)%boundaryElement, &
          & err,error,*999)        
      ENDDO !elementIdx
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Boundary nodes:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Number of nodes = ",meshNodes%numberOfNodes,err,error,*999)
      DO nodeIdx=1,meshNodes%numberOfNodes
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Node : ",nodeIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Boundary node = ",meshNodes%nodes(nodeIdx)%boundaryNode,err,error,*999)        
      ENDDO !nodeIdx            
    ENDIF
 
    EXITS("MeshTopology_BoundaryCalculate")
    RETURN
999 ERRORSEXITS("MeshTopology_BoundaryCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_BoundaryCalculate

  !
  !===============================================================================================================================
  !

  !>Calculates the degrees-of-freedom for a mesh topology. 
  SUBROUTINE MeshTopology_DofsCalculate(meshTopology,err,error,*)

    !Argument variables
    TYPE(MeshTopologyType), POINTER :: meshTopology !<A pointer to the mesh topology to calculate the dofs for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: derivativeIdx,nodeIdx,numberOfDofs,versionIdx
    TYPE(MeshDofsType), POINTER :: meshDofs
    TYPE(MeshNodesType), POINTER :: meshNodes

    ENTERS("MeshTopology_DofsCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(meshTopology)) CALL FlagError("Mesh topology is not associated.",err,error,*999)
    
    NULLIFY(meshNodes)
    CALL MeshTopology_MeshNodesGet(meshTopology,meshNodes,err,error,*999)
    NULLIFY(meshDofs)
    CALL MeshTopology_MeshDofsGet(meshTopology,meshDofs,err,error,*999)
    numberOfDofs=0
    DO nodeIdx=1,meshNodes%numberOfNodes
      DO derivativeIdx=1,meshNodes%nodes(nodeIdx)%numberOfDerivatives
        ALLOCATE(meshNodes%nodes(nodeIdx)%derivatives(derivativeIdx)%dofIndex( &
          & meshNodes%nodes(nodeIdx)%derivatives(derivativeIdx)%numberOfVersions),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate mesh topology node derivative version dof index.",err,error,*999)
        DO versionIdx=1,meshNodes%nodes(nodeIdx)%derivatives(derivativeIdx)%numberOfVersions
          numberOfDofs=numberOfDofs+1
          meshNodes%nodes(nodeIdx)%derivatives(derivativeIdx)%dofIndex(versionIdx)=numberOfDofs
        ENDDO !versionIdx
      ENDDO !derivativeIdx
    ENDDO !nodeIdx
    meshDofs%numberOfDofs=numberOfDofs
  
    EXITS("MeshTopology_DofsCalculate")
    RETURN
999 ERRORSEXITS("MeshTopology_DofsCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_DofsCalculate

  !
  !===============================================================================================================================
  !

  !>Finalises the dof data structures for a mesh topology and deallocates any memory. \todo pass in dofs
  SUBROUTINE MeshTopology_DofsFinalise(meshDofs,err,error,*)

    !Argument variables
    TYPE(MeshDofsType), POINTER :: meshDofs !<A pointer to the mesh topology to finalise the dofs for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshTopology_DofsFinalise",err,error,*999)

    IF(ASSOCIATED(meshDofs)) THEN
      DEALLOCATE(meshDofs)
    ENDIF
 
    EXITS("MeshTopology_DofsFinalise")
    RETURN
999 ERRORSEXITS("MeshTopology_DofsFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_DofsFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the dofs in a given mesh topology.
  SUBROUTINE MeshTopology_DofsInitialise(meshTopology,err,error,*)

    !Argument variables
    TYPE(MeshTopologyType), POINTER :: meshTopology !<A pointer to the mesh topology to initialise the dofs for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("MeshTopology_DofsInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(meshTopology)) CALL FlagError("Mesh topology is not associated.",err,error,*998)
    IF(ASSOCIATED(meshTopology%dofs)) CALL FlagError("Mesh already has topology dofs associated.",err,error,*998)
    
    ALLOCATE(meshTopology%dofs,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate mesh topology dofs.",err,error,*999)
    meshTopology%dofs%numberOfDofs=0
    meshTopology%dofs%meshTopology=>meshTopology
     
    EXITS("MeshTopology_DofsInitialise")
    RETURN
999 CALL MeshTopology_DofsFinalise(meshTopology%dofs,dummyErr,dummyError,*998)
998 ERRORSEXITS("MeshTopology_DofsInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_DofsInitialise

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating elements for a specified mesh component in a mesh topology. \see OpenCMISS::Iron::cmfe_MeshElements_CreateFinish
  SUBROUTINE MeshElements_CreateFinish(meshElements,err,error,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: meshElements !<A pointer to the mesh elements to finish creating
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementIdx
    TYPE(BasisType), POINTER :: basis
    TYPE(MeshType), POINTER :: mesh
    TYPE(MeshTopologyType), POINTER :: meshTopology
    TYPE(VARYING_STRING) :: localError

    ENTERS("MeshElements_CreateFinish",err,error,*999)

    CALL MeshElements_AssertNotFinished(meshElements,err,error,*999)

    meshElements%elementsFinished=.TRUE.
    
    IF(diagnostics1) THEN
      NULLIFY(meshTopology)
      CALL MeshElements_MeshTopologyGet(meshElements,meshTopology,err,error,*999)
      NULLIFY(mesh)
      CALL MeshTopology_MeshGet(meshTopology,mesh,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Number of global elements = ",mesh%numberOfElements,err,error,*999)
      DO elementIdx=1,mesh%numberOfElements
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Element = ",elementIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Global number        = ", &
          & meshElements%elements(elementIdx)%globalNumber,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    User number          = ", &
          & meshElements%elements(elementIdx)%userNumber,err,error,*999)
        NULLIFY(basis)
        CALL MeshElements_BasisGet(meshElements,elementIdx,basis,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Basis number         = ",basis%userNumber,err,error,*999)
        IF(.NOT.ALLOCATED(meshElements%elements(elementIdx)%userElementNodes)) THEN
          localError="User element nodes are not associated for element number "// &
            & TRIM(NumberToVString(elementIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,basis%numberOfNodes,8,8, &
          & meshElements%elements(elementIdx)%userElementNodes,'("    User element nodes   =",8(X,I6))','(26X,8(X,I6))', &
          & err,error,*999)
        IF(.NOT.ALLOCATED(meshElements%elements(elementIdx)%globalElementNodes)) THEN
          localError="Global element nodes are not associated for element number "// &
            & TRIM(NumberToVString(elementIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,basis%numberOfNodes,8,8, &
          & meshElements%elements(elementIdx)%globalElementNodes,'("    Global element nodes =",8(X,I6))','(26X,8(X,I6))', &
          & err,error,*999)
      ENDDO !elementIdx
    ENDIF
    
    EXITS("MeshElements_CreateFinish")
    RETURN
999 ERRORSEXITS("MeshElements_CreateFinish",err,error)
    RETURN 1
   
  END SUBROUTINE MeshElements_CreateFinish
    
  !
  !================================================================================================================================
  !

  !>Starts the process of creating elements in the mesh component identified by MESH and componentIdx. The elements will be created with a default basis of BASIS. ELEMENTS is the returned pointer to the MESH_ELEMENTS data structure. \see OpenCMISS::Iron::cmfe_MeshElementsCreateStart
  SUBROUTINE MeshElements_CreateStart(mesh,meshComponentNumber,basis,meshElements,err,error,*)

    !Argument variables
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh to start creating the elements on
    INTEGER(INTG), INTENT(IN) :: meshComponentNumber !<The mesh component number
    TYPE(BasisType), POINTER :: basis !<A pointer to the default basis to use
    TYPE(MeshElementsType), POINTER :: meshElements !<On return, a pointer to the created mesh elements
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,elementIdx,insertStatus
    TYPE(MeshTopologyType), POINTER :: meshTopology
    TYPE(VARYING_STRING) :: dummyError,localError
 
    ENTERS("MeshElements_CreateStart",err,error,*999)
    
    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Mesh is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)
    IF(ASSOCIATED(meshElements)) CALL FlagError("Mesh elements is already associated.",err,error,*999)
    IF(meshComponentNumber<=0.OR.meshComponentNumber>mesh%numberOfComponents) THEN
      localError="The specified mesh component number of "//TRIM(NumberToVString(meshComponentNumber,"*",err,error))// &
        & " is invalid. The component number must be between 1 and "// &
        & TRIM(NumberToVString(mesh%numberOfComponents,"*",err,error))
      CALL FlagError(localError,err,error,*998)
    ENDIF

    NULLIFY(meshTopology)
    CALL Mesh_MeshTopologyGet(mesh,meshComponentNumber,meshTopology,err,error,*999)
    NULLIFY(meshElements)
    CALL MeshTopology_MeshElementsGet(meshTopology,meshElements,err,error,*999)
    IF(ALLOCATED(meshElements%elements)) CALL FlagError("Mesh topology already has elements associated.",err,error,*999)
    meshTopology%meshComponentNumber=meshComponentNumber
    ALLOCATE(meshElements%elements(mesh%numberOfElements),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate individual elements.",err,error,*999)
    meshElements%numberOfElements=mesh%numberOfElements !Psuedo inheritance of the number of elements
    CALL Tree_CreateStart(meshElements%elementsTree,err,error,*999)
    CALL Tree_InsertTypeSet(meshElements%elementsTree,TREE_NO_DUPLICATES_ALLOWED,err,error,*999)
    CALL Tree_CreateFinish(meshElements%elementsTree,err,error,*999)
    meshElements%elementsFinished=.FALSE.
    !Set up the default values and allocate element structures
    DO elementIdx=1,meshElements%numberOfElements
      CALL MeshElement_Initialise(meshElements%elements(elementIdx),err,error,*999)
      meshElements%elements(elementIdx)%globalNumber=elementIdx
      meshElements%elements(elementIdx)%userNumber=elementIdx
      CALL Tree_ItemInsert(meshElements%elementsTree,elementIdx,elementIdx,insertStatus,err,error,*999)
      meshElements%elements(elementIdx)%basis=>basis
      ALLOCATE(meshElements%elements(elementIdx)%userElementNodes(basis%numberOfNodes),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate user element nodes.",err,error,*999)
      ALLOCATE(meshElements%elements(elementIdx)%globalElementNodes(basis%numberOfNodes),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate global element nodes.",err,error,*999)
      meshElements%elements(elementIdx)%userElementNodes=1
      meshElements%elements(elementIdx)%globalElementNodes=1
      ALLOCATE(meshElements%elements(elementIdx)%userElementNodeVersions(basis%maximumNumberOfDerivatives, &
        & basis%numberOfNodes),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate global element nodes versions.",err,error,*999)
      meshElements%elements(elementIdx)%userElementNodeVersions = 1
    ENDDO !elementIdx

    EXITS("MeshElements_CreateStart")
    RETURN
999 CALL MeshTopology_ElementsFinalise(meshElements,dummyErr,dummyError,*998)
998 ERRORSEXITS("MeshElements_CreateStart",err,error)
    RETURN 1
   
  END SUBROUTINE MeshElements_CreateStart

  !
  !================================================================================================================================
  !

  !>Destroys the elements in a mesh topology. \todo as this is a user routine it should take a mesh pointer like create start and finish? Split this into destroy and finalise?
  SUBROUTINE MeshElements_Destroy(meshElements,err,error,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: meshElements !<A pointer to the mesh elements to destroy 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshElements_Destroy",err,error,*999)

    IF(.NOT.ASSOCIATED(meshElements)) CALL FlagError("Mesh elements is not associated.",err,error,*999)
    
    CALL MeshTopology_ElementsFinalise(meshElements,err,error,*999)

    EXITS("MeshElements_Destroy")
    RETURN
999 ERRORSEXITS("MeshElements_Destroy",err,error)
    RETURN 1
    
  END SUBROUTINE MeshElements_Destroy
  
  !
  !================================================================================================================================
  !

  !>Finalises the given mesh topology element.
  SUBROUTINE MeshElement_Finalise(meshElement,err,error,*)
    
    !Argument variables
    TYPE(MeshElementType) :: meshElement !<The mesh element to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: adjacentElementIdx
    
    ENTERS("MeshElement_Finalise",err,error,*999)
    
    IF(ALLOCATED(meshElement%userElementNodeVersions)) DEALLOCATE(meshElement%userElementNodeVersions)
    IF(ALLOCATED(meshElement%userElementNodes)) DEALLOCATE(meshElement%userElementNodes)
    IF(ALLOCATED(meshElement%globalElementNodes)) DEALLOCATE(meshElement%globalElementNodes)
    IF(ALLOCATED(meshElement%meshElementNodes)) DEALLOCATE(meshElement%meshElementNodes)
    IF(ALLOCATED(meshElement%adjacentElements)) THEN
      DO adjacentElementIdx=LBOUND(meshElement%adjacentElements,1),UBOUND(meshElement%adjacentElements,1)
        CALL Mesh_AdjacentElementFinalise(meshElement%adjacentElements(adjacentElementIdx),err,error,*999)
      ENDDO !adjacentElementIdx
      DEALLOCATE(meshElement%adjacentElements)
    ENDIF
    
    EXITS("MeshElement_Finalise")
    RETURN
999 ERRORSEXITS("MeshElement_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE MeshElement_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the given mesh element.
  SUBROUTINE MeshElement_Initialise(meshElement,err,error,*)

    !Argument variables
    TYPE(MeshElementType) :: meshElement !<The mesh element to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshElement_Initialise",err,error,*999)

    meshElement%userNumber=0
    meshElement%globalNumber=0
    NULLIFY(meshElement%basis)
    meshElement%boundaryElement=.FALSE.
    
    EXITS("MeshElement_Initialise")
    RETURN
999 ERRORSEXITS("MeshElement_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE MeshElement_Initialise
  
  !
  !================================================================================================================================
  !

  !>Gets the basis for a mesh element identified by a given global number. \see OpenCMISS::Iron::cmfe_MeshElements_BasisGet
  SUBROUTINE MeshElements_ElementBasisGet(meshElements,userElementNumber,basis,err,error,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: meshElements !<A pointer to the mesh elements to get the basis for
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user number of the element to get the basis for
    TYPE(BasisType), POINTER :: basis !<On return, a pointer to the basis to get. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: globalElementNumber

    ENTERS("MeshElements_ElementBasisGet",err,error,*999)

    CALL MeshElements_AssertIsFinished(meshElements,err,error,*999)
    IF(ASSOCIATED(basis)) CALL FlagError("Basis is already associated.",err,error,*999)    
    CALL MeshElements_GlobalElementNumberGet(meshElements,userElementNumber,globalElementNumber,err,error,*999)
    
    CALL MeshElements_BasisGet(meshElements,globalElementNumber,basis,err,error,*999)
    
    EXITS("MeshElements_ElementBasisGet")
    RETURN
999 ERRORSEXITS("MeshElements_ElementBasisGet",err,error)    
    RETURN 1
    
  END SUBROUTINE MeshElements_ElementBasisGet

  !
  !================================================================================================================================
  !

  !>Changes/sets the basis for a mesh element identified by a given user number. 
  SUBROUTINE MeshElements_ElementBasisSet(meshElements,userElementNumber,basis,err,error,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: meshElements !<A pointer to the elements to set the basis for
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user number of the element to set the basis for
    TYPE(BasisType), POINTER :: basis !<A pointer to the basis to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: globalElementNumber,overlappingNumberNodes,overlappingNumberDerivatives
    INTEGER(INTG), ALLOCATABLE :: newUserElementNodes(:),newGlobalElementNodes(:),newUserElementNodeVersions(:,:)
    TYPE(BasisType), POINTER :: meshBasis
 
    ENTERS("MeshElements_ElementBasisSet",err,error,*999)

    CALL MeshElements_AssertNotFinished(meshElements,err,error,*999)
    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)
    CALL MeshElements_GlobalElementNumberGet(meshElements,userElementNumber,globalElementNumber,err,error,*999)
    
    NULLIFY(meshBasis)
    CALL MeshElements_BasisGet(meshElements,globalElementNumber,meshBasis,err,error,*999)
    IF(meshBasis%numberOfNodes/=basis%numberOfNodes.OR. &
      & meshBasis%maximumNumberOfDerivatives/=basis%maximumNumberOfDerivatives) THEN      
      !Allocate new user and global element nodes
      ALLOCATE(newUserElementNodes(basis%numberOfNodes),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new user element nodes",err,error,*999)
      ALLOCATE(newGlobalElementNodes(basis%numberOfNodes),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new user element nodes",err,error,*999)
      ALLOCATE(newUserElementNodeVersions(basis%maximumNumberOfDerivatives,basis%numberOfNodes),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate element node versions",err,error,*999)      
      overlappingNumberNodes=MIN(basis%numberOfNodes,meshBasis%numberOfNodes)
      overlappingNumberDerivatives=MIN(basis%maximumNumberOfDerivatives,meshBasis%maximumNumberOfDerivatives)
      !Set default values
      newUserElementNodeVersions=1
      newUserElementNodes(overlappingNumberNodes+1:)=0
      newGlobalElementNodes(overlappingNumberNodes+1:)=0
      !Copy previous values
      newUserElementNodes(1:overlappingNumberNodes)=meshElements%elements(globalElementNumber)% &
        & userElementNodes(1:overlappingNumberNodes)
      newGlobalElementNodes(1:overlappingNumberNodes)=meshElements%elements(globalElementNumber)% &
        & globalElementNodes(1:overlappingNumberNodes)
      newUserElementNodeVersions(1:overlappingNumberDerivatives,1:overlappingNumberNodes)= &
        & meshElements%elements(globalElementNumber)%userElementNodeVersions(1:overlappingNumberDerivatives, &
        & 1:overlappingNumberNodes)
      !Replace arrays with new ones
      CALL MOVE_ALLOC(newUserElementNodeVersions,meshElements%elements(globalElementNumber)%userElementNodeVersions)
      CALL MOVE_ALLOC(newUserElementNodes,meshElements%elements(globalElementNumber)%userElementNodes)
      CALL MOVE_ALLOC(newGlobalElementNodes,meshElements%elements(globalElementNumber)%globalElementNodes)
    ENDIF
    meshElements%elements(globalElementNumber)%basis=>basis
    
    EXITS("MeshElements_ElementBasisSet")
    RETURN
999 IF(ALLOCATED(newUserElementNodeVersions)) DEALLOCATE(newUserElementNodeVersions)
    IF(ALLOCATED(newUserElementNodes)) DEALLOCATE(newUserElementNodes)  
    IF(ALLOCATED(newGlobalElementNodes)) DEALLOCATE(newGlobalElementNodes)  
    ERRORSEXITS("MeshElements_ElementBasisSet",err,error)    
    RETURN 1
    
  END SUBROUTINE MeshElements_ElementBasisSet
  
  !
  !================================================================================================================================
  !

  !>Returns the adjacent element number for a mesh element identified by a user number. \see OpenCMISS::Iron::cmfe_MeshElements_AdjacentElementGet
  SUBROUTINE MeshElements_AdjacentElementGet(meshElements,userElementNumber,adjacentElementXi,adjacentUserNumber,err,error,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: meshElements !<A pointer to the elements of a mesh component from which to get the adjacent element from.
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user number of the element to get the adjacent element for
    INTEGER(INTG), INTENT(IN) :: adjacentElementXi !< The xi coordinate direction to get the adjacent element Note that -xiCoordinateDirection gives the adjacent element before the element in the xiCoordinateDirection'th direction and +xiCoordinateDirection gives the adjacent element after the element in the xiCoordinateDirection'th direction. The xiCoordinateDirection=0 index will give the information on the current element.
    INTEGER(INTG), INTENT(OUT) :: adjacentUserNumber !<On return, the adjacent element number in the specified xi coordinate direction. Return 0 if the specified element has no adjacent elements in the specified xi coordinate direction.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: globalElementNumber
    TYPE(BasisType), POINTER :: basis
    TYPE(VARYING_STRING) :: localError

    ENTERS("MeshElements_AdjacentElementGet",err,error,*999)

    CALL MeshElements_AssertIsFinished(meshElements,err,error,*999)
    CALL MeshElements_GlobalElementNumberGet(meshElements,userElementNumber,globalElementNumber,err,error,*999)
    NULLIFY(basis)
    CALL MeshElements_BasisGet(meshElements,globalElementNumber,basis,err,error,*999)      
    IF(adjacentElementXi<-basis%numberOfXi.OR.adjacentElementXi>basis%numberOfXi) THEN
      localError="The specified adjacent element xi is invalid. The supplied xi is "// &
        & TRIM(NumberToVString(adjacentElementXi,"*",err,error))//" and needs to be >= -"// &
        & TRIM(NumberToVString(basis%numberOfXi,"*",err,error))//" and <= "// &
        & TRIM(NumberToVString(basis%numberOfXi,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(meshElements%elements(globalElementNumber)%adjacentElements(adjacentElementXi)%numberOfAdjacentElements > 0) THEN !\todo Currently returns only the first adjacent element for now as the python binding require the output array size of the adjacent element to be known a-prior. Add routine to first output number of adjacent elements and then loop over all adjacent elements
      adjacentUserNumber=meshElements%elements(globalElementNumber)%adjacentElements(adjacentElementXi)%adjacentElements(1)
    ELSE !Return 0 indicating the specified element has no adjacent elements in the specified xi coordinate direction.
      adjacentUserNumber=0
    ENDIF

    EXITS("MeshElements_AdjacentElementGet")
    RETURN
999 ERRORSEXITS("MeshElements_AdjacentElementGet",err,error)
    RETURN 1

  END SUBROUTINE MeshElements_AdjacentElementGet

  !
  !================================================================================================================================
  !

  !>Gets the element nodes for a mesh element identified by a given global number. \see OpenCMISS::Iron::cmfe_MeshElements_NodesGet
  SUBROUTINE MeshElements_ElementNodesGet(meshElements,userElementNumber,userElementNodes,err,error,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: meshElements !<A pointer to the elements to get
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user number of the element to get the nodes for
    INTEGER(INTG), INTENT(OUT) :: userElementNodes(:) !<On return, userElementNodes(i). userElementNodes(i) is the i'th user node number for the element
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: globalElementNumber
     TYPE(VARYING_STRING) :: localError

    ENTERS("MeshElements_ElementNodesGet",err,error,*999)

    CALL MeshElements_AssertIsFinished(meshElements,err,error,*999)
    CALL MeshElements_GlobalElementNumberGet(meshElements,userElementNumber,globalElementNumber,err,error,*999)
    IF(SIZE(userElementNodes,1)<SIZE(meshElements%elements(globalElementNumber)%userElementNodes,1)) THEN
      localError="The size of user element nodes is too small. The supplied size is "// &
        & TRIM(NumberToVString(SIZE(userElementNodes,1),"*",err,error))//" and it needs to be >= "// &
        & TRIM(NumberToVString(SIZE(meshElements%elements(globalElementNumber)%userElementNodes,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    userElementNodes(1:SIZE(meshElements%elements(globalElementNumber)%userElementNodes,1))= &
      & meshElements%elements(globalElementNumber)%userElementNodes
    
    EXITS("MeshElements_ElementNodesGet")
    RETURN
999 ERRORSEXITS("MeshElements_ElementNodesGet",err,error)    
    RETURN 1
    
  END SUBROUTINE MeshElements_ElementNodesGet

  !
  !================================================================================================================================
  !

  !>Changes/sets the element nodes for a mesh element identified by a given global number. \see OpenCMISS::Iron::cmfe_MeshElements_NodesSet
  SUBROUTINE MeshElements_ElementNodesSet(meshElements,userElementNumber,userElementNodes,err,error,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: meshElements !<A pointer to the elements to set 
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user number of the element to set the nodes for
    INTEGER(INTG), INTENT(IN) :: userElementNodes(:) !<userElementNodes(i). userElementNodes(i) is the i'th user node number for the element
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: globalElementNumber,globalNodeNumber,localNodeIdx,numberOfBadNodes
    INTEGER(INTG), ALLOCATABLE :: globalElementNodes(:),badNodes(:)
    LOGICAL :: elementNodesOK,nodeExists
    TYPE(BasisType), POINTER :: basis
    TYPE(InterfaceType), POINTER :: interface
    TYPE(MeshType), POINTER :: mesh
    TYPE(MeshTopologyType), POINTER :: meshTopology
    TYPE(NodesType), POINTER :: nodes
    TYPE(RegionType), POINTER :: parentRegion,region
    TYPE(VARYING_STRING) :: localError

    ENTERS("MeshElements_ElementNodesSet",err,error,*999)

    CALL MeshElements_AssertNotFinished(meshElements,err,error,*999)
    CALL MeshElements_GlobalElementNumberGet(meshElements,userElementNumber,globalElementNumber,err,error,*999)
    NULLIFY(basis)
    CALL MeshElements_BasisGet(meshElements,globalElementNumber,basis,err,error,*999)      
    IF(SIZE(userElementNodes,1)<basis%numberOfNodes) THEN
      localError="The size of user element nodes is too small. The supplied size is "// &
        & TRIM(NumberToVString(SIZE(userElementNodes,1),"*",err,error))//" and it needs to be >= "// &
        & TRIM(NumberToVString(SIZE(meshElements%elements(globalElementNumber)%userElementNodes,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    NULLIFY(meshTopology)
    CALL MeshElements_MeshTopologyGet(meshElements,meshTopology,err,error,*999)
    NULLIFY(mesh)
    CALL MeshTopology_MeshGet(meshTopology,mesh,err,error,*999)
    region=>mesh%region
    IF(ASSOCIATED(region)) THEN
      nodes=>region%nodes
    ELSE
      interface=>mesh%interface
      IF(ASSOCIATED(interface)) THEN
        nodes=>interface%nodes
        parentRegion=>interface%parentRegion
        IF(.NOT.ASSOCIATED(parentRegion)) CALL FlagError("Mesh interface has no parent region.",err,error,*999)
      ELSE
        CALL FlagError("Elements mesh has no associated region or interface.",err,error,*999)
      ENDIF
    ENDIF
    IF(.NOT.ASSOCIATED(nodes)) CALL FlagError("Mesh has no associated nodes.",err,error,*999)
    elementNodesOK=.TRUE.    
    ALLOCATE(globalElementNodes(basis%numberOfNodes),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate global element nodes.",err,error,*999)
    ALLOCATE(badNodes(basis%numberOfNodes),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate bad nodes.",err,error,*999)
    numberOfBadNodes=0
    DO localNodeIdx=1,basis%numberOfNodes
      CALL Nodes_NodeCheckExists(nodes,userElementNodes(localNodeIdx),nodeExists,globalNodeNumber,err,error,*999)
      IF(nodeExists) THEN
        globalElementNodes(localNodeIdx)=globalNodeNumber
      ELSE
        numberOfBadNodes=numberOfBadNodes+1
        badNodes(numberOfBadNodes)=userElementNodes(localNodeIdx)
        elementNodesOK=.FALSE.
      ENDIF
    ENDDO !localNodeIdx
    IF(elementNodesOK) THEN
      meshElements%elements(globalElementNumber)%userElementNodes(1:basis%numberOfNodes)= &
        & userElementNodes(1:basis%numberOfNodes)
      meshElements%elements(globalElementNumber)%globalElementNodes(1:basis%numberOfNodes)= &
        & globalElementNodes(1:basis%numberOfNodes)      
    ELSE
      IF(numberOfBadNodes==1) THEN
        IF(ASSOCIATED(region)) THEN
          localError="The element user node number of "//TRIM(NumberToVString(badNodes(1),"*",err,error))// &
            & " is not defined in region "//TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
        ELSE
          localError="The element user node number of "//TRIM(NumberToVString(badNodes(1),"*",err,error))// &
            & " is not defined in interface number "// &
            & TRIM(NumberToVString(interface%userNumber,"*",err,error))// &
            & " of parent region number "//TRIM(NumberToVString(parentRegion%userNumber,"*",err,error))//"."
        ENDIF
      ELSE
        localError="The element user node number of "//TRIM(NumberToVString(badNodes(1),"*",err,error))
        DO localNodeIdx=2,numberOfBadNodes-1
          localError=localError//","//TRIM(NumberToVString(badNodes(localNodeIdx),"*",err,error))
        ENDDO !localNodeIdx
        IF(ASSOCIATED(region)) THEN
          localError=localError//" & "//TRIM(NumberToVString(badNodes(numberOfBadNodes),"*",err,error))// &
            & " are not defined in region number "//TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
        ELSE
          localError=localError//" & "//TRIM(NumberToVString(badNodes(numberOfBadNodes),"*",err,error))// &
            & " are not defined in interface number "// &
            & TRIM(NumberToVString(interface%userNumber,"*",err,error))//" of parent region number "// &
            & TRIM(NumberToVString(parentRegion%userNumber,"*",err,error))//"."
        ENDIF
      ENDIF
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("MeshElements_ElementNodesSet")
    RETURN
999 ERRORSEXITS("MeshElements_ElementNodesSet",err,error)    
    RETURN 1
    
  END SUBROUTINE MeshElements_ElementNodesSet

  !
  !================================================================================================================================
  !

  !>Changes/sets an element node's version for a mesh element identified by a given user number. \see OpenCMISS::Iron::cmfe_MeshElements_LocalElementNodeVersionSet
  SUBROUTINE MeshElements_ElementNodeVersionSet(meshElements,userElementNumber,versionNumber,derivativeNumber, &
    & userElementNodeIndex,err,error,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: meshElements !<A pointer to the elements to set
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user number of the element to set the nodes for
    INTEGER(INTG), INTENT(IN) :: versionNumber !<The version number of the specified element node to set.
    INTEGER(INTG), INTENT(IN) :: derivativeNumber !<The derivative number of the specified element node to set.
    INTEGER(INTG), INTENT(IN) :: userElementNodeIndex !< The node index of the specified element node to set a version for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: globalElementNumber
    TYPE(BasisType), POINTER :: basis
    TYPE(VARYING_STRING) :: localError

    ENTERS("MeshElements_ElementNodeVersionSet",err,error,*999)

    CALL MeshElements_AssertNotFinished(meshElements,err,error,*999)
    CALL MeshElements_GlobalElementNumberGet(meshElements,userElementNumber,globalElementNumber,err,error,*999)
    NULLIFY(basis)
    CALL MeshElements_BasisGet(meshElements,globalElementNumber,basis,err,error,*999)            
    IF(userElementNodeIndex<1.OR.userElementNodeIndex>basis%numberOfNodes) THEN
      localError="The specified element local node index of "//TRIM(NumberToVString(userElementNodeIndex,"*",err,error))// &
        & " is invalid. The element local node index should be between 1 and "// &
        & TRIM(NumberToVString(basis%numberOfNodes,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(derivativeNumber<1.OR.derivativeNumber>basis%numberOfDerivatives(userElementNodeIndex)) THEN 
      localError="The specified node derivative number of "//TRIM(NumberToVString(derivativeNumber,"*",err,error))// &
        & " is invalid. The element node derivative index should be between 1 and "// &
        & TRIM(NumberToVString(basis%numberOfDerivatives(userElementNodeIndex),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(versionNumber<1) THEN 
      localError="The specified node version number of "//TRIM(NumberToVString(versionNumber,"*",err,error))// &
        & " is invalid. The element node index should be >= 1."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    meshElements%elements(globalElementNumber)%userElementNodeVersions(derivativeNumber,userElementNodeIndex)=versionNumber
    !\todo : There is redunancy in userElementNodeVersions since it was allocated in MeshElement_CreateStart based on maximumNumberOfDerivatives for that elements basis:ALLOCATE(meshElements%elements(ne)%userElementNodeVersions(BASIS%maximumNumberOfDerivatives,BASIS%numberOfNodes),STAT=err)
    
    EXITS("MeshElements_ElementNodeVersionSet")
    RETURN
999 ERRORSEXITS("MeshElements_ElementNodeVersionSet",err,error)    
    RETURN 1
    
  END SUBROUTINE MeshElements_ElementNodeVersionSet

  !
  !================================================================================================================================
  !

  !>Changes/sets an element node's version for a mesh element identified by a given user number. \see OpenCMISS::Iron::cmfe_MeshElements_UserNodeVersionSet
  SUBROUTINE MeshElements_UserNodeVersionSet(meshElements,userElementNumber,versionNumber,derivativeNumber, &
    & userNodeNumber,err,error,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: meshElements !<A pointer to the elements to set
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user number of the element to set the nodes for
    INTEGER(INTG), INTENT(IN) :: versionNumber !<The version number of the specified element node to set.
    INTEGER(INTG), INTENT(IN) :: derivativeNumber !<The derivative number of the specified element node to set.
    INTEGER(INTG), INTENT(IN) :: userNodeNumber !< The user node number of the specified element node to set a version for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: globalElementNumber,localNodeIdx
    LOGICAL :: found
    TYPE(BasisType), POINTER :: basis
    TYPE(MeshNodesType), POINTER :: meshNodes
    TYPE(MeshTopologyType), POINTER :: meshTopology
    TYPE(VARYING_STRING) :: localError

    ENTERS("MeshElements_ElementNodeVersionSet",err,error,*999)

    CALL MeshElements_AssertNotFinished(meshElements,err,error,*999)
    CALL MeshElements_GlobalElementNumberGet(meshElements,userElementNumber,globalElementNumber,err,error,*999)
    NULLIFY(basis)    
    CALL MeshElements_BasisGet(meshElements,globalElementNumber,basis,err,error,*999)
    NULLIFY(meshTopology)
    CALL MeshElements_MeshTopologyGet(meshElements,meshTopology,err,error,*999)
    NULLIFY(meshNodes)
    CALL MeshTopology_MeshNodesGet(meshTopology,meshNodes,err,error,*999)
    found=.FALSE.
    DO localNodeIdx=1,basis%numberOfNodes
      IF(meshElements%elements(globalElementNumber)%userElementNodes(localNodeIdx)==userNodeNumber) THEN
        found=.TRUE.
        EXIT
      ENDIF
    ENDDO !localNodeIdx
    IF(found) THEN
      CALL MeshElements_ElementNodeVersionSet(meshElements,userElementNumber,versionNumber,derivativeNumber, &
        & localNodeIdx,err,error,*999)
    ELSE
      localError="User node number "//TRIM(NumberToVString(userNodeNumber,"*",err,error))// &
        & " does not exist in user element number "//TRIM(NumberToVString(userElementNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
        
    EXITS("MeshElements_UserNodeVersionSet")
    RETURN
999 ERRORSEXITS("MeshElements_UserNodeVersionSet",err,error)    
    RETURN 1
    
  END SUBROUTINE MeshElements_UserNodeVersionSet

  !
  !================================================================================================================================
  !

 !>Calculates the element numbers surrounding an element in a mesh topology.
  SUBROUTINE MeshTopology_ElementsAdjacentElementsCalculate(meshTopology,err,error,*)

    !Argument variables
    TYPE(MeshTopologyType), POINTER :: meshTopology !<A pointer to the mesh topology to calculate the elements adjacent to elements for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: j,elementIdx,elementIdx1,surroundingElementIdx,xiIdx,xiCoordIdx,localNodeIdx,localNodeIdx1,localNodeIdx2, &
      & localNodeIdx3,nodeIdx,nodeNumber,nodeNumber1,dummyErr,faceXi(2),faceXiC(3),nodePositionIndex(4)
    INTEGER(INTG) :: xiDirection,directionIndex,xiDirCheck,xiDirSearch,numberNodeMatches
    INTEGER(INTG) :: numberSurrounding,numberOfNodesXiC(4)
    INTEGER(INTG), ALLOCATABLE :: nodeMatches(:),adjacentElements(:)
    LOGICAL :: xiCollapsed,faceCollapsed(-3:3),subset
    TYPE(ListType), POINTER :: nodeMatchList
    TYPE(ListPtrType) :: adjacentElementsList(-4:4)
    TYPE(BasisType), POINTER :: basis
    TYPE(MeshElementsType), POINTER :: meshElements
    TYPE(MeshNodesType), POINTER :: meshNOdes
    TYPE(VARYING_STRING) :: dummyError,localError

    NULLIFY(nodeMatchList)
    DO xiCoordIdx=-4,4
      NULLIFY(adjacentElementsList(xiCoordIdx)%ptr)
    ENDDO !xiCoordIdx
    
    ENTERS("MeshTopology_ElementsAdjacentElementsCalculate",err,error,*999)
    
    IF(.NOT.ASSOCIATED(meshTopology)) CALL FlagError("Mesh topology is not associated.",err,error,*999)
    
    NULLIFY(meshNodes)
    CALL MeshTopology_MeshNodesGet(meshTopology,meshNodes,err,error,*999)
    NULLIFY(meshElements)
    CALL MeshTopology_MeshElementsGet(meshTopology,meshElements,err,error,*999)
    
    !Loop over the global elements in the mesh
    DO elementIdx=1,meshElements%numberOfElements
      NULLIFY(basis)
      CALL MeshElements_BasisGet(meshElements,elementIdx,basis,err,error,*999)
      !First initialize lists that are required to find the adjacent elements list
      DO xiCoordIdx=-basis%numberOfXiCoordinates,basis%numberOfXiCoordinates
        NULLIFY(adjacentElementsList(xiCoordIdx)%ptr)
        CALL List_CreateStart(adjacentElementsList(xiCoordIdx)%ptr,err,error,*999)
        CALL List_DataTypeSet(adjacentElementsList(xiCoordIdx)%ptr,LIST_INTG_TYPE,err,error,*999)
        CALL List_InitialSizeSet(adjacentElementsList(xiCoordIdx)%ptr,5,err,error,*999)
        CALL List_CreateFinish(adjacentElementsList(xiCoordIdx)%ptr,err,error,*999)
      ENDDO !xiCoordIdx
      numberOfNodesXiC=1
      numberOfNodesXiC(1:basis%numberOfXiCoordinates)=basis%numberOfNodesXiC(1:basis%numberOfXiCoordinates)
      !Place the current element in the surrounding list
      CALL List_ItemAdd(adjacentElementsList(0)%ptr,meshElements%elements(elementIdx)%globalNumber,err,error,*999)
      SELECT CASE(basis%TYPE)
      CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
        !Determine the collapsed "faces" if any
        nodePositionIndex=1
        !Loop over the face normals of the element
        DO xiIdx=1,basis%numberOfXi
          !Determine the xi directions that lie in this xi direction
          faceXi(1)=OTHER_XI_DIRECTIONS3(xiIdx,2,1)
          faceXi(2)=OTHER_XI_DIRECTIONS3(xiIdx,3,1)
          !Reset the node_position_index in this xi direction
          nodePositionIndex(xiIdx)=1
          !Loop over the two faces with this normal
          DO directionIndex=-1,1,2
            xiDirection=directionIndex*xiIdx
            faceCollapsed(xiDirection)=.FALSE.
            DO j=1,2
              xiDirCheck=faceXi(j)
              IF(xiDirCheck<=basis%numberOfXi) THEN
                xiDirSearch=faceXi(3-j)
                nodePositionIndex(xiDirSearch)=1
                xiCollapsed=.TRUE.
                DO WHILE(nodePositionIndex(xiDirSearch)<=numberOfNodesXiC(xiDirSearch).AND.xiCollapsed)
                  !Get the first local node along the xi check direction
                  nodePositionIndex(xiDirCheck)=1
                  localNodeIdx1=basis%nodePositionIndexInv(nodePositionIndex(1),nodePositionIndex(2),nodePositionIndex(3),1)
                  !Get the second local node along the xi check direction
                  nodePositionIndex(xiDirCheck)=2
                  localNodeIdx2=basis%nodePositionIndexInv(nodePositionIndex(1),nodePositionIndex(2),nodePositionIndex(3),1)
                  IF(localNodeIdx1/=0.AND.localNodeIdx2/=0) THEN
                    IF(meshElements%elements(elementIdx)%meshElementNodes(localNodeIdx1)/= &
                      & meshElements%elements(elementIdx)%meshElementNodes(localNodeIdx2)) xiCollapsed=.TRUE.
                  ENDIF
                  nodePositionIndex(xiDirSearch)=nodePositionIndex(xiDirSearch)+1
                ENDDO !xiDirSearch
                IF(xiCollapsed) faceCollapsed(xiDirection)=.TRUE.
              ENDIF
            ENDDO !j
            nodePositionIndex(xiIdx)=numberOfNodesXiC(xiIdx)
          ENDDO !directionIndex
        ENDDO !xiIdx
        !Loop over the xi directions and calculate the surrounding elements
        DO xiIdx=1,basis%numberOfXi
          !Determine the xi directions that lie in this xi direction
          faceXi(1)=OTHER_XI_DIRECTIONS3(xiIdx,2,1)
          faceXi(2)=OTHER_XI_DIRECTIONS3(xiIdx,3,1)
          !Loop over the two faces
          DO directionIndex=-1,1,2
            xiDirection=directionIndex*xiIdx
            !Find nodes in the element on the appropriate face/line/point
            NULLIFY(nodeMatchList)
            CALL List_CreateStart(nodeMatchList,err,error,*999)
            CALL List_DataTypeSet(nodeMatchList,LIST_INTG_TYPE,err,error,*999)            
            CALL List_InitialSizeSet(nodeMatchList,16,err,error,*999)
            CALL List_CreateFinish(nodeMatchList,err,error,*999)
            IF(directionIndex==-1) THEN
              nodePositionIndex(xiIdx)=1
            ELSE
              nodePositionIndex(xiIdx)=numberOfNodesXiC(xiIdx)
            ENDIF
            !If the face is collapsed then don't look in this xi direction. The exception is if the opposite face is also
            !collpased. This may indicate that we have a funny element in non-rc coordinates that goes around the central
            !axis back to itself
            IF(faceCollapsed(xiDirection).AND..NOT.faceCollapsed(-xiDirection)) THEN
              !Do nothing - the match lists are already empty
            ELSE
              !Find the nodes to match and add them to the node match list
              DO localNodeIdx1=1,numberOfNodesXiC(faceXi(1))
                nodePositionIndex(faceXi(1))=localNodeIdx1
                DO localNodeIdx2=1,numberOfNodesXiC(faceXi(2))
                  nodePositionIndex(faceXi(2))=localNodeIdx2
                  localNodeIdx=BASIS%nodePositionIndexInv(nodePositionIndex(1),nodePositionIndex(2),nodePositionIndex(3),1)
                  IF(localNodeIdx/=0) THEN
                    nodeNumber=meshElements%elements(elementIdx)%meshElementNodes(localNodeIdx)
                    CALL List_ItemAdd(nodeMatchList,nodeNumber,err,error,*999)
                  ENDIF
                ENDDO !localNodeIdx2
              ENDDO !localNodeIdx1
            ENDIF
            CALL List_RemoveDuplicates(nodeMatchList,err,error,*999)
            CALL List_DetachAndDestroy(nodeMatchList,numberNodeMatches,nodeMatches,err,error,*999)
            numberSurrounding=0
            IF(numberNodeMatches>0) THEN
              !Find list of elements surrounding those nodes
              nodeNumber1=nodeMatches(1)
              DO surroundingElementIdx=1,meshNodes%nodes(nodeNumber1)%numberOfSurroundingElements
                elementIdx1=meshNodes%nodes(nodeNumber1)%surroundingElements(surroundingElementIdx)
                IF(elementIdx1/=elementIdx) THEN !Don't want the current element
                  ! grab the nodes list for current and this surrouding elements
                  ! current face : nodeMatches
                  ! candidate elem : meshElements%elements(elementIdx1)%meshElementNodes ! should this be globalElementNodes?
                  ! if all of current face belongs to the candidate element, we will have found the neighbour
                  CALL List_SubsetOf(nodeMatches(1:numberNodeMatches),meshElements%elements(elementIdx1)%meshElementNodes, &
                    & subset,err,error,*999)
                  IF(subset) THEN
                    CALL List_ItemAdd(adjacentElementsList(xiDirection)%ptr,elementIdx1,err,error,*999)
                    numberSurrounding=numberSurrounding+1
                  ENDIF
                ENDIF
              ENDDO !surroundingElementIdx
            ENDIF
            IF(ALLOCATED(nodeMatches)) DEALLOCATE(nodeMatches)
          ENDDO !directionIndex
        ENDDO !xiIdx
      CASE(BASIS_SIMPLEX_TYPE)
        !Loop over the xi coordinates and calculate the surrounding elements
        DO xiCoordIdx=1,basis%numberOfXiCoordinates
          !Find the other coordinates of the face/line/point
          faceXiC(1)=OTHER_XI_DIRECTIONS4(xiCoordIdx,1)
          faceXiC(2)=OTHER_XI_DIRECTIONS4(xiCoordIdx,2)
          faceXiC(3)=OTHER_XI_DIRECTIONS4(xiCoordIdx,3)
          !Find nodes in the element on the appropriate face/line/point
          NULLIFY(nodeMatchList)
          CALL List_CreateStart(nodeMatchList,err,error,*999)
          CALL List_DataTypeSet(nodeMatchList,LIST_INTG_TYPE,err,error,*999)
          CALL List_InitialSizeSet(nodeMatchList,16,err,error,*999)
          CALL List_CreateFinish(nodeMatchList,err,error,*999)
          nodePositionIndex(xiCoordIdx)=1 !Furtherest away from node with the xiCoordIdx'th coordinate
          !Find the nodes to match and add them to the node match list
          DO localNodeIdx1=1,numberOfNodesXiC(faceXiC(1))
            nodePositionIndex(faceXiC(1))=localNodeIdx1
            DO localNodeIdx2=1,numberOfNodesXiC(faceXiC(2))
              nodePositionIndex(faceXiC(2))=localNodeIdx2
              DO localNodeIdx3=1,numberOfNodesXiC(faceXiC(3))
                nodePositionIndex(faceXiC(3))=localNodeIdx3
                localNodeIdx=basis%nodePositionIndexInv(nodePositionIndex(1),nodePositionIndex(2),nodePositionIndex(3), &
                  nodePositionIndex(4))
                IF(localNodeIdx/=0) THEN
                  nodeNumber=meshTopology%ELEMENTS%elements(elementIdx)%meshElementNodes(localNodeIdx)
                  CALL List_ItemAdd(nodeMatchList,nodeNumber,err,error,*999)
                ENDIF
              ENDDO !localNodeIdx3
            ENDDO !localNodeIdx2
          ENDDO !localNodeIdx1
          CALL List_RemoveDuplicates(nodeMatchList,err,error,*999)
          CALL List_DetachAndDestroy(nodeMatchList,numberNodeMatches,nodeMatches,err,error,*999)
          IF(numberNodeMatches>0) THEN
            !Find list of elements surrounding those nodes
            DO nodeIdx=1,numberNodeMatches
              nodeNumber1=nodeMatches(nodeIdx)
              DO surroundingElementIdx=1,meshNodes%nodes(nodeNumber1)%numberOfSurroundingElements
                elementIdx1=meshNodes%nodes(nodeNumber1)%surroundingElements(surroundingElementIdx)
                IF(elementIdx1/=elementIdx) THEN !Don't want the current element
                  ! grab the nodes list for current and this surrouding elements
                  ! current face : nodeMatches
                  ! candidate elem : meshElements%elements(elementIdx1)%meshElementNodes 
                  ! if all of current face belongs to the candidate element, we will have found the neighbour
                  CALL List_SubsetOf(nodeMatches(1:numberNodeMatches),meshElements%elements(elementIdx1)%meshElementNodes, &
                    & subset,err,error,*999)
                  IF(subset) THEN
                    CALL List_ItemAdd(adjacentElementsList(xiCoordIdx)%ptr,elementIdx1,err,error,*999)
                  ENDIF
                ENDIF
              ENDDO !surroundingElementIdx
            ENDDO !nodeIdx
          ENDIF
          IF(ALLOCATED(nodeMatches)) DEALLOCATE(nodeMatches)
        ENDDO !xiCoordIdx
      CASE(BASIS_SERENDIPITY_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(BASIS_AUXILLIARY_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(BASIS_B_SPLINE_TP_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(BASIS_FOURIER_LAGRANGE_HERMITE_TP_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(BASIS_EXTENDED_LAGRANGE_TP_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The basis type of "//TRIM(NumberToVString(basis%type,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Set the surrounding elements for this element
      ALLOCATE(meshElements%elements(elementIdx)%adjacentElements(-basis%numberOfXiCoordinates:basis%numberOfXiCoordinates), &
        & STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate adjacent elements.",err,error,*999)
      DO xiCoordIdx=-basis%numberOfXiCoordinates,basis%numberOfXiCoordinates
        CALL Mesh_AdjacentElementInitialise(meshElements%elements(elementIdx)%adjacentElements(xiCoordIdx),err,error,*999)
        CALL List_RemoveDuplicates(adjacentElementsList(xiCoordIdx)%ptr,err,error,*999)
        CALL List_DetachAndDestroy(adjacentElementsList(xiCoordIdx)%ptr,meshElements%elements(elementIdx)% &
          & adjacentElements(xiCoordIdx)%numberOfAdjacentElements,adjacentElements,err,error,*999)
        ALLOCATE(meshElements%elements(elementIdx)%adjacentElements(xiCoordIdx)%adjacentElements(meshElements% &
          & elements(elementIdx)%adjacentElements(xiCoordIdx)%numberOfAdjacentElements),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate element adjacent elements.",err,error,*999)
        meshElements%elements(elementIdx)%adjacentElements(xiCoordIdx)%adjacentElements(1:meshElements%elements(elementIdx)% &
          & adjacentElements(xiCoordIdx)%numberOfAdjacentElements)=adjacentElements(1:meshElements%elements(elementIdx)% &
          & adjacentElements(xiCoordIdx)%numberOfAdjacentElements)
        IF(ALLOCATED(adjacentElements)) DEALLOCATE(adjacentElements)
      ENDDO !xiCoordIdx
    ENDDO !elementIdx           
    
    IF(diagnostics1) THEN
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Number of elements = ",meshElements%numberOfElements,err,error,*999)
      DO elementIdx=1,meshElements%numberOfElements
        basis=>meshElements%elements(elementIdx)%basis
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Global element number : ",elementIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of xi coordinates = ",basis%numberOfXiCoordinates, &
          & err,error,*999)
        DO xiCoordIdx=-basis%numberOfXiCoordinates,basis%numberOfXiCoordinates
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Xi coordinate : ",xiCoordIdx,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Number of adjacent elements = ", &
            & meshElements%elements(elementIdx)%adjacentElements(xiCoordIdx)%numberOfAdjacentElements,err,error,*999)
          IF(meshElements%elements(elementIdx)%adjacentElements(xiCoordIdx)%numberOfAdjacentElements>0) THEN
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,meshElements%elements(elementIdx)% &
              & adjacentElements(xiCoordIdx)%numberOfAdjacentElements,8,8,meshElements%elements(elementIdx)% &
              & adjacentElements(xiCoordIdx)% adjacentElements, &
              & '("        Adjacent elements :",8(X,I8))','(30x,8(X,I8))',err,error,*999)
          ENDIF
        ENDDO !xiCoordIdx
      ENDDO !elementIdx
    ENDIF
    
    EXITS("MeshTopology_ElementsAdjacentElementsCalculate")
    RETURN
999 IF(ALLOCATED(nodeMatches)) DEALLOCATE(nodeMatches)
    IF(ALLOCATED(adjacentElements)) DEALLOCATE(adjacentElements)
    IF(ASSOCIATED(nodeMatchList)) CALL List_Destroy(nodeMatchList,dummyErr,dummyError,*998)
998 DO xiCoordIdx=-4,4
      IF(ASSOCIATED(adjacentElementsList(xiCoordIdx)%ptr)) &
        & CALL List_Destroy(adjacentElementsList(xiCoordIdx)%ptr,dummyErr,dummyError,*997)
    ENDDO !xiCoordIdx
997 ERRORS("MeshTopology_ElementsAdjacentElementsCalculate",err,error)
    EXITS("MeshTopology_ElementsAdjacentElementsCalculate")
    RETURN 1
    
  END SUBROUTINE MeshTopology_ElementsAdjacentElementsCalculate
  
  !
  !================================================================================================================================
  !

  !>Finalises the elements data structures for a mesh topology and deallocates any memory. 
  SUBROUTINE MeshTopology_ElementsFinalise(meshElements,err,error,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: meshElements !<A pointer to the mesh topology to finalise the elements for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementIdx

    ENTERS("MeshTopology_ElementsFinalise",err,error,*999)

    IF(ASSOCIATED(meshElements)) THEN
      DO elementIdx=1,meshElements%numberOfElements
        CALL MeshElement_Finalise(meshElements%elements(elementIdx),err,error,*999)
      ENDDO !ne
      DEALLOCATE(meshElements%elements)
      IF(ASSOCIATED(meshElements%elementsTree)) CALL Tree_Destroy(meshElements%elementsTree,err,error,*999)
      DEALLOCATE(meshElements)
    ENDIF
 
    EXITS("MeshTopology_ElementsFinalise")
    RETURN
999 ERRORSEXITS("MeshTopology_ElementsFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_ElementsFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the elements in a given mesh topology. 
  SUBROUTINE MeshTopology_ElementsInitialise(meshTopology,err,error,*)

    !Argument variables
    TYPE(MeshTopologyType), POINTER :: meshTopology !<A pointer to the mesh topology to initialise the elements for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("MeshTopology_ElementsInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(meshTopology)) CALL FlagError("Mesh topology is not associated.",err,error,*998)
    IF(ASSOCIATED(meshTopology%elements)) CALL FlagError("Mesh already has topology elements associated.",err,error,*998)
    
    ALLOCATE(meshTopology%elements,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate mesh topology elements",err,error,*999)
    meshTopology%elements%numberOfElements=0
    meshTopology%elements%meshTopology=>meshTopology
    NULLIFY(meshTopology%elements%elementsTree)
    
    EXITS("MeshTopology_ElementsInitialise")
    RETURN
999 CALL MeshTopology_ElementsFinalise(meshTopology%elements,dummyErr,dummyError,*998)
998 ERRORSEXITS("MeshTopology_ElementsInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_ElementsInitialise
  
  !
  !================================================================================================================================
  !

  !>Finalises mesh data points and deallocates all memory.
  SUBROUTINE MeshTopology_DataPointsFinalise(meshDataPoints,err,error,*)

    !Argument variables
    TYPE(MeshDataPointsType), POINTER :: meshDataPoints !<A pointer to the mesh data points to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshTopology_DataPointsFinalise",err,error,*999)

    IF(ASSOCIATED(meshDataPoints)) THEN
      IF(ALLOCATED(meshDataPoints%dataPoints)) DEALLOCATE(meshDataPoints%dataPoints)
      IF(ALLOCATED(meshDataPoints%elementDataPoints)) DEALLOCATE(meshDataPoints%elementDataPoints)
      DEALLOCATE(meshDataPoints)
    ENDIF
    
    EXITS("MeshTopology_DataPointsFinalise")
    RETURN
999 ERRORSEXITS("MeshTopology_DataPointsFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_DataPointsFinalise
  
  !
  !================================================================================================================================
  !

  !>Initialises the elements in a given mesh topology.
  SUBROUTINE MeshTopology_DataPointsInitialise(meshTopology,err,error,*)

    !Argument variables
    TYPE(MeshTopologyType), POINTER :: meshTopology !<A pointer to the mesh topology to initialise the elements for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("MeshTopology_DataPointsInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(meshTopology)) CALL FlagError("Mesh topology is not associated.",err,error,*998)
    IF(ASSOCIATED(meshTopology%dataPoints)) CALL FlagError("Mesh already has topology data points associated",err,error,*998)
      
    ALLOCATE(meshTopology%dataPoints,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate mesh topology data points.",err,error,*999)
    meshTopology%dataPoints%totalNumberOfProjectedData=0
    meshTopology%dataPoints%meshTopology=>meshTopology
   
    EXITS("MeshTopology_DataPointsInitialise")
    RETURN
999 CALL MeshTopology_DataPointsFinalise(meshTopology%dataPoints,dummyErr,dummyError,*998)
998 ERRORSEXITS("MeshTopology_DataPointsInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_DataPointsInitialise
  
  !
  !================================================================================================================================
  !

  !>Returns the user element number for a global element number. \see OpenCMISS::Iron::cmfe_MeshElements_UserNumberGet
  SUBROUTINE MeshElements_ElementUserNumberGet(meshElements,globalNumber,userNumber,err,error,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: meshElements !<A pointer to the elements to set the user number for 
    INTEGER(INTG), INTENT(IN) :: globalNumber !<The global number of the elements to get.
    INTEGER(INTG), INTENT(OUT) :: userNumber !<On exit, The user number of the element to get
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("MeshElements_ElementUserNumberGet",err,error,*999)

    CALL MeshElements_AssertIsFinished(meshElements,err,error,*999)
    IF(globalNumber<1.OR.globalNumber>meshElements%numberOfElements) THEN          
      localError="The specified global element number of "//TRIM(NumberToVString(globalNumber,"*",err,error))// &
        & " is invalid. The global element number must be between 1 and "// &
        & TRIM(NumberToVString(meshElements%numberOfElements,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    userNumber=meshElements%elements(globalNumber)%userNumber
    
    EXITS("MeshElements_ElementUserNumberGet")
    RETURN
999 ERRORSEXITS("MeshElements_ElementUserNumberGet",err,error)
    RETURN 1
  
  END SUBROUTINE MeshElements_ElementUserNumberGet

  !
  !================================================================================================================================
  !

  !>Changes/sets the user number for a global element identified by a given global number. \see OpenCMISS::Iron::cmfe_MeshElements_UserNumberSet
  SUBROUTINE MeshElements_ElementUserNumberSet(meshElements,globalNumber,userNumber,err,error,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: meshElements !<A pointer to the elements to set the user number for
    INTEGER(INTG), INTENT(IN) :: globalNumber !<The global number of the elements to set.
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the element to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: globalElementNumber,insertStatus
    TYPE(TreeNodeType), POINTER :: treeNode
    TYPE(VARYING_STRING) :: localError

    ENTERS("MeshElements_ElementUserNumberSet",err,error,*999)

    CALL MeshElements_AssertNotFinished(meshElements,err,error,*999)
    IF(globalNumber<1.OR.globalNumber>meshElements%numberOfElements) THEN
      localError="The specified global element number of "//TRIM(NumberToVString(globalNumber,"*",err,error))// &
        & " is invalid. The global element number must be between 1 and "// &
        & TRIM(NumberToVString(meshElements%numberOfElements,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF    
    NULLIFY(treeNode)
    CALL Tree_Search(meshElements%elementsTree,userNumber,treeNode,err,error,*999)
    IF(ASSOCIATED(treeNode)) THEN
      CALL Tree_NodeValueGet(meshElements%elementsTree,treeNode,globalElementNumber,err,error,*999)
      localError="Element user number "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " is already used by global element number "// &
        & TRIM(NumberToVString(globalElementNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    CALL Tree_ItemDelete(meshElements%elementsTree,meshElements%elements(globalNumber)%userNumber,err,error,*999)
    CALL Tree_ItemInsert(meshElements%elementsTree,userNumber,globalNumber,insertStatus,err,error,*999)
    meshElements%elements(globalNumber)%userNumber=userNumber
   
    EXITS("MeshElements_ElementUserNumberSet")
    RETURN
999 ERRORSEXITS("MeshElements_ElementUserNumberSet",err,error)    
    RETURN 1
   
  END SUBROUTINE MeshElements_ElementUserNumberSet
  
  !
  !================================================================================================================================
  !

  !>Changes/sets the user numbers for all elements.
  SUBROUTINE MeshElements_ElementsUserNumbersAllSet(meshElements,userNumbers,err,error,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: meshElements !<A pointer to the elements to set all the user numbers for 
    INTEGER(INTG), INTENT(IN) :: userNumbers(:) !<The user numbers to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementIdx,insertStatus
    TYPE(TreeType), POINTER :: newElementsTree
    TYPE(VARYING_STRING) :: localError

    NULLIFY(newElementsTree)

    ENTERS("MeshElements_ElementsUserNumbersAllSet",err,error,*999)

    CALL MeshElements_AssertNotFinished(meshElements,err,error,*999)    
    IF(SIZE(userNumbers,1)/=meshElements%numberOfElements) THEN
      localError="The number of specified element user numbers of "// &
        TRIM(NumberToVstring(SIZE(userNumbers,1),"*",err,error))// &
        " does not match number of elements of "// &
        TRIM(NumberToVstring(meshElements%numberOfElements,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    !Check the users numbers to ensure that there are no duplicates          
    CALL Tree_CreateStart(newElementsTree,err,error,*999)
    CALL Tree_InsertTypeSet(newElementsTree,TREE_NO_DUPLICATES_ALLOWED,err,error,*999)
    CALL Tree_CreateFinish(newElementsTree,err,error,*999)
    DO elementIdx=1,meshElements%numberOfElements
      CALL Tree_ItemInsert(newElementsTree,userNumbers(elementIdx),elementIdx,insertStatus,err,error,*999)
      IF(insertStatus/=TREE_NODE_INSERT_SUCESSFUL) THEN
        localError="The specified user number of "//TRIM(NumberToVstring(userNumbers(elementIdx),"*",err,error))// &
          & " for global element number "//TRIM(NumberToVString(elementIdx,"*",err,error))// &
          & " is a duplicate. The user element numbers must be unique."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDDO !elementIdx
    CALL Tree_Destroy(meshElements%elementsTree,err,error,*999)
    meshElements%elementsTree=>newElementsTree
    NULLIFY(newElementsTree)
    DO elementIdx=1,meshElements%numberOfElements
      meshElements%elements(elementIdx)%globalNumber=elementIdx
      meshElements%elements(elementIdx)%userNumber=userNumbers(elementIdx)
    ENDDO !elementIdx
   
    EXITS("MeshElements_ElementsUserNumbersAllSet")
    RETURN
999 IF(ASSOCIATED(newElementsTree)) CALL Tree_Destroy(newElementsTree,err,error,*998)
998 ERRORSEXITS("MeshElements_ElementsUserNumbersAllSet",err,error)    
    RETURN 1
   
  END SUBROUTINE MeshElements_ElementsUserNumbersAllSet
  
  !
  !================================================================================================================================
  !

  !>Calculates the data points in the given mesh topology.
  SUBROUTINE MeshTopology_DataPointsCalculateProjection(mesh,dataProjection,err,error,*)
  
    !Argument variables
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh topology to calcualte the data projection for
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointIdx,elementIdx,exitTag,countIdx,projectionNumber,globalCountIdx,elementNumber
    TYPE(DataPointsType), POINTER :: dataPoints
    TYPE(DataProjectionResultType), POINTER :: dataProjectionResult
    TYPE(MeshDataPointsType), POINTER :: meshDataPoints
    TYPE(MeshElementsType), POINTER :: meshElements
    TYPE(MeshTopologyType), POINTER :: meshTopology
 
    ENTERS("MeshTopology_DataPointsCalculateProjection",err,error,*999)

    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Mesh is not associated.",err,error,*999)
    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    NULLIFY(dataPoints)
    CALL DataProjection_DataPointsGet(dataProjection,dataPoints,err,error,*999)
    NULLIFY(meshTopology)
    !Default the first mesh component topology to contain data points and elments!
    !\TODO: need to be changed once the data points topology is moved under meshTopologyType.    
    CALL Mesh_MeshTopologyGet(mesh,1,meshTopology,err,error,*999)
    NULLIFY(meshDataPoints)
    CALL MeshTopology_MeshDataPointsGet(meshTopology,meshDataPoints,err,error,*999)    
    NULLIFY(meshElements)
    CALL MeshTopology_MeshElementsGet(meshTopology,meshElements,err,error,*999)
    !Extract the global number of the data projection 
    projectionNumber=dataProjection%globalNumber
    ALLOCATE(meshDataPoints%elementDataPoints(meshElements%numberOfElements),STAT=err)     
    IF(err/=0) CALL FlagError("Could not allocate data points mesh topology element data points.",err,error,*999)
    DO elementIdx=1,meshElements%numberOfElements
      meshDataPoints%elementDataPoints(elementIdx)%elementNumber=meshElements%elements(elementIdx)%globalNumber
      meshDataPoints%elementDataPoints(elementIdx)%numberOfProjectedData=0
    ENDDO !elementIdx
    !Calculate number of projected data points on an element
    DO dataPointIdx=1,dataPoints%numberOfDataPoints
      dataProjectionResult=>dataProjection%dataProjectionResults(dataPointIdx)
      elementNumber=dataProjectionResult%elementNumber
      exitTag=dataProjectionResult%exitTag
      IF(exitTag/=DATA_PROJECTION_CANCELLED) THEN
        DO elementIdx=1,meshElements%numberOfElements
          IF(meshDataPoints%elementDataPoints(elementIdx)%elementNumber==elementNumber) &
            & meshDataPoints%elementDataPoints(elementIdx)%numberOfProjectedData= &
            & meshDataPoints%elementDataPoints(elementIdx)%numberOfProjectedData+1;
        ENDDO !elementIdx
      ENDIF
    ENDDO !dataPointIdx      
    !Allocate memory to store data indices and initialise them to be zero   
    DO elementIdx=1,meshElements%numberOfElements
      ALLOCATE(meshDataPoints%elementDataPoints(elementIdx)%dataIndices(meshDataPoints% &
        & elementDataPoints(elementIdx)%numberOfProjectedData),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate data points topology element data points.",err,error,*999)
      DO countIdx=1,meshDataPoints%elementDataPoints(elementIdx)%numberOfProjectedData
        meshDataPoints%elementDataPoints(elementIdx)%dataIndices(countIdx)%userNumber=0
        meshDataPoints%elementDataPoints(elementIdx)%dataIndices(countIdx)%globalNumber=0
      ENDDO !countIdx
    ENDDO !elementIdx
    !Record the indices of the data that projected on the elements 
    globalCountIdx=0
    meshDataPoints%totalNumberOfProjectedData=0
    DO dataPointIdx=1,dataPoints%numberOfDataPoints 
      dataProjectionResult=>dataProjection%dataProjectionResults(dataPointIdx)
      elementNumber=dataProjectionResult%elementNumber
      exitTag=dataProjectionResult%exitTag
      IF(exitTag/=DATA_PROJECTION_CANCELLED) THEN
        DO elementIdx=1,meshElements%numberOfElements
          countIdx=1         
          IF(meshDataPoints%elementDataPoints(elementIdx)%elementNumber==elementNumber) THEN
            globalCountIdx=globalCountIdx+1
            !Find the next data point index in this element
            DO WHILE(meshDataPoints%elementDataPoints(elementIdx)%dataIndices(countIdx)%globalNumber/=0)
              countIdx=countIdx+1
            ENDDO
            meshDataPoints%elementDataPoints(elementIdx)%dataIndices(countIdx)%userNumber=dataPointIdx
            meshDataPoints%elementDataPoints(elementIdx)%dataIndices(countIdx)%globalNumber=dataPointIdx !globalCountIdx (used this if only projected data are taken into account)
            meshDataPoints%totalNumberOfProjectedData=meshDataPoints%totalNumberOfProjectedData+1
          ENDIF
        ENDDO !elementIdx
      ENDIF
    ENDDO !dataPointIdx
    !Allocate memory to store total data indices in ascending order and element map
    ALLOCATE(meshDataPoints%dataPoints(meshDataPoints%totalNumberOfProjectedData),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate data points topology data points.",err,error,*999)
    !The global number for the data points will be looping through elements.
    countIdx=1  
    DO elementIdx=1,meshElements%numberOfElements
      DO dataPointIdx=1,meshDataPoints%elementDataPoints(elementIdx)%numberOfProjectedData
        meshDataPoints%dataPoints(countIdx)%userNumber=meshDataPoints%elementDataPoints(elementIdx)% &
          & dataIndices(dataPointIdx)%userNumber
        meshDataPoints%dataPoints(countIdx)%globalNumber=meshDataPoints%elementDataPoints(elementIdx)% &
          & dataIndices(dataPointIdx)%globalNumber
        meshDataPoints%dataPoints(countIdx)%elementNumber=meshDataPoints%elementDataPoints(elementIdx)% &
          & elementNumber
        countIdx=countIdx+1
      ENDDO !dataPointIdx
    ENDDO !elementIdx                      
    
    EXITS("MeshTopology_DataPointsCalculateProjection")
    RETURN
999 ERRORSEXITS("MeshTopology_DataPointsCalculateProjection",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_DataPointsCalculateProjection

  !
  !================================================================================================================================
  !

  !>Finalises the topology in the given mesh.
  SUBROUTINE Mesh_TopologiesFinalise(mesh,err,error,*)

    !Argument variables
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh to finalise the topology for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx

    ENTERS("Mesh_TopologiesFinalise",err,error,*999)

    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Mesh is not associated.",err,error,*999)

    IF(ALLOCATED(mesh%topology)) THEN
      DO componentIdx=1,SIZE(mesh%topology,1)
        CALL MeshTopology_Finalise(mesh%topology(componentIdx)%ptr,err,error,*999)
      ENDDO !componentIdx
      DEALLOCATE(mesh%topology)
    ENDIF
      
    EXITS("Mesh_TopologiesFinalise")
    RETURN
999 ERRORSEXITS("Mesh_TopologiesFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE Mesh_TopologiesFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the topology for a given mesh.
  SUBROUTINE Mesh_TopologiesInitialise(mesh,err,error,*)

    !Argument variables
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh to initialise the mesh topology for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("Mesh_TopologiesInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Mesh is not associated.",err,error,*998)
    IF(ALLOCATED(mesh%topology)) CALL FlagError("Mesh topology has already been allocated.",err,error,*998)
    
    !Allocate mesh topology
    ALLOCATE(mesh%topology(mesh%numberOfComponents),STAT=err)
    IF(err/=0) CALL FlagError("Mesh topology could not be allocated.",err,error,*999)
    DO componentIdx=1,mesh%numberOfComponents
      NULLIFY(mesh%topology(componentIdx)%ptr)
      CALL MeshTopology_Initialise(mesh%topology(componentIdx)%ptr,err,error,*999)
      mesh%topology(componentIdx)%ptr%mesh=>mesh
    ENDDO !componentIdx
    
    EXITS("Mesh_TopologiesInitialise")
    RETURN
999 CALL Mesh_TopologiesFinalise(mesh,dummyErr,dummyError,*998)
998 ERRORSEXITS("Mesh_TopologiesInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE Mesh_TopologiesInitialise
  
  !
  !================================================================================================================================
  !

  !>Finalises the topology in the given mesh. 
  SUBROUTINE MeshTopology_Finalise(meshTopology,err,error,*)

    !Argument variables
    TYPE(MeshTopologyType), POINTER :: meshTopology !<A pointer to the mesh topology to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshTopology_Finalise",err,error,*999)

    IF(ASSOCIATED(meshTopology)) THEN
      CALL MeshTopology_NodesFinalise(meshTopology%nodes,err,error,*999)
      CALL MeshTopology_ElementsFinalise(meshTopology%elements,err,error,*999)
      CALL MeshTopology_DofsFinalise(meshTopology%dofs,err,error,*999)
      CALL MeshTopology_DataPointsFinalise(meshTopology%dataPoints,err,error,*999)
      DEALLOCATE(meshTopology)
    ENDIF
 
    EXITS("MeshTopology_Finalise")
    RETURN
999 ERRORSEXITS("MeshTopology_Finalise",err,error)
    RETURN 1
   
  END SUBROUTINE MeshTopology_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the topology for a given mesh.
  SUBROUTINE MeshTopology_Initialise(meshTopology,err,error,*)

    !Argument variables
    TYPE(MeshTopologyType), POINTER :: meshTopology !<A pointer to the mesh topology to initialise. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("MeshTopology_Initialise",err,error,*998)

    IF(ASSOCIATED(meshTopology)) CALL FlagError("Mesh topology is already associated.",err,error,*998)
    
    !Allocate mesh topology
    ALLOCATE(meshTopology,STAT=err)
    IF(err/=0) CALL FlagError("Mesh topology could not be allocated.",err,error,*999)
    NULLIFY(meshTopology%elements)
    NULLIFY(meshTopology%nodes)
    NULLIFY(meshTopology%dofs)
    NULLIFY(meshTopology%dataPoints)
    !Initialise the topology components
    CALL MeshTopology_ElementsInitialise(meshTopology,err,error,*999)
    CALL MeshTopology_NodesInitialise(meshTopology,err,error,*999)
    CALL MeshTopology_DofsInitialise(meshTopology,err,error,*999)
    CALL MeshTopology_DataPointsInitialise(meshTopology,err,error,*999)
    
    EXITS("MeshTopology_Initialise")
    RETURN
999 CALL MeshTopology_Finalise(meshTopology,dummyErr,dummyError,*998)
998 ERRORSEXITS("MeshTopology_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_Initialise
  
  !
  !================================================================================================================================
  !
  
  !>Returns if the element in a mesh is on the boundary or not
  SUBROUTINE MeshElements_ElementOnBoundaryGet(meshElements,userElementNumber,onBoundary,err,error,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: meshElements !<A pointer to the mesh element containing the element to get the boundary type for
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number to get the boundary type for
    INTEGER(INTG), INTENT(OUT) :: onBoundary !<On return, the boundary type of the specified user element number.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: globalElementNumber

    ENTERS("MeshElements_ElementOnBoundaryGet",err,error,*999)

    IF(.NOT.ASSOCIATED(meshElements)) CALL FlagError("Mesh elements is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(meshElements%elements)) CALL FlagError("Mesh elements elements is not associated.",err,error,*999)
    
    CALL MeshElements_GlobalElementNumberGet(meshElements,userElementNumber,globalElementNumber,err,error,*999)
    IF(meshElements%elements(globalElementNumber)%boundaryElement) THEN
      onBoundary=MESH_ON_DOMAIN_BOUNDARY
    ELSE
      onBoundary=MESH_OFF_DOMAIN_BOUNDARY
    ENDIF
    
    EXITS("MeshElements_ElementOnBoundaryGet")
    RETURN
999 onBoundary=MESH_OFF_DOMAIN_BOUNDARY
    ERRORSEXITS("MeshElements_ElementOnBoundaryGet",err,error)    
    RETURN 1
   
  END SUBROUTINE MeshElements_ElementOnBoundaryGet
  
  !
  !================================================================================================================================
  !

  !>Finalises the given mesh topology node. 
  SUBROUTINE MeshTopology_NodeFinalise(node,err,error,*)

    !Argument variables
    TYPE(MeshNodeType) :: node !<The mesh node to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: derivativeIdx

    ENTERS("MeshTopology_NodeFinalise",err,error,*999)

    IF(ALLOCATED(node%derivatives)) THEN
      DO derivativeIdx=1,node%numberOfDerivatives
        CALL MeshTopology_NodeDerivativeFinalise(node%derivatives(derivativeIdx),err,error,*999)
      ENDDO !derivativeIdx
      DEALLOCATE(node%derivatives)
    ENDIF
    IF(ALLOCATED(node%surroundingElements)) DEALLOCATE(node%surroundingElements)
  
    EXITS("MeshTopology_NodeFinalise")
    RETURN
999 ERRORSEXITS("MeshTopology_NodeFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_NodeFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the given mesh topology node.
  SUBROUTINE MeshTopology_NodeInitialise(node,err,error,*)

    !Argument variables
    TYPE(MeshNodeType) :: node !<The mesh node to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshTopology_NodeInitialise",err,error,*999)

    node%userNumber=0
    node%globalNumber=0
    node%numberOfSurroundingElements=0
    node%numberOfDerivatives=0
    node%boundaryNode=.FALSE.
    
    EXITS("MeshTopology_NodeInitialise")
    RETURN
999 ERRORSEXITS("MeshTopology_NodeInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_NodeInitialise

  !
  !================================================================================================================================
  !

  !>Calculates the nodes used the mesh identified by a given mesh topology.
  SUBROUTINE MeshTopology_NodesCalculate(meshTopology,err,error,*)

    !Argument variables
    TYPE(MeshTopologyType), POINTER :: meshTopology !<A pointer to the mesh topology
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,elementIdx,insertStatus,localNodeIdx,globalNode,meshNodeIdx,meshNode,numberOfNodes 
    INTEGER(INTG), POINTER :: globalNodeNumbers(:)
    TYPE(BasisType), POINTER :: basis
    TYPE(MeshType), POINTER :: mesh
    TYPE(MeshElementsType), POINTER :: meshElements
    TYPE(MeshNodesType), POINTER :: meshNodes
    TYPE(NodesType), POINTER :: nodes
    TYPE(TreeType), POINTER :: globalNodesTree
    TYPE(TreeNodeType), POINTER :: treeNode
    TYPE(VARYING_STRING) :: dummyError,localError

    NULLIFY(globalNodeNumbers)
    NULLIFY(globalNodesTree)
    
    ENTERS("MeshTopology_NodesCalculate",err,error,*998)

    IF(.NOT.ASSOCIATED(meshTopology)) CALL FlagError("Mesh topology is not associated.",err,error,*999)
    NULLIFY(mesh)
    CALL MeshTopology_MeshGet(meshTopology,mesh,err,error,*999)
    NULLIFY(meshNodes)
    CALL MeshTopology_MeshNodesGet(meshTopology,meshNodes,err,error,*999)
    IF(ALLOCATED(meshNodes%nodes)) THEN
      localError="Mesh number "//TRIM(NumberToVString(mesh%userNumber,"*",err,error))// &
        & " already has allocated mesh topology nodes."
      CALL FlagError(localError,err,error,*998)
    ENDIF
    NULLIFY(meshElements)
    CALL MeshTopology_MeshElementsGet(meshTopology,meshElements,err,error,*999)
    NULLIFY(nodes)
    CALL Mesh_NodesGet(mesh,nodes,err,error,*999)
    !Work out what nodes are in the mesh
    CALL Tree_CreateStart(globalNodesTree,err,error,*999)
    CALL Tree_InsertTypeSet(globalNodesTree,TREE_NO_DUPLICATES_ALLOWED,err,error,*999)
    CALL Tree_CreateFinish(globalNodesTree,err,error,*999)
    DO elementIdx=1,meshElements%numberOfElements
      NULLIFY(basis)
      CALL MeshElements_BasisGet(meshElements,elementIdx,basis,err,error,*999)
      DO localNodeIdx=1,basis%numberOfNodes
        globalNode=meshElements%elements(elementIdx)%globalElementNodes(localNodeIdx)
        CALL Tree_ItemInsert(globalNodesTree,globalNode,globalNode,insertStatus,err,error,*999)
      ENDDO !localNodeIdx
    ENDDO !elementIdx
    CALL Tree_DetachAndDestroy(globalNodesTree,numberOfNodes,globalNodeNumbers,err,error,*999)
    !Set up the mesh nodes.
    ALLOCATE(meshNodes%nodes(numberOfNodes),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate mesh topology nodes nodes.",err,error,*999)
    CALL Tree_CreateStart(meshNodes%nodesTree,err,error,*999)
    CALL Tree_InsertTypeSet(meshNodes%nodesTree,TREE_NO_DUPLICATES_ALLOWED,err,error,*999)
    CALL Tree_CreateFinish(meshNodes%nodesTree,err,error,*999) 
    DO meshNodeIdx=1,numberOfNodes
      CALL MeshTopology_NodeInitialise(meshNodes%nodes(meshNodeIdx),err,error,*999)
      meshNodes%nodes(meshNodeIdx)%meshNumber=meshNodeIdx
      meshNodes%nodes(meshNodeIdx)%globalNumber=globalNodeNumbers(meshNodeIdx)
      meshNodes%nodes(meshNodeIdx)%userNumber=nodes%nodes(globalNodeNumbers(meshNodeIdx))%userNumber
      CALL Tree_ItemInsert(meshNodes%nodesTree,globalNodeNumbers(meshNodeIdx),meshNodeIdx,insertStatus,err,error,*999)
    ENDDO !nodeIdx
    meshNodes%numberOfNodes=numberOfNodes
    IF(ASSOCIATED(globalNodeNumbers)) DEALLOCATE(globalNodeNumbers)
    !Now recalculate the mesh element nodes
    DO elementIdx=1,meshElements%numberOfElements
      NULLIFY(basis)
      CALL MeshElements_BasisGet(meshElements,elementIdx,basis,err,error,*999)
      ALLOCATE(meshElements%elements(elementIdx)%meshElementNodes(basis%numberOfNodes),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate mesh topology elements mesh element nodes.",err,error,*999)
      DO localNodeIdx=1,basis%numberOfNodes
        globalNode=meshElements%elements(elementIdx)%globalElementNodes(localNodeIdx)
        NULLIFY(treeNode)
        CALL Tree_Search(meshNodes%nodesTree,globalNode,treeNode,err,error,*999)
        IF(ASSOCIATED(treeNode)) THEN
          CALL Tree_NodeValueGet(meshNodes%nodesTree,treeNode,meshNode,err,error,*999)
          meshElements%elements(elementIdx)%meshElementNodes(localNodeIdx)=meshNode
        ELSE
          localError="Could not find global node "//TRIM(NumberToVString(globalNode,"*",err,error))//" (user node "// &
            & TRIM(NumberToVString(nodes%nodes(globalNode)%userNumber,"*",err,error))//") in the mesh nodes."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDDO !localNodeIdx
    ENDDO !elementIdx
    
    IF(diagnostics1) THEN
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Number of mesh nodes = ",meshNodes%numberOfNodes,err,error,*999)
      DO meshNodeIdx=1,meshNodes%numberOfNodes
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Mesh node number = ",meshNodeIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Global node number = ",meshNodes%nodes(meshNodeIdx)%globalNumber, &
          & err,error,*999)        
      ENDDO !meshNodeIdx
    ENDIF
    
    EXITS("MeshTopology_NodesCalculate")
    RETURN
999 IF(ASSOCIATED(globalNodeNumbers)) DEALLOCATE(globalNodeNumbers)
    IF(ASSOCIATED(globalNodesTree)) CALL Tree_Destroy(globalNodesTree,dummyErr,dummyError,*998)
998 ERRORSEXITS("MeshTopology_NodesCalculate",err,error)
    RETURN 1
   
  END SUBROUTINE MeshTopology_NodesCalculate

  !
  !================================================================================================================================
  !

  !>Destroys the nodes in a mesh topology.
  SUBROUTINE MeshTopology_NodesDestroy(meshNodes,err,error,*)

    !Argument variables
    TYPE(MeshNodesType), POINTER :: meshNodes !<A pointer to the mesh nodes to destroy 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshTopology_NodesDestroy",err,error,*999)

    IF(.NOT.ASSOCIATED(meshNodes)) CALL FlagError("Mesh topology nodes is not associated.",err,error,*999)
    
    CALL MeshTopology_NodesFinalise(meshNodes,err,error,*999)

    EXITS("MeshTopology_NodesDestroy")
    RETURN
999 ERRORSEXITS("MeshTopology_NodesDestroy",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_NodesDestroy
  
  !
  !================================================================================================================================
  !

  !>Finalises the given mesh topology node. 
  SUBROUTINE MeshTopology_NodeDerivativeFinalise(meshNodeDerivative,err,error,*)

    !Argument variables
    TYPE(MeshNodeDerivativeType) :: meshNodeDerivative !<The mesh node derivative to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshTopology_NodeDerivativeFinalise",err,error,*999)

    IF(ALLOCATED(meshNodeDerivative%userVersionNumbers)) DEALLOCATE(meshNodeDerivative%userVersionNumbers)
    IF(ALLOCATED(meshNodeDerivative%dofIndex)) DEALLOCATE(meshNodeDerivative%dofIndex)

    EXITS("MeshTopology_NodeDerivativeFinalise")
    RETURN
999 ERRORSEXITS("MeshTopology_NodeDerivativeFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_NodeDerivativeFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the given mesh topology node.
  SUBROUTINE MeshTopology_NodeDerivativeInitialise(meshNodeDerivative,err,error,*)

    !Argument variables
    TYPE(MeshNodeDerivativeType) :: meshNodeDerivative !<The mesh node derivative to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshTopology_NodeDerivativeInitialise",err,error,*999)

    meshNodeDerivative%numberOfVersions=0
    meshNodeDerivative%globalDerivativeIndex=0
    meshNodederivative%partialDerivativeIndex=0

    EXITS("MeshTopology_NodeDerivativeInitialise")
    RETURN
999 ERRORSEXITS("MeshTopology_NodeDerivativeInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_NodeDerivativeInitialise

  !
  !================================================================================================================================
  !

  !>Calculates the number of derivatives at each node in a mesh topology.
  SUBROUTINE MeshTopology_NodesDerivativesCalculate(meshTopology,err,error,*)

    !Argument variables
    TYPE(MeshTopologyType), POINTER :: meshTopology !<A pointer to the mesh topology to calculate the derivates at each node for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: derivativeIdx,dummyErr,element,elementIdx,globalDerivative,localNodeIdx,maxNumberOfDerivatives,nodeIdx, &
      & nodeUserNumber,numberOfDerivatives
    INTEGER(INTG), ALLOCATABLE :: derivatives(:)
    LOGICAL :: found
    TYPE(BasisType), POINTER :: basis
    TYPE(ListType), POINTER :: nodeDerivativeList
    TYPE(MeshElementsType), POINTER :: meshElements
    TYPE(MeshNodesType), POINTER :: meshNodes
    TYPE(VARYING_STRING) :: dummyError,localError
    
    ENTERS("MeshTopology_NodesDerivativesCalculate",err,ERROR,*999)

    IF(.NOT.ASSOCIATED(meshTopology)) CALL FlagError("Mesh topology is not associated.",err,error,*999)
    
    NULLIFY(meshElements)
    CALL MeshTopology_MeshElementsGet(meshTopology,meshElements,err,error,*999)
    NULLIFY(meshNodes)
    CALL MeshTopology_MeshNodesGet(meshTopology,meshNodes,err,error,*999)
    
    !Loop over the mesh nodes
    DO nodeIdx=1,meshNodes%numberOfNodes
      !Calculate the number of derivatives and versions at each node. This needs to be calculated by looking at the
      !mesh elements as we may have an adjacent element in another domain with a higher order basis also with versions.
      CALL MeshNodes_UserNodeNumberGet(meshNodes,nodeIdx,nodeUserNumber,err,error,*999)
      NULLIFY(nodeDerivativeList)
      CALL List_CreateStart(nodeDerivativeList,err,error,*999)
      CALL List_DataTypeSet(nodeDerivativeList,LIST_INTG_TYPE,err,error,*999)
      CALL List_InitialSizeSet(nodeDerivativeList,8,err,error,*999)
      CALL List_CreateFinish(nodeDerivativeList,err,error,*999)
      maxNumberOfDerivatives=-1
      DO elementIdx=1,meshNodes%nodes(nodeIdx)%numberOfSurroundingElements
        element=meshNodes%nodes(nodeIdx)%surroundingElements(elementIdx)
        NULLIFY(basis)
        CALL MeshElements_BasisGet(meshElements,element,basis,err,error,*999)
        !Find the local node corresponding to this node
        found=.FALSE.
        DO localNodeIdx=1,basis%numberOfNodes
          IF(meshElements%elements(element)%meshElementNodes(localNodeIdx)==nodeIdx) THEN
            found=.TRUE.
            EXIT
          ENDIF
        ENDDO !localNodeIdx
        IF(.NOT.found) THEN
          localError="Could not find user node "//TRIM(NumberToVString(nodeUserNumber,"*",err,error))// &
            & " in the list of local nodes for element "//TRIM(NumberToVString(element,"*",err,error))// &
            & "."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        DO derivativeIdx=1,basis%numberOfDerivatives(localNodeIdx)
          CALL List_ItemAdd(nodeDerivativeList,basis%partialDerivativeIndex(derivativeIdx,localNodeIdx),err,error,*999)
        ENDDO !derivativeIdx
        IF(basis%numberOfDerivatives(localNodeIdx)>maxNumberOfDerivatives) &
          & maxNumberOfDerivatives=basis%numberOfDerivatives(localNodeidx)
      ENDDO !elem_idx
      CALL List_RemoveDuplicates(nodeDerivativeList,err,error,*999)
      CALL List_DetachAndDestroy(nodeDerivativeList,numberOfDerivatives,derivatives,err,error,*999)
      IF(numberOfDerivatives/=maxNumberOfDerivatives) THEN
         localError="Invalid mesh configuration. User node "//TRIM(NumberToVstring(nodeUserNumber,"*",err,error))// &
          & " has inconsistent derivative directions."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      !Set up the node derivatives.
      ALLOCATE(meshNodes%nodes(nodeIdx)%derivatives(maxNumberOfDerivatives),STAT=err)
      meshNodes%nodes(nodeIdx)%numberOfDerivatives=maxNumberOfDerivatives
      DO derivativeIdx=1,numberOfDerivatives
        CALL MeshTopology_NodeDerivativeInitialise(meshNodes%nodes(nodeIdx)%derivatives(derivativeIdx),err,error,*999)
        meshNodes%nodes(nodeIdx)%derivatives(derivativeIdx)%partialDerivativeIndex=derivatives(derivativeIdx)
        globalDerivative=PARTIAL_DERIVATIVE_GLOBAL_DERIVATIVE_MAP(derivatives(derivativeIdx))
        IF(globalDerivative==0) THEN
          localError="The partial derivative index of "//TRIM(NumberToVstring(derivatives(derivativeIdx),"*", &
            & err,error))//" for derivative number "//TRIM(NumberToVstring(derivativeIdx,"*",err,error))// &
            & " does not have a corresponding global derivative."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        meshNodes%nodes(nodeIdx)%derivatives(derivativeIdx)%globalDerivativeIndex=globalDerivative
      ENDDO !derivativeIdx
      DEALLOCATE(derivatives)
    ENDDO !nodeIdx
   
    IF(diagnostics1) THEN
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Number of mesh global nodes = ",meshNodes%numberOfNodes,err,error,*999)
      DO nodeIdx=1,meshNodes%numberOfNodes
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Mesh global node number = ",nodeIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of derivatives = ",meshNodes%nodes(nodeIdx)%numberOfDerivatives, &
          & err,error,*999)
        DO derivativeIdx=1,meshNodes%nodes(nodeIdx)%numberOfDerivatives
          !TODO: change output string below so that it writes out derivativeIdx index as well
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Global derivative index(derivativeIdx) = ", &
            & meshNodes%nodes(nodeIdx)%derivatives(derivativeIdx)%globalDerivativeIndex,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Partial derivative index(derivativeIdx) = ", &
            & meshNodes%nodes(nodeIdx)%derivatives(derivativeIdx)%partialDerivativeIndex,err,error,*999)
        ENDDO !derivativeIdx
      ENDDO !nodeIdx
    ENDIF
    
    EXITS("MeshTopology_NodesDerivativesCalculate")
    RETURN
999 IF(ALLOCATED(derivatives)) DEALLOCATE(derivatives)
    IF(ASSOCIATED(nodeDerivativeList)) CALL List_Destroy(nodeDerivativeList,dummyErr,dummyError,*998)
998 ERRORSEXITS("MeshTopology_NodesDerivativesCalculate",err,error)
    RETURN 1
   
  END SUBROUTINE MeshTopology_NodesDerivativesCalculate

  !
  !================================================================================================================================
  !

  !>Calculates the number of versions at each node in a topology.
  SUBROUTINE MeshTopology_NodesVersionCalculate(meshTopology,err,error,*)

    !Argument variables
    TYPE(MeshTopologyType), POINTER :: meshTopology !<A pointer to the mesh topology to calculate the versions at each node for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementIdx,localNodeIdx,derivativeIdx,nodeIdx,numberOfVersions,versionIdx
    INTEGER(INTG), ALLOCATABLE :: versions(:)
    TYPE(BasisType), POINTER :: basis
    TYPE(ListPtrType), POINTER :: nodeVersionList(:,:)
    TYPE(MeshElementsType), POINTER :: meshElements
    TYPE(MeshNodesType), POINTER :: meshNodes
    
    ENTERS("MeshTopology_NodesVersionCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(meshTopology)) CALL FlagError("Mesh topology is not associated.",err,error,*999)
    
    NULLIFY(meshElements)
    CALL MeshTopology_MeshElementsGet(meshTopology,meshElements,err,error,*999)
    NULLIFY(meshNodes)
    CALL MeshTopology_MeshNodesGet(meshTopology,meshNodes,err,error,*999)
    
    !Loop over the mesh elements
    !Calculate the number of versions at each node. This needs to be calculated by looking at all the mesh elements
    !as we may have an adjacent elements in another domain with a higher order basis along with different versions
    !being assigned to its derivatives.
    !\todo : See if there are any constraints that can be applied to restrict the amount of memory being allocated here
    ALLOCATE(nodeVersionList(MAXIMUM_GLOBAL_DERIV_NUMBER,meshNodes%numberOfNodes),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate node version list.",err,error,*999)
    DO nodeIdx=1,meshNodes%numberOfNodes
      DO derivativeIdx=1,meshNodes%nodes(nodeIdx)%numberOfDerivatives
        NULLIFY(nodeVersionList(derivativeIdx,nodeIdx)%ptr)
        CALL List_CreateStart(nodeVersionList(derivativeIdx,nodeIdx)%ptr,err,error,*999)
        CALL List_DataTypeSet(nodeVersionList(derivativeIdx,nodeIdx)%ptr,LIST_INTG_TYPE,err,error,*999)
        CALL List_InitialSizeSet(nodeVersionList(derivativeIdx,nodeIdx)%ptr,8,err,error,*999)
        CALL List_CreateFinish(nodeVersionList(derivativeIdx,nodeIdx)%ptr,err,error,*999)
      ENDDO!derivativeIdx
    ENDDO!nodeIdx
    DO elementIdx=1,meshElements%numberOfElements
      NULLIFY(basis)
      CALL MeshElements_BasisGet(meshElements,elementIdx,basis,err,error,*999)
      DO localNodeIdx=1,basis%numberOfNodes
        DO derivativeIdx=1,basis%numberOfDerivatives(localNodeIdx)
          CALL List_ItemAdd(nodeVersionList(derivativeIdx,meshElements%elements(elementIdx)%meshElementNodes(localNodeIdx))%ptr, &
            & meshElements%elements(elementIdx)%userElementNodeVersions(derivativeIdx,localNodeIdx),err,error,*999)
        ENDDO!derivativeIdx
      ENDDO!localNodeIdx
    ENDDO!elementIdx
    DO nodeIdx=1,meshNodes%numberOfNodes
      DO derivativeIdx=1,meshNodes%nodes(nodeIdx)%numberOfDerivatives
        CALL List_RemoveDuplicates(nodeVersionList(derivativeIdx,nodeIdx)%ptr,err,error,*999)
        CALL List_DetachAndDestroy(nodeVersionList(derivativeIdx,nodeIdx)%ptr,numberOfVersions,versions,err,error,*999)
        meshNodes%nodes(nodeIdx)%derivatives(derivativeIdx)%numberOfVersions = MAXVAL(versions(1:numberOfVersions))
        ALLOCATE(meshNodes%nodes(nodeIdx)%derivatives(derivativeIdx)%userVersionNumbers(meshNodes%nodes(nodeIdx)% &
          & derivatives(derivativeIdx)%numberOfVersions),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate node global derivative index.",err,error,*999)
        DO versionIdx=1,meshNodes%nodes(nodeIdx)%derivatives(derivativeIdx)%numberOfVersions 
          meshNodes%nodes(nodeIdx)%derivatives(derivativeIdx)%userVersionNumbers(versionIdx)=versionIdx
        ENDDO !versionIdx
        DEALLOCATE(versions)
      ENDDO!derivativeIdx
    ENDDO!nodeIdx
    DEALLOCATE(nodeVersionList)
   
    IF(diagnostics1) THEN
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Number of mesh global nodes = ",meshNodes%numberOfNodes,err,error,*999)
      DO nodeIdx=1,meshNodes%numberOfNodes
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Mesh global node number = ",nodeIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of derivatives = ", &
          & meshNodes%nodes(nodeIdx)%numberOfDerivatives,err,error,*999)
        DO derivativeIdx=1,meshNodes%nodes(nodeIdx)%numberOfDerivatives
          !\todo : change output string below so that it writes out derivativeIdx index as well
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Global derivative index(derivativeIdx) = ", &
            & meshNodes%nodes(nodeIdx)%derivatives(derivativeIdx)%globalDerivativeIndex,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Partial derivative index(derivativeIdx) = ", &
            & meshNodes%nodes(nodeIdx)%derivatives(derivativeIdx)%partialDerivativeIndex,err,error,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1, &
            & meshNodes%nodes(nodeIdx)%derivatives(derivativeIdx)%numberOfVersions,8,8, &
            & meshNodes%nodes(nodeIdx)%derivatives(derivativeIdx)%userVersionNumbers, & 
            & '("    User Version index(derivativeIdx,:) :",8(X,I2))','(36X,8(X,I2))',err,error,*999)
        ENDDO !derivativeIdx
      ENDDO !nodeIdx
    ENDIF
    
    EXITS("MeshTopology_NodesVersionCalculate")
    RETURN
999 IF(ALLOCATED(versions)) DEALLOCATE(versions)
    IF(ASSOCIATED(nodeVersionList)) THEN
      DO nodeIdx=1,SIZE(nodeVersionList,1)
        DO derivativeIdx=1,SIZE(nodeVersionList,2)
          CALL List_Destroy(nodeVersionList(derivativeIdx,nodeIdx)%ptr,err,error,*998)
        ENDDO !derivativeIdx
      ENDDO !nodeIdx
      DEALLOCATE(nodeVersionList)
    ENDIF
998 ERRORSEXITS("MeshTopology_NodesVersionCalculate",err,error)
    RETURN 1
   
  END SUBROUTINE MeshTopology_NodesVersionCalculate

  !
  !================================================================================================================================
  !

  !>Calculates the element numbers surrounding a node for a mesh.
  SUBROUTINE MeshTopology_SurroundingElementsCalculate(meshTopology,err,error,*)

    !Argument variables
    TYPE(MeshTopologyType), POINTER :: meshTopology !<A pointer to the mesh topology to calculate the elements surrounding each node for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: element,elementIdx,insertPosition,localNodeIdx,node,surroundingElementNumber
    INTEGER(INTG), ALLOCATABLE :: newSurroundingElements(:)
    LOGICAL :: foundElement
    TYPE(BasisType), POINTER :: basis
    TYPE(MeshElementsType), POINTER :: meshElements
    TYPE(MeshNodesType), POINTER :: meshNodes

    ENTERS("MeshTopology_SurroundingElementsCalculate",err,error,*999)
    
    IF(.NOT.ASSOCIATED(meshTopology)) CALL FlagError("Mesh topology is not associated.",err,error,*999)
    
    NULLIFY(meshElements)
    CALL MeshTopology_MeshElementsGet(meshTopology,meshElements,err,error,*999)
    NULLIFY(meshNodes)
    CALL MeshTopology_MeshNodesGet(meshTopology,meshNodes,err,error,*999)
    IF(.NOT.ALLOCATED(meshNodes%nodes)) CALL FlagError("Mesh topology nodes nodes have not been allocated.",err,error,*999)

    DO elementIdx=1,meshElements%numberOfElements
      NULLIFY(basis)
      CALL MeshElements_BasisGet(meshElements,elementIdx,basis,err,error,*999)
      DO localNodeIdx=1,basis%numberOfNodes
        node=meshElements%elements(elementIdx)%meshElementNodes(localNodeIdx)
        foundElement=.FALSE.
        element=1
        insertPosition=1
        DO WHILE(element<=meshNodes%nodes(node)%numberOfSurroundingElements.AND..NOT.foundElement)
          surroundingElementNumber=meshNodes%nodes(node)%surroundingElements(element)
          IF(surroundingElementNumber==elementIdx) foundElement=.TRUE.
          element=element+1
          IF(elementIdx>=surroundingElementNumber) insertPosition=element
        ENDDO
        IF(.NOT.foundElement) THEN
          !Insert element into surrounding elements
          ALLOCATE(newSurroundingElements(meshNodes%nodes(node)%numberOfSurroundingElements+1),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate new surrounding elements.",err,error,*999)
          IF(ALLOCATED(meshNodes%nodes(node)%surroundingElements)) THEN
            newSurroundingElements(1:insertPosition-1)=meshNodes%nodes(node)%surroundingElements(1:insertPosition-1)
            newSurroundingElements(insertPosition)=elementIdx
            newSurroundingElements(insertPosition+1:meshNodes%nodes(node)%numberOfSurroundingElements+1)= &
              & meshNodes%nodes(node)%surroundingElements(insertPosition:meshNodes%nodes(node)%numberOfSurroundingElements)
          ELSE
            newSurroundingElements(1)=elementIdx
          ENDIF
          CALL MOVE_ALLOC(newSurroundingElements,meshNodes%nodes(node)%surroundingElements)
          meshNodes%nodes(node)%numberOfSurroundingElements=meshNodes%nodes(node)%numberOfSurroundingElements+1
        ENDIF
      ENDDO !localNodeIdx
    ENDDO !elementIdx

    EXITS("MeshTopology_SurroundingElementsCalculate")
    RETURN
999 IF(ALLOCATED(newSurroundingElements)) DEALLOCATE(newSurroundingElements)
    ERRORSEXITS("MeshTopology_SurroundingElementsCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_SurroundingElementsCalculate
  
  !
  !===============================================================================================================================
  !

  !>Finalises the nodes data structures for a mesh topology and deallocates any memory.
  SUBROUTINE MeshTopology_NodesFinalise(meshNodes,err,error,*)

    !Argument variables
    TYPE(MeshNodesType), POINTER :: meshNodes !<A pointer to the mesh topology nodes to finalise the nodes for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: nodeIdx

    ENTERS("MeshTopology_NodesFinalise",err,error,*999)

    IF(ASSOCIATED(meshNodes)) THEN
      IF(ALLOCATED(meshNodes%nodes)) THEN
        DO nodeIdx=1,SIZE(meshNodes%nodes,1)
          CALL MeshTopology_NodeFinalise(meshNodes%nodes(nodeIdx),err,error,*999)
        ENDDO !nodesIdx
        DEALLOCATE(meshNodes%nodes)
      ENDIF
      IF(ASSOCIATED(meshNodes%nodesTree)) CALL Tree_Destroy(meshNodes%nodesTree,err,error,*999)
      DEALLOCATE(meshNodes)
    ENDIF
 
    EXITS("MeshTopology_NodesFinalise")
    RETURN
999 ERRORSEXITS("MeshTopology_NodesFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_NodesFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the nodes in a given mesh topology. \todo finalise on errors
  SUBROUTINE MeshTopology_NodesInitialise(meshTopology,err,error,*)

    !Argument variables
    TYPE(MeshTopologyType), POINTER :: meshTopology !<A pointer to the mesh topology to initialise the nodes for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("MeshTopology_NodesInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(meshTopology)) CALL FlagError("Mesh topology is not associated.",err,error,*998)
    IF(ASSOCIATED(meshTopology%nodes)) CALL FlagError("Mesh already has topology nodes associated.",err,error,*998)
    
    ALLOCATE(meshTopology%nodes,STAT=err)    
    IF(err/=0) CALL FlagError("Could not allocate topology nodes.",err,error,*999)
    meshTopology%nodes%numberOfNodes=0
    meshTopology%nodes%meshTopology=>meshTopology
    NULLIFY(meshTopology%nodes%nodesTree)
   
    EXITS("MeshTopology_NodesInitialise")
    RETURN
999 CALL MeshTopology_NodesFinalise(meshTopology%nodes,dummyErr,dummyError,*998)
998 ERRORSEXITS("MeshTopology_NodesInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_NodesInitialise

  !
  !================================================================================================================================
  !

  !>Finalises the meshes and deallocates all memory
  SUBROUTINE Meshes_Finalise(meshes,err,error,*)

   !Argument variables
    TYPE(MeshesType), POINTER :: meshes !<A pointer to the meshes to finalise the meshes for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(MeshType), POINTER :: mesh
 
    ENTERS("Meshes_Finalise",err,error,*999)

    IF(ASSOCIATED(meshes)) THEN
      DO WHILE(meshes%numberOfMeshes>0)
        mesh=>meshes%meshes(1)%ptr
        CALL Mesh_Destroy(mesh,err,error,*999)
      ENDDO !meshIdx
      DEALLOCATE(meshes)
    ENDIF
 
    EXITS("Meshes_Finalise")
    RETURN
999 ERRORSEXITS("Meshes_Finalise",err,error)
    RETURN 1
   
  END SUBROUTINE Meshes_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the generic meshes.
  SUBROUTINE Meshes_InitialiseGeneric(meshes,err,error,*)

    !Argument variables
    TYPE(MeshesType), POINTER :: meshes !<A pointer to the meshes to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("Meshes_InitialiseGeneric",err,error,*998)

    IF(ASSOCIATED(meshes)) CALL FlagError("Meshes is already associated.",err,error,*998)
    
    ALLOCATE(meshes,STAT=err)
    IF(err/=0) CALL FlagError("Meshes could not be allocated",err,error,*999)
    NULLIFY(meshes%region)
    NULLIFY(meshes%interface)
    meshes%numberOfMeshes=0
    
    EXITS("Meshes_InitialiseGeneric")
    RETURN
999 CALL Meshes_Finalise(meshes,dummyErr,dummyError,*998)
998 ERRORSEXITS("Meshes_InitialiseGeneric",err,error)
    RETURN 1
    
  END SUBROUTINE Meshes_InitialiseGeneric

  !
  !================================================================================================================================
  !

  !>Initialises the meshes for the given interface.
  SUBROUTINE Meshes_InitialiseInterface(interface,err,error,*)

    !Argument variables
    TYPE(InterfaceType), POINTER :: interface !<A pointer to the interface to initialise the meshes for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Meshes_InitialiseInterface",err,error,*999)

    IF(.NOT.ASSOCIATED(interface))  CALL FlagError("Interface is not associated.",err,error,*999)
    IF(ASSOCIATED(interface%meshes)) THEN
      localError="Interface number "//TRIM(NumberToVString(interface%userNumber,"*",err,error))// &
        & " already has a mesh associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    CALL Meshes_InitialiseGeneric(interface%meshes,err,error,*999)
    interface%meshes%interface=>interface
     
    EXITS("Meshes_InitialiseInterface")
    RETURN
999 ERRORSEXITS("Meshes_InitialiseInterface",err,error)
    RETURN 1
    
  END SUBROUTINE Meshes_InitialiseInterface

  !
  !================================================================================================================================
  !

  !>Initialises the meshes for the given region.
  SUBROUTINE Meshes_InitialiseRegion(region,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: region !<A pointer to the region to initialise the meshes for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Meshes_InitialiseRegion",err,error,*999)

    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated.",err,error,*999)
    IF(ASSOCIATED(region%meshes)) THEN
      localError="Region number "//TRIM(NumberToVString(region%userNumber,"*",err,error))// &
        & " already has a mesh associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    CALL Meshes_InitialiseGeneric(region%meshes,err,error,*999)
    region%meshes%region=>region
   
    EXITS("Meshes_InitialiseRegion")
    RETURN
999 ERRORSEXITS("Meshes_InitialiseRegion",err,error)
    RETURN 1
    
  END SUBROUTINE Meshes_InitialiseRegion

  !
  !================================================================================================================================
  !

  !>Initialises the embedded meshes.
  SUBROUTINE EMBEDDED_MESH_INITIALISE(MESH_EMBEDDING,err,error,*)

    !Argument variables
    !TYPE(MESH_EMBEDDING_TYPE), INTENT(INOUT) :: MESH_EMBEDDING !<Mesh embedding to initialise
    TYPE(MESH_EMBEDDING_TYPE), POINTER :: MESH_EMBEDDING !<Mesh embedding to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EMBEDDED_MESH_INITIALISE",err,error,*998)
    
    ALLOCATE(MESH_EMBEDDING,STAT=err)
    NULLIFY(MESH_EMBEDDING%PARENT_MESH)
    NULLIFY(MESH_EMBEDDING%CHILD_MESH)
    
    EXITS("EMBEDDED_MESH_INITIALISE")
    RETURN
!999 CALL EMBEDDED_Mesh_Finalise(MESH_EMBEDDING,dummyErr,dummyError,*998)
998 ERRORSEXITS("EMBEDDED_MESH_INITIALISE",err,error)
    RETURN 1
  END SUBROUTINE EMBEDDED_MESH_INITIALISE

  !
  !================================================================================================================================
  !

  !>Creates an embedding of one mesh in another
  SUBROUTINE MESH_EMBEDDING_CREATE(MESH_EMBEDDING, PARENT_MESH, CHILD_MESH,err,error,*)
!    TYPE(MESH_EMBEDDING_TYPE), INTENT(INOUT) :: MESH_EMBEDDING !<Mesh embedding to create.
    TYPE(MESH_EMBEDDING_TYPE), POINTER :: MESH_EMBEDDING !<Mesh embedding to create.
    TYPE(MeshType), POINTER, INTENT(IN) :: PARENT_MESH !<The parent mesh
    TYPE(MeshType), POINTER, INTENT(IN) :: CHILD_MESH  !<The child mesh
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: NGP = 0, ne

    ENTERS("MESH_EMBEDDING_CREATE",err,error,*999) 
    
    WRITE(*,*) 'parent mesh', PARENT_MESH%numberOfElements
    WRITE(*,*) 'child mesh', child_MESH%numberOfElements
    CALL EMBEDDED_MESH_INITIALISE(MESH_EMBEDDING,err,error,*999)

    DO ne=1,PARENT_MESH%numberOfElements
      NGP = MAX(NGP,PARENT_MESH%TOPOLOGY(1)%ptr%ELEMENTS%elements(ne)%BASIS%QUADRATURE%&
        & QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr%NUMBER_OF_GAUSS)
    ENDDO !ne

    MESH_EMBEDDING%PARENT_MESH => PARENT_MESH
    MESH_EMBEDDING%CHILD_MESH  => CHILD_MESH
    ALLOCATE(MESH_EMBEDDING%CHILD_NODE_XI_POSITION(PARENT_MESH%numberOfElements),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate child node positions.",err,error,*999)
    ALLOCATE(MESH_EMBEDDING%GAUSS_POINT_XI_POSITION(NGP,PARENT_MESH%numberOfElements),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate gauss point positions.",err,error,*999)
    
    EXITS("MESH_EMBEDDING_CREATE")
    RETURN 

999 ERRORSEXITS("MESH_EMBEDDING_CREATE",err,error)
    RETURN 1
  END SUBROUTINE MESH_EMBEDDING_CREATE

  !
  !================================================================================================================================
  !

  !>Sets the positions of nodes in the child mesh for one element in the parent mesh
  SUBROUTINE MESH_EMBEDDING_SET_CHILD_NODE_POSITION(MESH_EMBEDDING, ELEMENT_NUMBER, NODE_NUMBERS, XI_COORDS,err,error,*)
    TYPE(MESH_EMBEDDING_TYPE), INTENT(INOUT) :: MESH_EMBEDDING !<The mesh embedding object
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER  !<Element number in the parent mesh
    INTEGER(INTG), INTENT(IN) :: NODE_NUMBERS(:) !<NODE_NUMBERS(nodeIdx) Node numbers in child mesh for the nodeIdx'th embedded node in the ELEMENT_NUMBER'th element of the parent mesh
    REAL(DP), INTENT(IN)      :: XI_COORDS(:,:)  !<XI_COORDS(:,nodeIdx) Xi coordinates of the nodeIdx'th embedded node in the ELEMENT_NUMBER'th

    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("MESH_EMBEDDING_SET_CHILD_NODE_POSITION",err,error,*999)

    IF(ELEMENT_NUMBER<1 .OR. ELEMENT_NUMBER > MESH_EMBEDDING%PARENT_MESH%numberOfElements) THEN
      CALL FlagError("Element number out of range",err,error,*999)
    ENDIF

    MESH_EMBEDDING%CHILD_NODE_XI_POSITION(ELEMENT_NUMBER)%NUMBER_OF_NODES = SIZE(NODE_NUMBERS)

    ALLOCATE(MESH_EMBEDDING%CHILD_NODE_XI_POSITION(ELEMENT_NUMBER)%NODE_NUMBERS(SIZE(NODE_NUMBERS)))
    MESH_EMBEDDING%CHILD_NODE_XI_POSITION(ELEMENT_NUMBER)%NODE_NUMBERS(1:SIZE(NODE_NUMBERS)) = NODE_NUMBERS(1:SIZE(NODE_NUMBERS))

    ALLOCATE(MESH_EMBEDDING%CHILD_NODE_XI_POSITION(ELEMENT_NUMBER)%XI_COORDS(SIZE(XI_COORDS,1),SIZE(XI_COORDS,2)))
    MESH_EMBEDDING%CHILD_NODE_XI_POSITION(ELEMENT_NUMBER)%XI_COORDS(1:SIZE(XI_COORDS,1),1:SIZE(XI_COORDS,2)) = &
      & XI_COORDS(1:SIZE(XI_COORDS,1),1:SIZE(XI_COORDS,2))
    
    EXITS("MESH_EMBEDDING_SET_CHILD_NODE_POSITION")
    RETURN
999 ERRORSEXITS("MESH_EMBEDDING_SET_CHILD_NODE_POSITION",err,error)
    RETURN 1
  END SUBROUTINE MESH_EMBEDDING_SET_CHILD_NODE_POSITION

  !
  !================================================================================================================================
  !

  !>Sets the positions of a Gauss point of the parent mesh in terms of element/xi coordinate in the child mesh
  SUBROUTINE MESH_EMBEDDING_SET_GAUSS_POINT_DATA(MESH_EMBEDDING, PARENT_ELEMENT_NUMBER, GAUSSPT_NUMBER,&
    & PARENT_XI_COORD, CHILD_ELEMENT_NUMBER, CHILD_XI_COORD,err,error,*)
    TYPE(MESH_EMBEDDING_TYPE), INTENT(INOUT) :: MESH_EMBEDDING   !<The mesh embedding object
    INTEGER(INTG), INTENT(IN) :: PARENT_ELEMENT_NUMBER           !<Element number in the parent mesh
    INTEGER(INTG), INTENT(IN) :: GAUSSPT_NUMBER                  !<Gauss point number in this element
    REAL(DP), INTENT(IN) :: PARENT_XI_COORD(:)              !<Xi coordinate in parent element

    INTEGER(INTG), INTENT(IN) :: CHILD_ELEMENT_NUMBER !<Element number in the child mesh
    REAL(DP), INTENT(IN) :: CHILD_XI_COORD(:)    !<Xi coordinate in child element

    INTEGER(INTG), INTENT(OUT) :: err           !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error  !<The error string

    ENTERS("MESH_EMBEDDING_SET_GAUSS_POINT_DATA",err,error,*999)

    IF(PARENT_ELEMENT_NUMBER<1 .OR. PARENT_ELEMENT_NUMBER > MESH_EMBEDDING%PARENT_MESH%numberOfElements) THEN
      CALL FlagError("Parent element number out of range",err,error,*999)
    ENDIF
    IF(CHILD_ELEMENT_NUMBER<1 .OR. CHILD_ELEMENT_NUMBER > MESH_EMBEDDING%CHILD_MESH%numberOfElements) THEN
      CALL FlagError("Child element number out of range",err,error,*999)
    ENDIF
    IF(GAUSSPT_NUMBER<1 .OR. GAUSSPT_NUMBER > SIZE(MESH_EMBEDDING%GAUSS_POINT_XI_POSITION,1)) THEN
      CALL FlagError("Gauss point number out of range",err,error,*999)
    ENDIF

    ALLOCATE(MESH_EMBEDDING%GAUSS_POINT_XI_POSITION(GAUSSPT_NUMBER,PARENT_ELEMENT_NUMBER)&
     & %PARENT_XI_COORD(SIZE(PARENT_XI_COORD)))
    ALLOCATE(MESH_EMBEDDING%GAUSS_POINT_XI_POSITION(GAUSSPT_NUMBER,PARENT_ELEMENT_NUMBER)&
     & %CHILD_XI_COORD(SIZE(CHILD_XI_COORD)))


    MESH_EMBEDDING%GAUSS_POINT_XI_POSITION(GAUSSPT_NUMBER,PARENT_ELEMENT_NUMBER)%PARENT_XI_COORD(1:SIZE(PARENT_XI_COORD)) = &
      & PARENT_XI_COORD(1:SIZE(PARENT_XI_COORD))
    MESH_EMBEDDING%GAUSS_POINT_XI_POSITION(GAUSSPT_NUMBER,PARENT_ELEMENT_NUMBER)%CHILD_XI_COORD(1:SIZE(CHILD_XI_COORD)) = &
      & CHILD_XI_COORD(1:SIZE(CHILD_XI_COORD))
    MESH_EMBEDDING%GAUSS_POINT_XI_POSITION(GAUSSPT_NUMBER,PARENT_ELEMENT_NUMBER)%ELEMENT_NUMBER = CHILD_ELEMENT_NUMBER

    EXITS("MESH_EMBEDDING_SET_GAUSS_POINT_DATA")
    RETURN
999 ERRORSEXITS("MESH_EMBEDDING_SET_GAUSS_POINT_DATA",err,error)
    RETURN 1
  END SUBROUTINE MESH_EMBEDDING_SET_GAUSS_POINT_DATA


  !
  !================================================================================================================================
  !

END MODULE MeshRoutines


