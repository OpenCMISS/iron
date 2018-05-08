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
  USE DataProjectionAccessRoutines
  USE DecompositionRoutines
  USE Kinds
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE Lists
  USE MeshAccessRoutines
#ifndef NOMPIMOD
  USE MPI
#endif
  USE NodeRoutines
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
  
  !>Starts the process of creating a mesh
  INTERFACE MESH_CREATE_START
    MODULE PROCEDURE MESH_CREATE_START_INTERFACE
    MODULE PROCEDURE MESH_CREATE_START_REGION
  END INTERFACE !MESH_CREATE_START

  !>Initialises the meshes for a region or interface.
  INTERFACE MESHES_INITIALISE
    MODULE PROCEDURE MESHES_INITIALISE_INTERFACE
    MODULE PROCEDURE MESHES_INITIALISE_REGION
  END INTERFACE !MESHES_INITIALISE

  INTERFACE MeshTopology_ElementCheckExists
    MODULE PROCEDURE MeshTopology_ElementCheckExistsMesh
    MODULE PROCEDURE MeshTopology_ElementCheckExistsMeshElements
  END INTERFACE MeshTopology_ElementCheckExists
  
  INTERFACE MeshTopology_ElementGet
    MODULE PROCEDURE MeshTopology_ElementGetMeshElements
  END INTERFACE MeshTopology_ElementGet
  
  INTERFACE MeshTopology_NodeCheckExists
    MODULE PROCEDURE MeshTopology_NodeCheckExistsMesh
    MODULE PROCEDURE MeshTopology_NodeCheckExistsMeshNodes    
  END INTERFACE MeshTopology_NodeCheckExists  
  
  INTERFACE MeshTopology_NodeGet
    MODULE PROCEDURE MeshTopology_NodeGetMeshNodes    
  END INTERFACE MeshTopology_NodeGet
  
  PUBLIC MESH_ON_DOMAIN_BOUNDARY,MESH_OFF_DOMAIN_BOUNDARY

  PUBLIC MeshTopology_ElementCheckExists,MeshTopology_NodeCheckExists
  
  PUBLIC MESH_CREATE_START,MESH_CREATE_FINISH

  PUBLIC MESH_DESTROY
  
  PUBLIC MESH_NUMBER_OF_COMPONENTS_GET,MESH_NUMBER_OF_COMPONENTS_SET

  PUBLIC MESH_NUMBER_OF_ELEMENTS_GET,MESH_NUMBER_OF_ELEMENTS_SET

  PUBLIC MESH_TOPOLOGY_ELEMENTS_CREATE_START,MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH

  PUBLIC MESH_TOPOLOGY_ELEMENTS_DESTROY

  PUBLIC MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_GET,MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET

  PUBLIC MESH_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENT_GET

  PUBLIC MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_GET

  PUBLIC MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET,MeshElements_ElementNodeVersionSet

  PUBLIC MeshTopology_DataPointsCalculateProjection

  PUBLIC MeshTopology_ElementOnBoundaryGet

  PUBLIC MeshElements_ElementUserNumberGet,MeshElements_ElementUserNumberSet
  
  PUBLIC MeshTopology_ElementsUserNumbersAllSet

  PUBLIC MeshTopology_NodeOnBoundaryGet

  PUBLIC MeshTopology_NodeDerivativesGet

  PUBLIC MeshTopology_NodeNumberOfDerivativesGet

  PUBLIC MeshTopology_NodeNumberOfVersionsGet

  PUBLIC MeshTopology_NodesNumberOfNodesGet

  PUBLIC MeshTopology_NodesDestroy
  
  PUBLIC MESH_SURROUNDING_ELEMENTS_CALCULATE_SET

  PUBLIC MESH_EMBEDDING_CREATE,MESH_EMBEDDING_SET_CHILD_NODE_POSITION

  PUBLIC MESH_EMBEDDING_SET_GAUSS_POINT_DATA

  PUBLIC MESHES_INITIALISE,MESHES_FINALISE

CONTAINS

  !
  !================================================================================================================================
  !

  !>Finialises the mesh adjacent elements information and deallocates all memory
  SUBROUTINE MESH_ADJACENT_ELEMENT_FINALISE(MESH_ADJACENT_ELEMENT,ERR,ERROR,*)
    
    !Argument variables
    TYPE(MeshAdjacentElementType) :: MESH_ADJACENT_ELEMENT !<The mesh adjacent element to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("MESH_ADJACENT_ELEMENT_FINALISE",ERR,ERROR,*999)

    MESH_ADJACENT_ELEMENT%numberOfAdjacentElements=0
    IF(ALLOCATED(MESH_ADJACENT_ELEMENT%adjacentElements)) DEALLOCATE(MESH_ADJACENT_ELEMENT%adjacentElements)
       
    EXITS("MESH_ADJACENT_ELEMENT_FINALISE")
    RETURN
999 ERRORSEXITS("MESH_ADJACENT_ELEMENT_FINALISE",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MESH_ADJACENT_ELEMENT_FINALISE

  !
  !================================================================================================================================
  !
  !>Initalises the mesh adjacent elements information.
  SUBROUTINE MESH_ADJACENT_ELEMENT_INITIALISE(MESH_ADJACENT_ELEMENT,ERR,ERROR,*)
    
    !Argument variables
    TYPE(MeshAdjacentElementType) :: MESH_ADJACENT_ELEMENT !<The mesh adjacent element to initialise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("MESH_ADJACENT_ELEMENT_INITIALISE",ERR,ERROR,*999)

    MESH_ADJACENT_ELEMENT%numberOfAdjacentElements=0
       
    EXITS("MESH_ADJACENT_ELEMENT_INITIALISE")
    RETURN
999 ERRORSEXITS("MESH_ADJACENT_ELEMENT_INITIALISE",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MESH_ADJACENT_ELEMENT_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a mesh. \see OpenCMISS::Iron::cmfe_MeshCreateFinish
  SUBROUTINE MESH_CREATE_FINISH(MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(MeshType), POINTER :: MESH !<A pointer to the mesh to finish creating
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx
    LOGICAL :: FINISHED
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(ASSOCIATED(MESH%TOPOLOGY)) THEN
        !Check that the mesh component elements have been finished
        FINISHED=.TRUE.
        DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
          IF(ASSOCIATED(MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS)) THEN
            IF(.NOT.MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%elementsFinished) THEN
              LOCAL_ERROR="The elements for mesh component "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
                & " have not been finished"
              FINISHED=.FALSE.
              EXIT
            ENDIF
          ELSE
            LOCAL_ERROR="The elements for mesh topology component "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
              & " are not associated"
            FINISHED=.FALSE.
            EXIT
          ENDIF
        ENDDO !component_idx
        IF(.NOT.FINISHED) CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        MESH%meshFinished=.TRUE.
        !Calulcate the mesh topology
        DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
          CALL MeshTopology_Calculate(MESH%TOPOLOGY(component_idx)%PTR,ERR,ERROR,*999)
        ENDDO !component_idx
      ELSE
        CALL FlagError("Mesh topology is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Mesh user number       = ",MESH%userNumber,ERR,ERROR,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Global number        = ",MESH%globalNumber,ERR,ERROR,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of dimensions = ",MESH%numberOfDimensions,ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_CREATE_FINISH")
    RETURN
999 ERRORSEXITS("MESH_CREATE_FINISH",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MESH_CREATE_FINISH
        
  !
  !================================================================================================================================
  !

  !>Starts the process of creating a mesh 
  SUBROUTINE MESH_CREATE_START_GENERIC(MESHES,USER_NUMBER,numberOfDimensions,MESH,ERR,ERROR,*)
    
    !Argument variables
    TYPE(MeshesType), POINTER :: MESHES !<The pointer to the meshes
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the mesh to create
    INTEGER(INTG), INTENT(IN) :: numberOfDimensions !<The number of dimensions in the mesh.
    TYPE(MeshType), POINTER :: MESH !<On return, a pointer to the mesh. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,mesh_idx
    TYPE(MeshType), POINTER :: NEW_MESH
    TYPE(MeshPtrType), POINTER :: NEW_MESHES(:)
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    NULLIFY(NEW_MESH)
    NULLIFY(NEW_MESHES)

    ENTERS("MESH_CREATE_START_GENERIC",ERR,ERROR,*997)

    IF(ASSOCIATED(MESHES)) THEN
      IF(ASSOCIATED(MESH)) THEN
        CALL FlagError("Mesh is already associated.",ERR,ERROR,*997)
      ELSE
        CALL MESH_INITIALISE(NEW_MESH,ERR,ERROR,*999)
        !Set default mesh values
        NEW_MESH%userNumber=USER_NUMBER
        NEW_MESH%globalNumber=MESHES%numberOfMeshes+1
        NEW_MESH%MESHES=>MESHES
        NEW_MESH%numberOfDimensions=numberOfDimensions
        NEW_MESH%NUMBER_OF_COMPONENTS=1
        NEW_MESH%surroundingElementsCalculate=.true. !default true
        !Initialise mesh topology and decompositions
        CALL MeshTopology_Initialise(NEW_MESH,ERR,ERROR,*999)
        CALL DECOMPOSITIONS_INITIALISE(NEW_MESH,ERR,ERROR,*999)
        !Add new mesh into list of meshes 
        ALLOCATE(NEW_MESHES(MESHES%numberOfMeshes+1),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate new meshes",ERR,ERROR,*999)
        DO mesh_idx=1,MESHES%numberOfMeshes
          NEW_MESHES(mesh_idx)%PTR=>MESHES%MESHES(mesh_idx)%PTR
        ENDDO !mesh_idx
        NEW_MESHES(MESHES%numberOfMeshes+1)%PTR=>NEW_MESH
        IF(ASSOCIATED(MESHES%MESHES)) DEALLOCATE(MESHES%MESHES)
        MESHES%MESHES=>NEW_MESHES
        MESHES%numberOfMeshes=MESHES%numberOfMeshes+1
        MESH=>NEW_MESH
      ENDIF
    ELSE
      CALL FlagError("Meshes is not associated.",ERR,ERROR,*997)
    ENDIF
      
    EXITS("MESH_CREATE_START_GENERIC")
    RETURN
999 CALL MESH_FINALISE(NEW_MESH,DUMMY_ERR,DUMMY_ERROR,*998)
998 IF(ASSOCIATED(NEW_MESHES)) DEALLOCATE(NEW_MESHES)
    NULLIFY(MESH)    
997 ERRORSEXITS("MESH_CREATE_START_GENERIC",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MESH_CREATE_START_GENERIC

  !
  !================================================================================================================================
  !

  !>Starts the process of creating a mesh defined by a user number with the specified numberOfDimensions in an interface. \see OpenCMISS::Iron::cmfe_MeshCreateStart
  !>Default values set for the MESH's attributes are:
  !>- NUMBER_OF_COMPONENTS: 1
  SUBROUTINE MESH_CREATE_START_INTERFACE(USER_NUMBER,INTERFACE,numberOfDimensions,MESH,ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the mesh to create
    TYPE(InterfaceType), POINTER :: INTERFACE !<A pointer to the interface to create the mesh on
    INTEGER(INTG), INTENT(IN) :: numberOfDimensions !<The number of dimensions in the mesh.
    TYPE(MeshType), POINTER :: MESH !<On exit, a pointer to the created mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(RegionType), POINTER :: PARENT_REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_CREATE_START_INTERFACE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(ASSOCIATED(MESH)) THEN
        CALL FlagError("Mesh is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(MESH)
        IF(ASSOCIATED(INTERFACE%MESHES)) THEN
          CALL Mesh_UserNumberFindGeneric(USER_NUMBER,INTERFACE%MESHES,MESH,ERR,ERROR,*999)
          IF(ASSOCIATED(MESH)) THEN
            LOCAL_ERROR="Mesh number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
              & " has already been created on interface number "// &
              & TRIM(NUMBER_TO_VSTRING(INTERFACE%userNumber,"*",ERR,ERROR))//"."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            IF(ASSOCIATED(INTERFACE%INTERFACES)) THEN
              PARENT_REGION=>INTERFACE%INTERFACES%parentRegion
              IF(ASSOCIATED(PARENT_REGION)) THEN
                IF(ASSOCIATED(PARENT_REGION%coordinateSystem)) THEN                  
                  IF(numberOfDimensions>0) THEN
                    IF(numberOfDimensions<=PARENT_REGION%coordinateSystem%numberOfDimensions) THEN
                      CALL MESH_CREATE_START_GENERIC(INTERFACE%MESHES,USER_NUMBER,numberOfDimensions,MESH,ERR,ERROR,*999)
                      MESH%INTERFACE=>INTERFACE
                    ELSE
                      LOCAL_ERROR="Number of mesh dimensions ("//TRIM(NUMBER_TO_VSTRING(numberOfDimensions,"*",ERR,ERROR))// &
                        & ") must be <= number of parent region dimensions ("// &
                        & TRIM(NUMBER_TO_VSTRING(PARENT_REGION%coordinateSystem%numberOfDimensions,"*",ERR,ERROR))//")."
                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Number of mesh dimensions must be > 0.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Parent region coordinate system is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Interfaces parent region is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Interface interfaces is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDIF
        ELSE
          LOCAL_ERROR="The meshes on interface number "//TRIM(NUMBER_TO_VSTRING(INTERFACE%userNumber,"*",ERR,ERROR))// &
            & " are not associated."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Interface is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_CREATE_START_INTERFACE")
    RETURN
999 ERRORSEXITS("MESH_CREATE_START_INTERFACE",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MESH_CREATE_START_INTERFACE

  !
  !================================================================================================================================
  !

  !>Starts the process of creating a mesh defined by a user number with the specified numberOfDimensions in the region identified by REGION. \see OpenCMISS::Iron::cmfe_MeshCreateStart
  !>Default values set for the MESH's attributes are:
  !>- NUMBER_OF_COMPONENTS: 1
  SUBROUTINE MESH_CREATE_START_REGION(USER_NUMBER,REGION,numberOfDimensions,MESH,ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the mesh to create
    TYPE(RegionType), POINTER :: REGION !<A pointer to the region to create the mesh on
    INTEGER(INTG), INTENT(IN) :: numberOfDimensions !<The number of dimensions in the mesh.
    TYPE(MeshType), POINTER :: MESH !<On exit, a pointer to the created mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_CREATE_START_REGION",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(MESH)) THEN
        CALL FlagError("Mesh is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(MESH)
        IF(ASSOCIATED(REGION%MESHES)) THEN
          CALL Mesh_UserNumberFindGeneric(USER_NUMBER,REGION%MESHES,MESH,ERR,ERROR,*999)
          IF(ASSOCIATED(MESH)) THEN
            LOCAL_ERROR="Mesh number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
              & " has already been created on region number "//TRIM(NUMBER_TO_VSTRING(REGION%userNumber,"*",ERR,ERROR))//"."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            IF(ASSOCIATED(REGION%coordinateSystem)) THEN
              IF(numberOfDimensions>0) THEN
                IF(numberOfDimensions<=REGION%coordinateSystem%numberOfDimensions) THEN
                  CALL MESH_CREATE_START_GENERIC(REGION%MESHES,USER_NUMBER,numberOfDimensions,MESH,ERR,ERROR,*999)
                  MESH%REGION=>REGION
                ELSE
                  LOCAL_ERROR="Number of mesh dimensions ("//TRIM(NUMBER_TO_VSTRING(numberOfDimensions,"*",ERR,ERROR))// &
                    & ") must be <= number of region dimensions ("// &
                    & TRIM(NUMBER_TO_VSTRING(REGION%coordinateSystem%numberOfDimensions,"*",ERR,ERROR))//")."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Number of mesh dimensions must be > 0.",ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The coordinate system on region number "//TRIM(NUMBER_TO_VSTRING(REGION%userNumber,"*",ERR,ERROR))// &
                & " are not associated."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ENDIF
        ELSE
          LOCAL_ERROR="The meshes on region number "//TRIM(NUMBER_TO_VSTRING(REGION%userNumber,"*",ERR,ERROR))// &
            & " are not associated."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Region is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("MESH_CREATE_START_REGION")
    RETURN
999 ERRORSEXITS("MESH_CREATE_START_REGION",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MESH_CREATE_START_REGION

  !
  !================================================================================================================================
  !

  !>Destroys the mesh identified by a user number on the given region and deallocates all memory. \see OpenCMISS::Iron::cmfe_MeshDestroy
  SUBROUTINE MESH_DESTROY_NUMBER(USER_NUMBER,REGION,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the mesh to destroy
    TYPE(RegionType), POINTER :: REGION !<A pointer to the region containing the mesh
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: mesh_idx,mesh_position
    LOGICAL :: FOUND
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(MeshType), POINTER :: MESH
    TYPE(MeshPtrType), POINTER :: NEW_MESHES(:)

    NULLIFY(NEW_MESHES)

    ENTERS("MESH_DESTROY_NUMBER",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%MESHES)) THEN

!!TODO: have a mesh_destory_ptr and mesh_destroy_number
        
        !Find the problem identified by the user number
        FOUND=.FALSE.
        mesh_position=0
        DO WHILE(mesh_position<REGION%MESHES%numberOfMeshes.AND..NOT.FOUND)
          mesh_position=mesh_position+1
          IF(REGION%MESHES%MESHES(mesh_position)%PTR%userNumber==USER_NUMBER) FOUND=.TRUE.
        ENDDO
        
        IF(FOUND) THEN
          
          MESH=>REGION%MESHES%MESHES(mesh_position)%PTR

          CALL MESH_FINALISE(MESH,ERR,ERROR,*999)

          !Remove the mesh from the list of meshes
          IF(REGION%MESHES%numberOfMeshes>1) THEN
            ALLOCATE(NEW_MESHES(REGION%MESHES%numberOfMeshes-1),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate new meshes",ERR,ERROR,*999)
            DO mesh_idx=1,REGION%MESHES%numberOfMeshes
              IF(mesh_idx<mesh_position) THEN
                NEW_MESHES(mesh_idx)%PTR=>REGION%MESHES%MESHES(mesh_idx)%PTR
              ELSE IF(mesh_idx>mesh_position) THEN
                REGION%MESHES%MESHES(mesh_idx)%PTR%globalNumber=REGION%MESHES%MESHES(mesh_idx)%PTR%globalNumber-1
                NEW_MESHES(mesh_idx-1)%PTR=>REGION%MESHES%MESHES(mesh_idx)%PTR
              ENDIF
            ENDDO !mesh_idx
            DEALLOCATE(REGION%MESHES%MESHES)
            REGION%MESHES%MESHES=>NEW_MESHES
            REGION%MESHES%numberOfMeshes=REGION%MESHES%numberOfMeshes-1
          ELSE
            DEALLOCATE(REGION%MESHES%MESHES)
            REGION%MESHES%numberOfMeshes=0
          ENDIF
          
        ELSE
          LOCAL_ERROR="Mesh number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
            & " has not been created on region number "//TRIM(NUMBER_TO_VSTRING(REGION%userNumber,"*",ERR,ERROR))
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="The meshes on region number "//TRIM(NUMBER_TO_VSTRING(REGION%userNumber,"*",ERR,ERROR))// &
          & " are not associated"
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Region is not associated",ERR,ERROR,*999)
    ENDIF

    EXITS("MESH_DESTROY_NUMBER")
    RETURN
999 IF(ASSOCIATED(NEW_MESHES)) DEALLOCATE(NEW_MESHES)
    ERRORSEXITS("MESH_DESTROY_NUMBER",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE MESH_DESTROY_NUMBER

  !
  !================================================================================================================================
  !

  !>Destroys the mesh and deallocates all memory. \see OpenCMISS::Iron::cmfe_MeshDestroy
  SUBROUTINE MESH_DESTROY(MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(MeshType), POINTER :: MESH !<A pointer to the mesh to destroy.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: mesh_idx,mesh_position
    TYPE(MeshesType), POINTER :: MESHES
    TYPE(MeshPtrType), POINTER :: NEW_MESHES(:)

    NULLIFY(NEW_MESHES)

    ENTERS("MESH_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      MESHES=>MESH%MESHES
      IF(ASSOCIATED(MESHES)) THEN
        mesh_position=MESH%globalNumber
          
        CALL MESH_FINALISE(MESH,ERR,ERROR,*999)

        !Remove the mesh from the list of meshes
        IF(MESHES%numberOfMeshes>1) THEN
          ALLOCATE(NEW_MESHES(MESHES%numberOfMeshes-1),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate new meshes.",ERR,ERROR,*999)
          DO mesh_idx=1,MESHES%numberOfMeshes
            IF(mesh_idx<mesh_position) THEN
              NEW_MESHES(mesh_idx)%PTR=>MESHES%MESHES(mesh_idx)%PTR
            ELSE IF(mesh_idx>mesh_position) THEN
              MESHES%MESHES(mesh_idx)%PTR%globalNumber=MESHES%MESHES(mesh_idx)%PTR%globalNumber-1
              NEW_MESHES(mesh_idx-1)%PTR=>MESHES%MESHES(mesh_idx)%PTR
            ENDIF
          ENDDO !mesh_idx
          DEALLOCATE(MESHES%MESHES)
          MESHES%MESHES=>NEW_MESHES
          MESHES%numberOfMeshes=MESHES%numberOfMeshes-1
        ELSE
          DEALLOCATE(MESHES%MESHES)
          MESHES%numberOfMeshes=0
        ENDIF
      ELSE
        CALL FlagError("The mesh meshes is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("MESH_DESTROY")
    RETURN
999 IF(ASSOCIATED(NEW_MESHES)) DEALLOCATE(NEW_MESHES)
    ERRORSEXITS("MESH_DESTROY",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE MESH_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalises a mesh and deallocates all memory.
  SUBROUTINE MESH_FINALISE(MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(MeshType), POINTER :: MESH !<A pointer to the mesh to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("MESH_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      CALL MeshTopology_Finalise(MESH,ERR,ERROR,*999)
      CALL DECOMPOSITIONS_FINALISE(MESH,ERR,ERROR,*999)
!      IF(ASSOCIATED(MESH%INTF)) CALL INTERFACE_MESH_FINALISE(MESH,ERR,ERROR,*999)  ! <<??>>
      DEALLOCATE(MESH)
    ENDIF
 
    EXITS("MESH_FINALISE")
    RETURN
999 ERRORSEXITS("MESH_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE MESH_FINALISE

  !
  !================================================================================================================================
  !

  !>Returns a region nodes pointer corresponding to the mesh global nodes accounting for interfaces.
  SUBROUTINE MeshGlobalNodesGet(mesh,nodes,err,error,*)

    !Argument variables
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh to get the global nodes for
    TYPE(NodesType), POINTER :: nodes !<On return, the nodes pointer corresponding to the global nodes for the mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(InterfaceType), POINTER :: interface
    TYPE(RegionType), POINTER :: region
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("MeshGlobalNodesGet",err,error,*999)
    
    IF(ASSOCIATED(mesh)) THEN
      IF(ASSOCIATED(nodes)) THEN
        CALL FlagError("Nodes is already associated.",err,error,*999)
      ELSE
        NULLIFY(nodes)
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
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated.",err,error,*999)
    ENDIF
        
    EXITS("MeshGlobalNodesGet")
    RETURN
999 ERRORSEXITS("MeshGlobalNodesGet",err,error)
    RETURN 1
    
  END SUBROUTINE MeshGlobalNodesGet
  
  !
  !================================================================================================================================
  !

  !>Initialises a mesh.
  SUBROUTINE MESH_INITIALISE(MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(MeshType), POINTER :: MESH !<A pointer to the mesh to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("MESH_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      CALL FlagError("Mesh is already associated.",ERR,ERROR,*999)
    ELSE
      ALLOCATE(MESH,STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate new mesh.",ERR,ERROR,*999)
      MESH%userNumber=0
      MESH%globalNumber=0
      MESH%meshFinished=.FALSE.
      NULLIFY(MESH%MESHES)
      NULLIFY(MESH%REGION)
      NULLIFY(MESH%INTERFACE)
      NULLIFY(MESH%generatedMesh)
      MESH%numberOfDimensions=0
      MESH%NUMBER_OF_COMPONENTS=0
      MESH%MESH_EMBEDDED=.FALSE.
      NULLIFY(MESH%EMBEDDING_MESH)
      MESH%NUMBER_OF_EMBEDDED_MESHES=0
      NULLIFY(MESH%EMBEDDED_MESHES)
      MESH%numberOfElements=0
      NULLIFY(MESH%TOPOLOGY)
      NULLIFY(MESH%DECOMPOSITIONS)
    ENDIF
    
    EXITS("MESH_INITIALISE")
    RETURN
999 ERRORSEXITS("MESH_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESH_INITIALISE

  !
  !================================================================================================================================
  !
  
  !>Gets the number of mesh components for a mesh identified by a pointer. \see OpenCMISS::Iron::cmfe_MeshNumberOfComponentsGet
  SUBROUTINE MESH_NUMBER_OF_COMPONENTS_GET(MESH,NUMBER_OF_COMPONENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MeshType), POINTER :: MESH !<A pointer to the mesh to get the number of components for
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_COMPONENTS !<On return, the number of components in the specified mesh.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("MESH_NUMBER_OF_COMPONENTS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(MESH%meshFinished) THEN
        NUMBER_OF_COMPONENTS=MESH%NUMBER_OF_COMPONENTS
      ELSE
        CALL FlagError("Mesh has not finished",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_NUMBER_OF_COMPONENTS_GET")
    RETURN
999 ERRORSEXITS("MESH_NUMBER_OF_COMPONENTS_GET",ERR,ERROR)    
    RETURN
  END SUBROUTINE MESH_NUMBER_OF_COMPONENTS_GET

  !
  !================================================================================================================================
  !

  !>Changes/sets the number of mesh components for a mesh. \see OpenCMISS::Iron::cmfe_MeshNumberOfComponentsSet
  SUBROUTINE MESH_NUMBER_OF_COMPONENTS_SET(MESH,NUMBER_OF_COMPONENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MeshType), POINTER :: MESH !<A pointer to the mesh to set the number of components for
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_COMPONENTS !<The number of components to set.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(MeshComponentTopologyPtrType), POINTER :: NEW_TOPOLOGY(:)

    NULLIFY(NEW_TOPOLOGY)
    
    ENTERS("MESH_NUMBER_OF_COMPONENTS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(NUMBER_OF_COMPONENTS>0) THEN
        IF(MESH%meshFinished) THEN
          CALL FlagError("Mesh has been finished",ERR,ERROR,*999)
        ELSE
          IF(NUMBER_OF_COMPONENTS/=MESH%NUMBER_OF_COMPONENTS) THEN
            ALLOCATE(NEW_TOPOLOGY(NUMBER_OF_COMPONENTS),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate new topology",ERR,ERROR,*999)
            IF(NUMBER_OF_COMPONENTS<MESH%NUMBER_OF_COMPONENTS) THEN
              DO component_idx=1,NUMBER_OF_COMPONENTS
                NEW_TOPOLOGY(component_idx)%PTR=>MESH%TOPOLOGY(component_idx)%PTR
              ENDDO !component_idx
            ELSE !NUMBER_OF_COMPONENTS>MESH%NUMBER_OF_COMPONENTS
              DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
                NEW_TOPOLOGY(component_idx)%PTR=>MESH%TOPOLOGY(component_idx)%PTR
              ENDDO !component_idx
!!TODO \todo sort out mesh_topology initialise/finalise so that they allocate and deal with this below then call that routine
              DO component_idx=MESH%NUMBER_OF_COMPONENTS+1,NUMBER_OF_COMPONENTS
                ALLOCATE(NEW_TOPOLOGY(component_idx)%PTR,STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate new topology component",ERR,ERROR,*999)
                NEW_TOPOLOGY(component_idx)%PTR%mesh=>mesh
                NEW_TOPOLOGY(component_idx)%PTR%meshComponentNumber=component_idx
                NULLIFY(NEW_TOPOLOGY(component_idx)%PTR%elements)
                NULLIFY(NEW_TOPOLOGY(component_idx)%PTR%nodes)
                NULLIFY(NEW_TOPOLOGY(component_idx)%PTR%dofs)
                NULLIFY(NEW_TOPOLOGY(component_idx)%PTR%dataPoints)
                !Initialise the topology components
                CALL MESH_TOPOLOGY_ELEMENTS_INITIALISE(NEW_TOPOLOGY(component_idx)%PTR,ERR,ERROR,*999)
                CALL MeshTopology_NodesInitialise(NEW_TOPOLOGY(component_idx)%PTR,ERR,ERROR,*999)
                CALL MeshTopology_DofsInitialise(NEW_TOPOLOGY(component_idx)%PTR,ERR,ERROR,*999)
                CALL MESH_TOPOLOGY_DATA_POINTS_INITIALISE(NEW_TOPOLOGY(component_idx)%PTR,ERR,ERROR,*999)
              ENDDO !component_idx
            ENDIF
            IF(ASSOCIATED(MESH%TOPOLOGY)) DEALLOCATE(MESH%TOPOLOGY)
            MESH%TOPOLOGY=>NEW_TOPOLOGY
            MESH%NUMBER_OF_COMPONENTS=NUMBER_OF_COMPONENTS
          ENDIF
        ENDIF
      ELSE
        LOCAL_ERROR="The specified number of mesh components ("//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
          & ") is illegal. You must have >0 mesh components"
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF      
    ELSE
      CALL FlagError("Mesh is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_NUMBER_OF_COMPONENTS_SET")
    RETURN
!!TODO: tidy up memory deallocation on error
999 ERRORSEXITS("MESH_NUMBER_OF_COMPONENTS_SET",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MESH_NUMBER_OF_COMPONENTS_SET

  !
  !================================================================================================================================
  !
  
  !>Gets the number of elements for a mesh identified by a pointer. \see OpenCMISS::Iron::cmfe_MeshNumberOfElementsGet
  SUBROUTINE MESH_NUMBER_OF_ELEMENTS_GET(MESH,numberOfElements,ERR,ERROR,*)

    !Argument variables
    TYPE(MeshType), POINTER :: MESH !<A pointer to the mesh to get the number of elements for
    INTEGER(INTG), INTENT(OUT) :: numberOfElements !<On return, the number of elements in the specified mesh
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("MESH_NUMBER_OF_ELEMENTS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(MESH%meshFinished) THEN
        numberOfElements=MESH%numberOfElements
      ELSE
        CALL FlagError("Mesh has not been finished",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_NUMBER_OF_ELEMENTS_GET")
    RETURN
999 ERRORSEXITS("MESH_NUMBER_OF_ELEMENTS_GET",ERR,ERROR)    
    RETURN
  END SUBROUTINE MESH_NUMBER_OF_ELEMENTS_GET

  !
  !================================================================================================================================
  !

  !>Changes/sets the number of elements for a mesh. \see OpenCMISS::Iron::cmfe_MeshNumberOfElementsSet
  SUBROUTINE MESH_NUMBER_OF_ELEMENTS_SET(MESH,numberOfElements,ERR,ERROR,*)

    !Argument variables
    TYPE(MeshType), POINTER :: MESH !<A pointer to the mesh to set the number of elements for
    INTEGER(INTG), INTENT(IN) :: numberOfElements !<The number of elements to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_NUMBER_OF_ELEMENTS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(numberOfElements>0) THEN
        IF(MESH%meshFinished) THEN
          CALL FlagError("Mesh has been finished.",ERR,ERROR,*999)
        ELSE
          IF(numberOfElements/=MESH%numberOfElements) THEN
            IF(ASSOCIATED(MESH%TOPOLOGY)) THEN
              DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
                IF(ASSOCIATED(MESH%TOPOLOGY(component_idx)%PTR)) THEN
                  IF(ASSOCIATED(MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS)) THEN
                    IF(MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%numberOfElements>0) THEN
!!TODO: Reallocate the elements and copy information. 
                      CALL FlagError("Not implemented.",ERR,ERROR,*999)
                    ENDIF
                  ENDIF
                ELSE
                  CALL FlagError("Mesh topology component pointer is not associated.",ERR,ERROR,*999)
                ENDIF
              ENDDO !component_idx
            ELSE
              CALL FlagError("Mesh topology is not associated.",ERR,ERROR,*999)
            ENDIF
            MESH%numberOfElements=numberOfElements
          ENDIF
        ENDIF
      ELSE
        LOCAL_ERROR="The specified number of elements ("//TRIM(NUMBER_TO_VSTRING(numberOfElements,"*",ERR,ERROR))// &
          & ") is invalid. You must have > 0 elements."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF      
    ELSE
      CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_NUMBER_OF_ELEMENTS_SET")
    RETURN
999 ERRORSEXITS("MESH_NUMBER_OF_ELEMENTS_SET",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MESH_NUMBER_OF_ELEMENTS_SET

  !
  !================================================================================================================================
  !

  !>Changes/sets the surrounding elements calculate flag. \see OpenCMISS::Iron::cmfe_MeshSurroundingElementsCalculateSet
  SUBROUTINE MESH_SURROUNDING_ELEMENTS_CALCULATE_SET(MESH,SURROUNDING_ELEMENTS_CALCULATE_FLAG,ERR,ERROR,*)

    !Argument variables
    TYPE(MeshType), POINTER :: MESH !<A pointer to the mesh to set the surrounding elements calculate flag for
    LOGICAL, INTENT(IN) :: SURROUNDING_ELEMENTS_CALCULATE_FLAG !<The surrounding elements calculate flag
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    ENTERS("MESH_SURROUNDING_ELEMENTS_CALCULATE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(MESH%meshFinished) THEN
        CALL FlagError("Mesh has been finished.",ERR,ERROR,*999)
      ELSE
        MESH%surroundingElementsCalculate=SURROUNDING_ELEMENTS_CALCULATE_FLAG
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_SURROUNDING_ELEMENTS_CALCULATE_SET")
    RETURN
999 ERRORSEXITS("MESH_SURROUNDING_ELEMENTS_CALCULATE_SET",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MESH_SURROUNDING_ELEMENTS_CALCULATE_SET

  !
  !================================================================================================================================
  !

  !>Calculates the mesh topology.
  SUBROUTINE MeshTopology_Calculate(topology,err,error,*)    

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: topology !<A pointer to the mesh topology to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshTopology_Calculate",err,error,*999)

    IF(ASSOCIATED(topology)) THEN
      !Calculate the nodes used in the mesh
      CALL MeshTopology_NodesCalculate(topology,err,error,*999)
      !Calculate the elements surrounding the nodes in a mesh
      CALL MeshTopology_SurroundingElementsCalculate(topology,err,error,*999)
      !Calculate the number of derivatives at each node in a mesh
      CALL MeshTopology_NodesDerivativesCalculate(topology,err,error,*999)
      !Calculate the number of versions for each derivative at each node in a mesh
      CALL MeshTopology_NodesVersionCalculate(topology,err,error,*999)
      !Calculate the elements surrounding the elements in the mesh
      CALL MeshTopology_ElementsAdjacentElementsCalculate(topology,err,error,*999)
      !Calculate the boundary nodes and elements in the mesh
      CALL MeshTopology_BoundaryCalculate(topology,err,error,*999)
      !Calculate the elements surrounding the elements in the mesh
      CALL MeshTopology_DofsCalculate(topology,err,error,*999)
    ELSE
      CALL FlagError("Topology is not associated",err,error,*999)
    ENDIF
    
    EXITS("MeshTopology_Calculate")
    RETURN
999 ERRORSEXITS("MeshTopology_Calculate",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_Calculate
  
  !
  !===============================================================================================================================
  !

  !>Calculates the boundary nodes and elements for a mesh topology. 
  SUBROUTINE MeshTopology_BoundaryCalculate(topology,err,error,*)

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: topology !<A pointer to the mesh topology to calculate the boundary for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementIdx,localNodeIdx,matchIndex,nodeIdx,xiCoordIdx,xiDirection
    TYPE(BasisType), POINTER :: basis
    TYPE(MeshElementsType), POINTER :: elements
    TYPE(MeshNodesType), POINTER :: nodes
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("MeshTopology_BoundaryCalculate",err,error,*999)

    IF(ASSOCIATED(topology)) THEN
      nodes=>topology%nodes
      IF(ASSOCIATED(nodes)) THEN
        elements=>topology%elements
        IF(ASSOCIATED(elements)) THEN
          DO elementIdx=1,elements%numberOfElements
            basis=>elements%elements(elementIdx)%basis
            SELECT CASE(basis%type)
            CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
              DO xiCoordIdx=-basis%numberOfXiCoordinates,basis%numberOfXiCoordinates
                IF(xiCoordIdx/=0) THEN
                  IF(elements%elements(elementIdx)%adjacentElements(xiCoordIdx)%numberOfAdjacentElements==0) THEN
                    elements%elements(elementIdx)%boundaryElement=.TRUE.
                    IF(xiCoordIdx<0) THEN
                      xiDirection=-xiCoordIdx
                      matchIndex=1
                    ELSE
                      xiDirection=xiCoordIdx
                      matchIndex=BASIS%numberOfNodesXiC(xiCoordIdx)
                    ENDIF                    
                    DO localNodeIdx=1,BASIS%numberOfNodes
                      IF(basis%nodePositionIndex(localNodeIdx,XIDIRECTION)==matchIndex) THEN
                        nodeIdx=elements%elements(elementIdx)%meshElementNodes(localNodeIdx)
                        nodes%nodes(nodeIdx)%boundaryNode=.TRUE.
                      ENDIF
                    ENDDO !nn
                  ENDIF
                ENDIF
              ENDDO !xiCoordIdx            
            CASE(BASIS_SIMPLEX_TYPE)
              elements%elements(elementIdx)%boundaryElement=.FALSE.
              DO xiCoordIdx=1,basis%numberOfXiCoordinates
                elements%elements(elementIdx)%boundaryElement=elements%elements(elementIdx)%boundaryElement.OR. &
                  & elements%elements(elementIdx)%adjacentElements(xiCoordIdx)%numberOfAdjacentElements==0
                IF(elements%elements(elementIdx)%adjacentElements(xiCoordIdx)%numberOfAdjacentElements==0) THEN
                  DO localNodeIdx=1,basis%numberOfNodes
                    IF(basis%nodePositionIndex(localNodeIdx,xiCoordIdx)==1) THEN
                      nodeIdx=elements%elements(elementIdx)%meshElementNodes(localNodeIdx)
                      nodes%nodes(nodeIdx)%boundaryNode=.TRUE.
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
              localError="The basis type of "//TRIM(NumberToVString(basis%TYPE,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ENDDO !elementIdx
        ELSE
          CALL FlagError("Topology elements is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Topology nodes is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",err,error,*999)
    ENDIF

    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Boundary elements:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Number of elements = ",elements%numberOfElements,err,error,*999)
      DO elementIdx=1,elements%numberOfElements
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Element : ",elementIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Boundary element = ",elements%elements(elementIdx)%boundaryElement, &
          & err,error,*999)        
      ENDDO !elementIdx
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Boundary nodes:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Number of nodes = ",nodes%numberOfNodes,err,error,*999)
      DO nodeIdx=1,nodes%numberOfNodes
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Node : ",nodeIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Boundary node = ",nodes%nodes(nodeIdx)%boundaryNode,err,error,*999)        
      ENDDO !elementIdx            
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
  SUBROUTINE MeshTopology_DofsCalculate(topology,err,error,*)

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: topology !<A pointer to the mesh topology to calculate the dofs for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: derivativeIdx,nodeIdx,numberOfDofs,versionIdx
    TYPE(MeshDofsType), POINTER :: dofs
    TYPE(MeshNodesType), POINTER :: nodes

    ENTERS("MeshTopology_DofsCalculate",err,error,*999)

    IF(ASSOCIATED(topology)) THEN
      nodes=>topology%nodes
      IF(ASSOCIATED(nodes)) THEN
        dofs=>topology%dofs
        IF(ASSOCIATED(dofs)) THEN
          numberOfDofs=0
          DO nodeIdx=1,nodes%numberOfNodes
            DO derivativeIdx=1,nodes%nodes(nodeIdx)%numberOfDerivatives
              ALLOCATE(nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%dofIndex( &
                & nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%numberOfVersions),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate mesh topology node derivative version dof index.",err,error,*999)
              DO versionIdx=1,nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%numberOfVersions
                numberOfDofs=numberOfDofs+1
                nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%dofIndex(versionIdx)=numberOfDofs
              ENDDO !versionIdx
            ENDDO !derivativeIdx
          ENDDO !nodeIdx
          dofs%numberOfDofs=numberOfDofs
        ELSE
          CALL FlagError("Topology dofs is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Topology nodes is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",err,error,*999)
    ENDIF
 
    EXITS("MeshTopology_DofsCalculate")
    RETURN
999 ERRORSEXITS("MeshTopology_DofsCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_DofsCalculate

  !
  !===============================================================================================================================
  !

  !>Finalises the dof data structures for a mesh topology and deallocates any memory. \todo pass in dofs
  SUBROUTINE MeshTopology_DofsFinalise(dofs,err,error,*)

    !Argument variables
    TYPE(MeshDofsType), POINTER :: dofs !<A pointer to the mesh topology to finalise the dofs for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshTopology_DofsFinalise",err,error,*999)

    IF(ASSOCIATED(dofs)) THEN
      DEALLOCATE(dofs)
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
  SUBROUTINE MeshTopology_DofsInitialise(topology,err,error,*)

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: topology !<A pointer to the mesh topology to initialise the dofs for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("MeshTopology_DofsInitialise",err,error,*998)

    IF(ASSOCIATED(topology)) THEN
      IF(ASSOCIATED(topology%dofs)) THEN
        CALL FlagError("Mesh already has topology dofs associated",err,error,*998)
      ELSE
        ALLOCATE(topology%dofs,STAT=err)
        IF(ERR/=0) CALL FlagError("Could not allocate topology dofs",err,error,*999)
        topology%dofs%numberOfDofs=0
        topology%dofs%meshComponentTopology=>topology
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",err,error,*998)
    ENDIF
    
    EXITS("MeshTopology_DofsInitialise")
    RETURN
999 CALL MeshTopology_DofsFinalise(topology%dofs,dummyErr,dummyError,*998)
998 ERRORSEXITS("MeshTopology_DofsInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_DofsInitialise

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating elements for a specified mesh component in a mesh topology. \see OpenCMISS::Iron::cmfe_MeshElementsCreateFinish
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH(ELEMENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the mesh elements to finish creating
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ne
    TYPE(MeshType), POINTER :: mesh
    TYPE(MeshComponentTopologyType), POINTER :: meshComponentTopology

    ENTERS("MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(ELEMENTS%elementsFinished) THEN
        CALL FlagError("Mesh elements have already been finished.",ERR,ERROR,*999)
      ELSE        
        ELEMENTS%elementsFinished=.TRUE.
      ENDIF
    ELSE
      CALL FlagError("Mesh elements is not associated.",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      meshComponentTopology=>elements%meshComponentTopology
      IF(ASSOCIATED(meshComponentTopology)) THEN
        mesh=>meshComponentTopology%mesh
        IF(ASSOCIATED(MESH)) THEN
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Number of global elements = ",MESH%numberOfElements, &
            & ERR,ERROR,*999)
          DO ne=1,MESH%numberOfElements
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Element = ",ne,ERR,ERROR,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Global number        = ",ELEMENTS%ELEMENTS(ne)%globalNumber, &
              & ERR,ERROR,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    User number          = ",ELEMENTS%ELEMENTS(ne)%userNumber, &
              & ERR,ERROR,*999)
            IF(ASSOCIATED(ELEMENTS%ELEMENTS(ne)%BASIS)) THEN
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Basis number         = ",ELEMENTS%ELEMENTS(ne)%BASIS% &
                & userNumber,ERR,ERROR,*999)
            ELSE
              CALL FlagError("Basis is not associated.",ERR,ERROR,*999)
            ENDIF
            IF(ALLOCATED(ELEMENTS%ELEMENTS(ne)%userElementNodes)) THEN
              CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS%ELEMENTS(ne)% BASIS%numberOfNodes,8,8, &
                & ELEMENTS%ELEMENTS(ne)%userElementNodes,'("    User element nodes   =",8(X,I6))','(26X,8(X,I6))', &
                & ERR,ERROR,*999)
            ELSE
              CALL FlagError("User element nodes are not associated.",ERR,ERROR,*999)
            ENDIF
            IF(ALLOCATED(ELEMENTS%ELEMENTS(ne)%globalElementNodes)) THEN
              CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS%ELEMENTS(ne)%BASIS%numberOfNodes,8,8, &
                & ELEMENTS%ELEMENTS(ne)%globalElementNodes,'("    Global element nodes =",8(X,I6))','(26X,8(X,I6))', &
                & ERR,ERROR,*999)
            ELSE
              CALL FlagError("Global element nodes are not associated.",ERR,ERROR,*999)
            ENDIF
          ENDDO !ne
        ELSE
          CALL FlagError("Mesh component topology mesh is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Mesh elements mesh component topology is not associated.",ERR,ERROR,*999)
      ENDIF
    ENDIF
    
    EXITS("MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH
    
  !
  !================================================================================================================================
  !

  !>Starts the process of creating elements in the mesh component identified by MESH and component_idx. The elements will be created with a default basis of BASIS. ELEMENTS is the returned pointer to the MESH_ELEMENTS data structure. \see OpenCMISS::Iron::cmfe_MeshElementsCreateStart
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_CREATE_START(MESH,MESH_COMPONENT_NUMBER,BASIS,ELEMENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MeshType), POINTER :: MESH !<A pointer to the mesh to start creating the elements on
    INTEGER(INTG), INTENT(IN) :: MESH_COMPONENT_NUMBER !<The mesh component number
    TYPE(BasisType), POINTER :: BASIS !<A pointer to the default basis to use
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<On return, a pointer to the created mesh elements
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,INSERT_STATUS,ne
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR
 
    ENTERS("MESH_TOPOLOGY_ELEMENTS_CREATE_START",ERR,ERROR,*999)
    
    IF(ASSOCIATED(MESH)) THEN     
      IF(MESH_COMPONENT_NUMBER>0.AND.MESH_COMPONENT_NUMBER<=MESH%NUMBER_OF_COMPONENTS) THEN
        IF(ASSOCIATED(ELEMENTS)) THEN
          CALL FlagError("Elements is already associated.",ERR,ERROR,*999)
        ELSE
          IF(ASSOCIATED(MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR)) THEN
            IF(ASSOCIATED(MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR%ELEMENTS)) THEN
              ELEMENTS=>MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR%ELEMENTS
              IF(ASSOCIATED(ELEMENTS%ELEMENTS)) THEN
                CALL FlagError("Mesh topology already has elements associated",ERR,ERROR,*998)
              ELSE
                IF(ASSOCIATED(BASIS)) THEN
                  MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR%meshComponentNumber=MESH_COMPONENT_NUMBER
                  ALLOCATE(ELEMENTS%ELEMENTS(MESH%numberOfElements),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate individual elements",ERR,ERROR,*999)
                  ELEMENTS%numberOfElements=MESH%numberOfElements !Psuedo inheritance of the number of elements
                  CALL TREE_CREATE_START(ELEMENTS%elementsTree,ERR,ERROR,*999)
                  CALL TREE_INSERT_TYPE_SET(ELEMENTS%elementsTree,TREE_NO_DUPLICATES_ALLOWED,ERR,ERROR,*999)
                  CALL TREE_CREATE_FINISH(ELEMENTS%elementsTree,ERR,ERROR,*999)
                  ELEMENTS%elementsFinished=.FALSE.
                  !Set up the default values and allocate element structures
                  DO ne=1,ELEMENTS%numberOfElements
                    CALL MESH_TOPOLOGY_ELEMENT_INITIALISE(ELEMENTS%ELEMENTS(ne),ERR,ERROR,*999)
                    ELEMENTS%ELEMENTS(ne)%globalNumber=ne
                    ELEMENTS%ELEMENTS(ne)%userNumber=ne
                    CALL TREE_ITEM_INSERT(ELEMENTS%elementsTree,ne,ne,INSERT_STATUS,ERR,ERROR,*999)
                    ELEMENTS%ELEMENTS(ne)%BASIS=>BASIS
                    ALLOCATE(ELEMENTS%ELEMENTS(ne)%userElementNodes(BASIS%numberOfNodes),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate user element nodes",ERR,ERROR,*999)
                    ALLOCATE(ELEMENTS%ELEMENTS(ne)%globalElementNodes(BASIS%numberOfNodes),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate global element nodes",ERR,ERROR,*999)
                    ELEMENTS%ELEMENTS(ne)%userElementNodes=1
                    ELEMENTS%ELEMENTS(ne)%globalElementNodes=1
                    ALLOCATE(ELEMENTS%ELEMENTS(ne)%userElementNodeVersions(BASIS%maximumNumberOfDerivatives, &
                      & BASIS%numberOfNodes),STAT=ERR)
                    IF(ERR/=0) CALL FlagError("Could not allocate global element nodes versions",ERR,ERROR,*999)
                    ELEMENTS%ELEMENTS(ne)%userElementNodeVersions = 1
                  ENDDO !ne
                ELSE
                  CALL FlagError("Basis is not associated",ERR,ERROR,*999)
                ENDIF
              ENDIF
            ELSE
              CALL FlagError("Mesh topology elements is not associated",ERR,ERROR,*998)
            ENDIF
          ELSE
            CALL FlagError("Mesh topology is not associated",ERR,ERROR,*998)
          ENDIF
        ENDIF
      ELSE
        LOCAL_ERROR="The specified mesh component number of "//TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT_NUMBER,"*",ERR,ERROR))// &
          & " is invalid. The component number must be between 1 and "// &
          & TRIM(NUMBER_TO_VSTRING(MESH%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated",ERR,ERROR,*998)
    ENDIF

    EXITS("MESH_TOPOLOGY_ELEMENTS_CREATE_START")
    RETURN
999 CALL MESH_TOPOLOGY_ELEMENTS_FINALISE(ELEMENTS,DUMMY_ERR,DUMMY_ERROR,*998)
998 NULLIFY(ELEMENTS)
    ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_CREATE_START",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroys the elements in a mesh topology. \todo as this is a user routine it should take a mesh pointer like create start and finish? Split this into destroy and finalise?
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_DESTROY(ELEMENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the mesh elements to destroy 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("MESH_TOPOLOGY_ELEMENTS_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      CALL MESH_TOPOLOGY_ELEMENTS_FINALISE(ELEMENTS,ERR,ERROR,*999)
    ELSE
      CALL FlagError("Mesh topology is not associated",ERR,ERROR,*999)
    ENDIF

    EXITS("MESH_TOPOLOGY_ELEMENTS_DESTROY")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_DESTROY",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_DESTROY
  
  !
  !================================================================================================================================
  !

  !>Finalises the given mesh topology element.
  SUBROUTINE MESH_TOPOLOGY_ELEMENT_FINALISE(ELEMENT,ERR,ERROR,*)
    
    !Argument variables
    TYPE(MeshElementType) :: ELEMENT !<The mesh element to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nic
    
    ENTERS("MESH_TOPOLOGY_ELEMENT_FINALISE",ERR,ERROR,*999)
    
    IF(ALLOCATED(ELEMENT%userElementNodeVersions)) DEALLOCATE(ELEMENT%userElementNodeVersions)
    IF(ALLOCATED(ELEMENT%userElementNodes)) DEALLOCATE(ELEMENT%userElementNodes)
    IF(ALLOCATED(ELEMENT%globalElementNodes)) DEALLOCATE(ELEMENT%globalElementNodes)
    IF(ALLOCATED(ELEMENT%meshElementNodes)) DEALLOCATE(ELEMENT%meshElementNodes)
    IF(ALLOCATED(ELEMENT%adjacentElements)) THEN
      DO nic=LBOUND(ELEMENT%adjacentElements,1),UBOUND(ELEMENT%adjacentElements,1)
        CALL MESH_ADJACENT_ELEMENT_FINALISE(ELEMENT%adjacentElements(nic),ERR,ERROR,*999)
      ENDDO !nic
      DEALLOCATE(ELEMENT%adjacentElements)
    ENDIF
    
    EXITS("MESH_TOPOLOGY_ELEMENT_FINALISE")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENT_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESH_TOPOLOGY_ELEMENT_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the given mesh topology element.
  SUBROUTINE MESH_TOPOLOGY_ELEMENT_INITIALISE(ELEMENT,ERR,ERROR,*)

    !Argument variables
    TYPE(MeshElementType) :: ELEMENT !<The mesh element to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("MESH_TOPOLOGY_ELEMENT_INITIALISE",ERR,ERROR,*999)

    ELEMENT%userNumber=0
    ELEMENT%globalNumber=0
    NULLIFY(ELEMENT%BASIS)
    ELEMENT%boundaryElement=.FALSE.
    
    EXITS("MESH_TOPOLOGY_ELEMENT_INITIALISE")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENT_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESH_TOPOLOGY_ELEMENT_INITIALISE
  
  !
  !================================================================================================================================
  !

!!MERGE: Take user number
  
  !>Gets the basis for a mesh element identified by a given global number. \todo should take user number \see OpenCMISS::Iron::cmfe_MeshElementsBasisGet
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_GET(globalNumber,ELEMENTS,BASIS,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: globalNumber !<The global number of the element to get the basis for
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the elements to get the basis for \todo before number?
    TYPE(BasisType), POINTER :: BASIS !<On return, a pointer to the basis to get
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(MeshElementType), POINTER :: ELEMENT
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(.NOT.ELEMENTS%elementsFinished) THEN
        CALL FlagError("Elements have been finished",ERR,ERROR,*999)
      ELSE
        IF(globalNumber>=1.AND.globalNumber<=ELEMENTS%numberOfElements) THEN
          ELEMENT=>ELEMENTS%ELEMENTS(globalNumber)
          BASIS=>ELEMENT%BASIS
        ELSE
          LOCAL_ERROR="Global element number "//TRIM(NUMBER_TO_VSTRING(globalNumber,"*",ERR,ERROR))// &
            & " is invalid. The limits are 1 to "//TRIM(NUMBER_TO_VSTRING(ELEMENTS%numberOfElements,"*",ERR,ERROR))
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Elements is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_GET")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_GET",ERR,ERROR)    
    RETURN 1
    
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_GET

  !
  !================================================================================================================================
  !

  !>Changes/sets the basis for a mesh element identified by a given global number. \todo should take user number
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET(globalNumber,ELEMENTS,BASIS,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: globalNumber !<The global number of the element to set the basis for
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the elements to set the basis for \todo before number?
    TYPE(BasisType), POINTER :: BASIS !<A pointer to the basis to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG), ALLOCATABLE :: NEW_userElementNodes(:),NEW_GLOBAL_ELEMENT_NODES(:),NEW_userElementNodeVersions(:,:)
    INTEGER(INTG) :: OVERLAPPING_NUMBER_NODES,OVERLAPPING_NUMBER_DERIVATIVES
    TYPE(MeshElementType), POINTER :: ELEMENT
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(ELEMENTS%elementsFinished) THEN
        CALL FlagError("Elements have been finished",ERR,ERROR,*999)
      ELSE
        IF(globalNumber>=1.AND.globalNumber<=ELEMENTS%numberOfElements) THEN
          IF(ASSOCIATED(BASIS)) THEN
            ELEMENT=>ELEMENTS%ELEMENTS(globalNumber)
            IF(ELEMENT%BASIS%numberOfNodes/=BASIS%numberOfNodes.OR. &
                & ELEMENT%BASIS%maximumNumberOfDerivatives/=BASIS%maximumNumberOfDerivatives) THEN
              !Allocate new user and global element nodes
              ALLOCATE(NEW_userElementNodes(BASIS%numberOfNodes),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate new user element nodes",ERR,ERROR,*999)
              ALLOCATE(NEW_GLOBAL_ELEMENT_NODES(BASIS%numberOfNodes),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate new user element nodes",ERR,ERROR,*999)
              ALLOCATE(NEW_userElementNodeVersions(BASIS%maximumNumberOfDerivatives, &
                & BASIS%numberOfNodes),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate element node versions",ERR,ERROR,*999)

              OVERLAPPING_NUMBER_NODES=MIN(BASIS%numberOfNodes,ELEMENT%BASIS%numberOfNodes)
              OVERLAPPING_NUMBER_DERIVATIVES=MIN(BASIS%maximumNumberOfDerivatives,ELEMENT%BASIS%maximumNumberOfDerivatives)

              !Set default values
              NEW_userElementNodeVersions=1
              NEW_userElementNodes(OVERLAPPING_NUMBER_NODES+1:)=0
              NEW_GLOBAL_ELEMENT_NODES(OVERLAPPING_NUMBER_NODES+1:)=0
              !Copy previous values
              NEW_userElementNodes(1:OVERLAPPING_NUMBER_NODES)=ELEMENT%userElementNodes(1:OVERLAPPING_NUMBER_NODES)
              NEW_GLOBAL_ELEMENT_NODES(1:OVERLAPPING_NUMBER_NODES)=ELEMENT%globalElementNodes(1:OVERLAPPING_NUMBER_NODES)
              NEW_userElementNodeVersions(1:OVERLAPPING_NUMBER_DERIVATIVES,1:OVERLAPPING_NUMBER_NODES)= &
                & ELEMENT%userElementNodeVersions(1:OVERLAPPING_NUMBER_DERIVATIVES,1:OVERLAPPING_NUMBER_NODES)

              !Replace arrays with new ones
              CALL MOVE_ALLOC(NEW_userElementNodeVersions,ELEMENT%userElementNodeVersions)
              CALL MOVE_ALLOC(NEW_userElementNodes,ELEMENT%userElementNodes)
              CALL MOVE_ALLOC(NEW_GLOBAL_ELEMENT_NODES,ELEMENT%globalElementNodes)
            ENDIF            
            ELEMENT%BASIS=>BASIS
          ELSE
            CALL FlagError("Basis is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="Global element number "//TRIM(NUMBER_TO_VSTRING(globalNumber,"*",ERR,ERROR))// &
            & " is invalid. The limits are 1 to "//TRIM(NUMBER_TO_VSTRING(ELEMENTS%numberOfElements,"*",ERR,ERROR))
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Elements is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET",ERR,ERROR)    
    RETURN 1
    
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET
  
  !
  !================================================================================================================================
  !

  !>Returns the adjacent element number for a mesh element identified by a global number. \todo specify by user number not global number \see OpenCMISS::Iron::cmfe_MeshElementsNo
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENT_GET(globalNumber,ELEMENTS,ADJACENT_ELEMENT_XI,ADJACENT_ELEMENT_NUMBER, &
    & ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: globalNumber !<The global number of the element to get the adjacent element for
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the elements of a mesh component from which to get the adjacent element from.
    INTEGER(INTG), INTENT(IN) :: ADJACENT_ELEMENT_XI !< The xi coordinate direction to get the adjacent element Note that -xiCoordinateDirection gives the adjacent element before the element in the xiCoordinateDirection'th direction and +xiCoordinateDirection gives the adjacent element after the element in the xiCoordinateDirection'th direction. The xiCoordinateDirection=0 index will give the information on the current element.
    INTEGER(INTG), INTENT(OUT) :: ADJACENT_ELEMENT_NUMBER !<On return, the adjacent element number in the specified xi coordinate direction. Return 0 if the specified element has no adjacent elements in the specified xi coordinate direction.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENT_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(.NOT.ELEMENTS%elementsFinished) THEN
        CALL FlagError("Elements have not been finished",ERR,ERROR,*999)
      ELSE
        IF(globalNumber>=1.AND.globalNumber<=ELEMENTS%numberOfElements) THEN
          IF(ADJACENT_ELEMENT_XI>=-ELEMENTS%ELEMENTS(globalNumber)%BASIS%numberOfXi .AND. &
            & ADJACENT_ELEMENT_XI<=ELEMENTS%ELEMENTS(globalNumber)%BASIS%numberOfXi) THEN
            IF(ELEMENTS%ELEMENTS(globalNumber)%adjacentElements(ADJACENT_ELEMENT_XI)%numberOfAdjacentElements > 0) THEN !\todo Currently returns only the first adjacent element for now as the python binding require the output array size of the adjacent element to be known a-prior. Add routine to first output number of adjacent elements and then loop over all adjacent elements
              ADJACENT_ELEMENT_NUMBER=ELEMENTS%ELEMENTS(globalNumber)%adjacentElements(ADJACENT_ELEMENT_XI)%adjacentElements(1)
            ELSE !Return 0 indicating the specified element has no adjacent elements in the specified xi coordinate direction.
              ADJACENT_ELEMENT_NUMBER=0
            ENDIF
          ELSE
            LOCAL_ERROR="The specified adjacent element xi is invalid. The supplied xi is "// &
            & TRIM(NUMBER_TO_VSTRING(ADJACENT_ELEMENT_XI,"*",ERR,ERROR))//" and needs to be >=-"// &
            & TRIM(NUMBER_TO_VSTRING(ELEMENTS%ELEMENTS(globalNumber)%BASIS%numberOfXi,"*",ERR,ERROR))//" and <="// &
            & TRIM(NUMBER_TO_VSTRING(ELEMENTS%ELEMENTS(globalNumber)%BASIS%numberOfXi,"*",ERR,ERROR))//"."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="Global element number "//TRIM(NUMBER_TO_VSTRING(globalNumber,"*",ERR,ERROR))// &
            & " is invalid. The limits are 1 to "//TRIM(NUMBER_TO_VSTRING(ELEMENTS%numberOfElements,"*",ERR,ERROR))
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Elements is not associated",ERR,ERROR,*999)
    ENDIF

    EXITS("MESH_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENT_GET")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENT_GET",ERR,ERROR)
    RETURN 1

  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENT_GET

  !
  !================================================================================================================================
  !

  !>Gets the element nodes for a mesh element identified by a given global number. \todo specify by user number not global number \see OpenCMISS::Iron::cmfe_MeshElementsNodesGet
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_GET(globalNumber,ELEMENTS,userElementNodes,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: globalNumber !<The global number of the element to set the nodes for
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the elements to set \todo before number?
    INTEGER(INTG), INTENT(OUT) :: userElementNodes(:) !<On return, userElementNodes(i). userElementNodes(i) is the i'th user node number for the element
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(.NOT.ELEMENTS%elementsFinished) THEN
        CALL FlagError("Elements have not been finished",ERR,ERROR,*999)
      ELSE
        IF(globalNumber>=1.AND.globalNumber<=ELEMENTS%numberOfElements) THEN
          IF(SIZE(userElementNodes,1)>=SIZE(ELEMENTS%ELEMENTS(globalNumber)%userElementNodes,1)) THEN
            userElementNodes=ELEMENTS%ELEMENTS(globalNumber)%userElementNodes
          ELSE
            LOCAL_ERROR="The size of userElementNodes is too small. The supplied size is "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(userElementNodes,1),"*",ERR,ERROR))//" and it needs to be >= "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(ELEMENTS%ELEMENTS(globalNumber)%userElementNodes,1),"*",ERR,ERROR))//"."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="Global element number "//TRIM(NUMBER_TO_VSTRING(globalNumber,"*",ERR,ERROR))// &
            & " is invalid. The limits are 1 to "//TRIM(NUMBER_TO_VSTRING(ELEMENTS%numberOfElements,"*",ERR,ERROR))
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Elements is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_GET")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_GET",ERR,ERROR)    
    RETURN 1
    
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_GET

  !
  !================================================================================================================================
  !

  !>Changes/sets the element nodes for a mesh element identified by a given global number. \todo specify by user number not global number \see OpenCMISS::Iron::cmfe_MeshElementsNodesSet
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(globalNumber,ELEMENTS,userElementNodes,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: globalNumber !<The global number of the element to set the nodes for
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the elements to set \todo before number?
    INTEGER(INTG), INTENT(IN) :: userElementNodes(:) !<userElementNodes(i). userElementNodes(i) is the i'th user node number for the element
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nn,NUMBER_OF_BAD_NODES,GLOBAL_NODE_NUMBER
    INTEGER(INTG), ALLOCATABLE :: GLOBAL_ELEMENT_NODES(:),BAD_NODES(:)
    LOGICAL :: ELEMENT_NODES_OK,NODE_EXISTS
    TYPE(InterfaceType), POINTER :: INTERFACE
    TYPE(MeshType), POINTER :: mesh
    TYPE(MeshComponentTopologyType), POINTER :: meshComponentTopology
    TYPE(NodesType), POINTER :: NODES
    TYPE(RegionType), POINTER :: PARENT_REGION,REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(ELEMENTS%elementsFinished) THEN
        CALL FlagError("Elements have been finished.",ERR,ERROR,*999)
      ELSE
        IF(globalNumber>=1.AND.globalNumber<=ELEMENTS%numberOfElements) THEN
          IF(SIZE(userElementNodes,1)==ELEMENTS%ELEMENTS(globalNumber)%BASIS%numberOfNodes) THEN
            meshComponentTopology=>elements%meshComponentTopology
            IF(ASSOCIATED(meshComponentTopology)) THEN
              mesh=>meshComponentTopology%mesh
              IF(ASSOCIATED(mesh)) THEN
                REGION=>MESH%REGION
                IF(ASSOCIATED(REGION)) THEN
                  NODES=>REGION%NODES
                ELSE
                  INTERFACE=>MESH%INTERFACE
                  IF(ASSOCIATED(INTERFACE)) THEN
                    NODES=>INTERFACE%NODES
                    PARENT_REGION=>INTERFACE%parentRegion
                    IF(.NOT.ASSOCIATED(PARENT_REGION)) CALL FlagError("Mesh interface has no parent region.",ERR,ERROR,*999)
                  ELSE
                    CALL FlagError("Elements mesh has no associated region or interface.",ERR,ERROR,*999)
                  ENDIF
                ENDIF
                IF(ASSOCIATED(NODES)) THEN
                  ELEMENT_NODES_OK=.TRUE.
                  ALLOCATE(GLOBAL_ELEMENT_NODES(ELEMENTS%ELEMENTS(globalNumber)%BASIS%numberOfNodes),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate global element nodes.",ERR,ERROR,*999)
                  ALLOCATE(BAD_NODES(ELEMENTS%ELEMENTS(globalNumber)%BASIS%numberOfNodes),STAT=ERR)
                  IF(ERR/=0) CALL FlagError("Could not allocate bad nodes.",ERR,ERROR,*999)
                  NUMBER_OF_BAD_NODES=0
                  DO nn=1,ELEMENTS%ELEMENTS(globalNumber)%BASIS%numberOfNodes
                    CALL NODE_CHECK_EXISTS(NODES,userElementNodes(nn),NODE_EXISTS,GLOBAL_NODE_NUMBER,ERR,ERROR,*999)
                    IF(NODE_EXISTS) THEN
                      GLOBAL_ELEMENT_NODES(nn)=GLOBAL_NODE_NUMBER
                    ELSE
                      NUMBER_OF_BAD_NODES=NUMBER_OF_BAD_NODES+1
                      BAD_NODES(NUMBER_OF_BAD_NODES)=userElementNodes(nn)
                      ELEMENT_NODES_OK=.FALSE.
                    ENDIF
                  ENDDO !nn
                  IF(ELEMENT_NODES_OK) THEN
                    ELEMENTS%ELEMENTS(globalNumber)%userElementNodes=userElementNodes
                    ELEMENTS%ELEMENTS(globalNumber)%globalElementNodes=GLOBAL_ELEMENT_NODES
                  ELSE
                    IF(NUMBER_OF_BAD_NODES==1) THEN
                      IF(ASSOCIATED(REGION)) THEN
                        LOCAL_ERROR="The element user node number of "//TRIM(NUMBER_TO_VSTRING(BAD_NODES(1),"*",ERR,ERROR))// &
                          & " is not defined in region "//TRIM(NUMBER_TO_VSTRING(REGION%userNumber,"*",ERR,ERROR))//"."
                      ELSE
                        LOCAL_ERROR="The element user node number of "//TRIM(NUMBER_TO_VSTRING(BAD_NODES(1),"*",ERR,ERROR))// &
                          & " is not defined in interface number "// &
                          & TRIM(NUMBER_TO_VSTRING(INTERFACE%userNumber,"*",ERR,ERROR))// &
                          & " of parent region number "//TRIM(NUMBER_TO_VSTRING(PARENT_REGION%userNumber,"*",ERR,ERROR))//"."
                      ENDIF
                    ELSE
                      LOCAL_ERROR="The element user node number of "//TRIM(NUMBER_TO_VSTRING(BAD_NODES(1),"*",ERR,ERROR))
                      DO nn=2,NUMBER_OF_BAD_NODES-1
                        LOCAL_ERROR=LOCAL_ERROR//","//TRIM(NUMBER_TO_VSTRING(BAD_NODES(nn),"*",ERR,ERROR))
                      ENDDO !nn
                      IF(ASSOCIATED(REGION)) THEN
                        LOCAL_ERROR=LOCAL_ERROR//" & "//TRIM(NUMBER_TO_VSTRING(BAD_NODES(NUMBER_OF_BAD_NODES),"*",ERR,ERROR))// &
                          & " are not defined in region number "//TRIM(NUMBER_TO_VSTRING(REGION%userNumber,"*",ERR,ERROR))//"."
                      ELSE
                        LOCAL_ERROR=LOCAL_ERROR//" & "//TRIM(NUMBER_TO_VSTRING(BAD_NODES(NUMBER_OF_BAD_NODES),"*",ERR,ERROR))// &
                          & " are not defined in interface number "// &
                          & TRIM(NUMBER_TO_VSTRING(INTERFACE%userNumber,"*",ERR,ERROR))//" of parent region number "// &
                          &  TRIM(NUMBER_TO_VSTRING(PARENT_REGION%userNumber,"*",ERR,ERROR))//"."
                      ENDIF
                    ENDIF
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  IF(ASSOCIATED(REGION)) THEN                   
                    CALL FlagError("The elements mesh region does not have any associated nodes.",ERR,ERROR,*999)
                  ELSE
                    CALL FlagError("The elements mesh interface does not have any associated nodes.",ERR,ERROR,*999) 
                  ENDIF
                ENDIF
              ELSE
                CALL FlagError("The mesh component topology mesh is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("The elements mesh component topology is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Number of element nodes does not match number of basis nodes for this element.",ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified global element number of "//TRIM(NUMBER_TO_VSTRING(globalNumber,"*",ERR,ERROR))// &
            & " is invalid. The global element number should be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(ELEMENTS%numberOfElements,"*",ERR,ERROR))//"."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Elements is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET",ERR,ERROR)    
    RETURN 1
    
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET

  !
  !================================================================================================================================
  !

  !>Changes/sets an element node's version for a mesh element identified by a given global number. \todo specify by user number not global number \see OpenCMISS::Iron::cmfe_MeshElementsNodesSet
  SUBROUTINE MeshElements_ElementNodeVersionSet(globalNumber,ELEMENTS,VERSION_NUMBER,DERIVATIVE_NUMBER, &
      & USER_ELEMENT_NODE_INDEX,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: globalNumber !<The global number of the element to set the nodes for
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the elements to set \todo before number?
    INTEGER(INTG), INTENT(IN) :: VERSION_NUMBER !<The version number of the specified element node to set.
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<The derivative number of the specified element node to set.
    INTEGER(INTG), INTENT(IN) :: USER_ELEMENT_NODE_INDEX !< The node index of the specified element node to set a version for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(InterfaceType), POINTER :: INTERFACE
    TYPE(MeshType), POINTER :: MESH
    TYPE(MeshComponentTopologyType), POINTER :: meshComponentTopology
    TYPE(NodesType), POINTER :: NODES
    TYPE(RegionType), POINTER :: PARENT_REGION,REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MeshElements_ElementNodeVersionSet",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(ELEMENTS%elementsFinished) THEN
        CALL FlagError("Elements have been finished.",ERR,ERROR,*999)
      ELSE
        IF(globalNumber>=1.AND.globalNumber<=ELEMENTS%numberOfElements) THEN
          IF(USER_ELEMENT_NODE_INDEX>=1.AND.USER_ELEMENT_NODE_INDEX<=ELEMENTS%ELEMENTS(globalNumber)%BASIS%numberOfNodes) THEN
            meshComponentTopology=>elements%meshComponentTopology
            IF(ASSOCIATED(meshComponentTopology)) THEN
              mesh=>meshComponentTopology%mesh
              IF(ASSOCIATED(mesh)) THEN
                REGION=>mesh%REGION
                IF(ASSOCIATED(REGION)) THEN
                  NODES=>REGION%NODES
                ELSE
                  INTERFACE=>mesh%INTERFACE
                  IF(ASSOCIATED(INTERFACE)) THEN
                    NODES=>INTERFACE%NODES
                    PARENT_REGION=>INTERFACE%parentRegion
                    IF(.NOT.ASSOCIATED(PARENT_REGION)) CALL FlagError("Mesh interface has no parent region.",ERR,ERROR,*999)
                  ELSE
                    CALL FlagError("Elements mesh has no associated region or interface.",ERR,ERROR,*999)
                  ENDIF
                ENDIF
                IF(DERIVATIVE_NUMBER>=1.AND.DERIVATIVE_NUMBER<=ELEMENTS%ELEMENTS(globalNumber)%BASIS% &
                  & numberOfDerivatives(USER_ELEMENT_NODE_INDEX)) THEN !Check if the specified derivative exists
                  IF(VERSION_NUMBER>=1) THEN !Check if the specified version is greater than 1
                    ELEMENTS%ELEMENTS(globalNumber)%userElementNodeVersions(DERIVATIVE_NUMBER,USER_ELEMENT_NODE_INDEX) & 
                      & = VERSION_NUMBER
                    !\todo : There is redunancy in userElementNodeVersions since it was allocated in MESH_TOPOLOGY_ELEMENTS_CREATE_START based on maximumNumberOfDerivatives for that elements basis:ALLOCATE(ELEMENTS%ELEMENTS(ne)%userElementNodeVersions(BASIS%maximumNumberOfDerivatives,BASIS%numberOfNodes),STAT=ERR)
                  ELSE
                    LOCAL_ERROR="The specified node version number of "//TRIM(NUMBER_TO_VSTRING(VERSION_NUMBER,"*", & 
                      & ERR,ERROR))//" is invalid. The element node index should be greater than 1."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The specified node derivative number of "//TRIM(NUMBER_TO_VSTRING(DERIVATIVE_NUMBER,"*", & 
                    & ERR,ERROR))//" is invalid. The element node derivative index should be between 1 and "// &
                    & TRIM(NUMBER_TO_VSTRING(ELEMENTS%ELEMENTS(globalNumber)%BASIS%numberOfDerivatives( &
                    & USER_ELEMENT_NODE_INDEX),"*",ERR,ERROR))//"."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("The mesh component topology mesh is not associated.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("The elements mesh component topology is not associated.",err,error,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The specified element node index of "//TRIM(NUMBER_TO_VSTRING(USER_ELEMENT_NODE_INDEX,"*",ERR,ERROR))// &
              & " is invalid. The element node index should be between 1 and "// &
              & TRIM(NUMBER_TO_VSTRING(ELEMENTS%ELEMENTS(globalNumber)%BASIS%numberOfNodes,"*",ERR,ERROR))//"."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified global element number of "//TRIM(NUMBER_TO_VSTRING(globalNumber,"*",ERR,ERROR))// &
            & " is invalid. The global element number should be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(ELEMENTS%numberOfElements,"*",ERR,ERROR))//"."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Elements is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MeshElements_ElementNodeVersionSet")
    RETURN
999 ERRORSEXITS("MeshElements_ElementNodeVersionSet",ERR,ERROR)    
    RETURN 1
    
  END SUBROUTINE MeshElements_ElementNodeVersionSet

  !
  !================================================================================================================================
  !

 !>Calculates the element numbers surrounding an element in a mesh topology.
  SUBROUTINE MeshTopology_ElementsAdjacentElementsCalculate(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: TOPOLOGY !<A pointer to the mesh topology to calculate the elements adjacent to elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: j,ne,ne1,nep1,ni,nic,nn,nn1,nn2,nn3,node_idx,np,np1,DUMMY_ERR,FACE_XI(2),FACE_XIC(3),nodePositionIndex(4)
    INTEGER(INTG) :: xi_direction,direction_index,xi_dir_check,xi_dir_search,NUMBER_NODE_MATCHES
    INTEGER(INTG) :: NUMBER_SURROUNDING,NUMBER_OF_NODES_XIC(4)
    INTEGER(INTG), ALLOCATABLE :: NODE_MATCHES(:),adjacentElements(:)
    LOGICAL :: XI_COLLAPSED,FACE_COLLAPSED(-3:3),SUBSET
    TYPE(LIST_TYPE), POINTER :: NODE_MATCH_LIST
    TYPE(LIST_PTR_TYPE) :: ADJACENT_ELEMENTS_LIST(-4:4)
    TYPE(BasisType), POINTER :: BASIS
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    NULLIFY(NODE_MATCH_LIST)
    DO nic=-4,4
      NULLIFY(ADJACENT_ELEMENTS_LIST(nic)%PTR)
    ENDDO !nic
    
    ENTERS("MeshTopology_ElementsAdjacentElementsCalculate",ERR,ERROR,*999)
    
    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%NODES)) THEN
        IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
          !Loop over the global elements in the mesh
          DO ne=1,TOPOLOGY%ELEMENTS%numberOfElements
            !%%%% first we initialize lists that are required to find the adjacent elements list
            BASIS=>TOPOLOGY%ELEMENTS%ELEMENTS(ne)%BASIS
            DO nic=-BASIS%numberOfXiCoordinates,BASIS%numberOfXiCoordinates
              NULLIFY(ADJACENT_ELEMENTS_LIST(nic)%PTR)
              CALL LIST_CREATE_START(ADJACENT_ELEMENTS_LIST(nic)%PTR,ERR,ERROR,*999)
              CALL LIST_DATA_TYPE_SET(ADJACENT_ELEMENTS_LIST(nic)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
              CALL LIST_INITIAL_SIZE_SET(ADJACENT_ELEMENTS_LIST(nic)%PTR,5,ERR,ERROR,*999)
              CALL LIST_CREATE_FINISH(ADJACENT_ELEMENTS_LIST(nic)%PTR,ERR,ERROR,*999)
            ENDDO !ni
            NUMBER_OF_NODES_XIC=1
            NUMBER_OF_NODES_XIC(1:BASIS%numberOfXiCoordinates)=BASIS%numberOfNodesXiC(1:BASIS%numberOfXiCoordinates)
            !Place the current element in the surrounding list
            CALL LIST_ITEM_ADD(ADJACENT_ELEMENTS_LIST(0)%PTR,TOPOLOGY%ELEMENTS%ELEMENTS(ne)%globalNumber,ERR,ERROR,*999)
            SELECT CASE(BASIS%TYPE)
            CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
              !Determine the collapsed "faces" if any
              nodePositionIndex=1
              !Loop over the face normals of the element
              DO ni=1,BASIS%numberOfXi
                !Determine the xi directions that lie in this xi direction
                FACE_XI(1)=OTHER_XI_DIRECTIONS3(ni,2,1)
                FACE_XI(2)=OTHER_XI_DIRECTIONS3(ni,3,1)
                !Reset the node_position_index in this xi direction
                nodePositionIndex(ni)=1
                !Loop over the two faces with this normal
                DO direction_index=-1,1,2
                  xi_direction=direction_index*ni
                  FACE_COLLAPSED(xi_direction)=.FALSE.
                  DO j=1,2
                    xi_dir_check=FACE_XI(j)
                    IF(xi_dir_check<=BASIS%numberOfXi) THEN
                      xi_dir_search=FACE_XI(3-j)
                      nodePositionIndex(xi_dir_search)=1
                      XI_COLLAPSED=.TRUE.
                      DO WHILE(nodePositionIndex(xi_dir_search)<=NUMBER_OF_NODES_XIC(xi_dir_search).AND.XI_COLLAPSED)
                        !Get the first local node along the xi check direction
                        nodePositionIndex(xi_dir_check)=1
                        nn1=BASIS%nodePositionIndexInv(nodePositionIndex(1),nodePositionIndex(2),nodePositionIndex(3),1)
                        !Get the second local node along the xi check direction
                        nodePositionIndex(xi_dir_check)=2
                        nn2=BASIS%nodePositionIndexInv(nodePositionIndex(1),nodePositionIndex(2),nodePositionIndex(3),1)
                        IF(nn1/=0.AND.nn2/=0) THEN
                          IF(TOPOLOGY%ELEMENTS%ELEMENTS(ne)%meshElementNodes(nn1)/= &
                            & TOPOLOGY%ELEMENTS%ELEMENTS(ne)%meshElementNodes(nn2)) XI_COLLAPSED=.TRUE.
                        ENDIF
                        nodePositionIndex(xi_dir_search)=nodePositionIndex(xi_dir_search)+1
                      ENDDO !xi_dir_search
                      IF(XI_COLLAPSED) FACE_COLLAPSED(xi_direction)=.TRUE.
                    ENDIF
                  ENDDO !j
                  nodePositionIndex(ni)=NUMBER_OF_NODES_XIC(ni)
                ENDDO !direction_index
              ENDDO !ni
              !Loop over the xi directions and calculate the surrounding elements
              DO ni=1,BASIS%numberOfXi
                !Determine the xi directions that lie in this xi direction
                FACE_XI(1)=OTHER_XI_DIRECTIONS3(ni,2,1)
                FACE_XI(2)=OTHER_XI_DIRECTIONS3(ni,3,1)
                !Loop over the two faces
                DO direction_index=-1,1,2
                  xi_direction=direction_index*ni
                  !Find nodes in the element on the appropriate face/line/point
                  NULLIFY(NODE_MATCH_LIST)
                  CALL LIST_CREATE_START(NODE_MATCH_LIST,ERR,ERROR,*999)
                  CALL LIST_DATA_TYPE_SET(NODE_MATCH_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)

                  CALL LIST_INITIAL_SIZE_SET(NODE_MATCH_LIST,16,ERR,ERROR,*999)
                  CALL LIST_CREATE_FINISH(NODE_MATCH_LIST,ERR,ERROR,*999)
                  IF(direction_index==-1) THEN
                    nodePositionIndex(ni)=1
                  ELSE
                    nodePositionIndex(ni)=NUMBER_OF_NODES_XIC(ni)
                  ENDIF
                  !If the face is collapsed then don't look in this xi direction. The exception is if the opposite face is also
                  !collpased. This may indicate that we have a funny element in non-rc coordinates that goes around the central
                  !axis back to itself
                  IF(FACE_COLLAPSED(xi_direction).AND..NOT.FACE_COLLAPSED(-xi_direction)) THEN
                    !Do nothing - the match lists are already empty
                  ELSE
                    !Find the nodes to match and add them to the node match list
                    DO nn1=1,NUMBER_OF_NODES_XIC(FACE_XI(1))
                      nodePositionIndex(FACE_XI(1))=nn1
                      DO nn2=1,NUMBER_OF_NODES_XIC(FACE_XI(2))
                        nodePositionIndex(FACE_XI(2))=nn2
                        nn=BASIS%nodePositionIndexInv(nodePositionIndex(1),nodePositionIndex(2),nodePositionIndex(3),1)
                        IF(nn/=0) THEN
                          np=TOPOLOGY%ELEMENTS%ELEMENTS(ne)%meshElementNodes(nn)
                          CALL LIST_ITEM_ADD(NODE_MATCH_LIST,np,ERR,ERROR,*999)
                        ENDIF
                      ENDDO !nn2
                    ENDDO !nn1
                  ENDIF
                  CALL LIST_REMOVE_DUPLICATES(NODE_MATCH_LIST,ERR,ERROR,*999)
                  CALL LIST_DETACH_AND_DESTROY(NODE_MATCH_LIST,NUMBER_NODE_MATCHES,NODE_MATCHES,ERR,ERROR,*999)
                  NUMBER_SURROUNDING=0
                  IF(NUMBER_NODE_MATCHES>0) THEN
                    !Find list of elements surrounding those nodes
                    np1=NODE_MATCHES(1)
                    DO nep1=1,TOPOLOGY%NODES%NODES(np1)%numberOfSurroundingElements
                      ne1=TOPOLOGY%NODES%NODES(np1)%surroundingElements(nep1)
                      IF(ne1/=ne) THEN !Don't want the current element
                        ! grab the nodes list for current and this surrouding elements
                        ! current face : NODE_MATCHES
                        ! candidate elem : TOPOLOGY%ELEMENTS%ELEMENTS(ne1)%meshElementNodes ! should this be GLOBAL_ELEMENT_NODES?
                        ! if all of current face belongs to the candidate element, we will have found the neighbour
                        CALL LIST_SUBSET_OF(NODE_MATCHES(1:NUMBER_NODE_MATCHES),TOPOLOGY%ELEMENTS%ELEMENTS(ne1)% &
                          & meshElementNodes,SUBSET,ERR,ERROR,*999)
                        IF(SUBSET) THEN
                          CALL LIST_ITEM_ADD(ADJACENT_ELEMENTS_LIST(xi_direction)%PTR,ne1,ERR,ERROR,*999)
                          NUMBER_SURROUNDING=NUMBER_SURROUNDING+1
                        ENDIF
                      ENDIF
                    ENDDO !nep1
                  ENDIF
                  IF(ALLOCATED(NODE_MATCHES)) DEALLOCATE(NODE_MATCHES)
                ENDDO !direction_index
              ENDDO !ni
            CASE(BASIS_SIMPLEX_TYPE)
              !Loop over the xi coordinates and calculate the surrounding elements
              DO nic=1,BASIS%numberOfXiCoordinates
                !Find the other coordinates of the face/line/point
                FACE_XIC(1)=OTHER_XI_DIRECTIONS4(nic,1)
                FACE_XIC(2)=OTHER_XI_DIRECTIONS4(nic,2)
                FACE_XIC(3)=OTHER_XI_DIRECTIONS4(nic,3)
                !Find nodes in the element on the appropriate face/line/point
                NULLIFY(NODE_MATCH_LIST)
                CALL LIST_CREATE_START(NODE_MATCH_LIST,ERR,ERROR,*999)
                CALL LIST_DATA_TYPE_SET(NODE_MATCH_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                CALL LIST_INITIAL_SIZE_SET(NODE_MATCH_LIST,16,ERR,ERROR,*999)
                CALL LIST_CREATE_FINISH(NODE_MATCH_LIST,ERR,ERROR,*999)
                nodePositionIndex(nic)=1 !Furtherest away from node with the nic'th coordinate
                !Find the nodes to match and add them to the node match list
                DO nn1=1,NUMBER_OF_NODES_XIC(FACE_XIC(1))
                  nodePositionIndex(FACE_XIC(1))=nn1
                  DO nn2=1,NUMBER_OF_NODES_XIC(FACE_XIC(2))
                    nodePositionIndex(FACE_XIC(2))=nn2
                    DO nn3=1,NUMBER_OF_NODES_XIC(FACE_XIC(3))
                      nodePositionIndex(FACE_XIC(3))=nn3
                      nn=BASIS%nodePositionIndexInv(nodePositionIndex(1),nodePositionIndex(2),nodePositionIndex(3), &
                        nodePositionIndex(4))
                      IF(nn/=0) THEN
                        np=TOPOLOGY%ELEMENTS%ELEMENTS(ne)%meshElementNodes(nn)
                        CALL LIST_ITEM_ADD(NODE_MATCH_LIST,np,ERR,ERROR,*999)
                      ENDIF
                    ENDDO !nn3
                  ENDDO !nn2
                ENDDO !nn1
                CALL LIST_REMOVE_DUPLICATES(NODE_MATCH_LIST,ERR,ERROR,*999)
                CALL LIST_DETACH_AND_DESTROY(NODE_MATCH_LIST,NUMBER_NODE_MATCHES,NODE_MATCHES,ERR,ERROR,*999)
                IF(NUMBER_NODE_MATCHES>0) THEN
                  !Find list of elements surrounding those nodes
                  DO node_idx=1,NUMBER_NODE_MATCHES
                    np1=NODE_MATCHES(node_idx)
                    DO nep1=1,TOPOLOGY%NODES%NODES(np1)%numberOfSurroundingElements
                      ne1=TOPOLOGY%NODES%NODES(np1)%surroundingElements(nep1)
                      IF(ne1/=ne) THEN !Don't want the current element
                        ! grab the nodes list for current and this surrouding elements
                        ! current face : NODE_MATCHES
                        ! candidate elem : TOPOLOGY%ELEMENTS%ELEMENTS(ne1)%meshElementNodes 
                        ! if all of current face belongs to the candidate element, we will have found the neighbour
                        CALL LIST_SUBSET_OF(NODE_MATCHES(1:NUMBER_NODE_MATCHES),TOPOLOGY%ELEMENTS%ELEMENTS(ne1)% &
                          & meshElementNodes,SUBSET,ERR,ERROR,*999)
                        IF(SUBSET) THEN
                          CALL LIST_ITEM_ADD(ADJACENT_ELEMENTS_LIST(nic)%PTR,ne1,ERR,ERROR,*999)
                        ENDIF
                      ENDIF
                    ENDDO !nep1
                  ENDDO !node_idx
                ENDIF
                IF(ALLOCATED(NODE_MATCHES)) DEALLOCATE(NODE_MATCHES)
              ENDDO !nic
            CASE(BASIS_SERENDIPITY_TYPE)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE(BASIS_AUXILLIARY_TYPE)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE(BASIS_B_SPLINE_TP_TYPE)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE(BASIS_FOURIER_LAGRANGE_HERMITE_TP_TYPE)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE(BASIS_EXTENDED_LAGRANGE_TP_TYPE)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The basis type of "//TRIM(NUMBER_TO_VSTRING(BASIS%TYPE,"*",ERR,ERROR))// &
                & " is invalid."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
            !Set the surrounding elements for this element
            ALLOCATE(TOPOLOGY%ELEMENTS%ELEMENTS(ne)%adjacentElements(-BASIS%numberOfXiCoordinates: &
              & BASIS%numberOfXiCoordinates),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate adjacent elements.",ERR,ERROR,*999)
            DO nic=-BASIS%numberOfXiCoordinates,BASIS%numberOfXiCoordinates
              CALL MESH_ADJACENT_ELEMENT_INITIALISE(TOPOLOGY%ELEMENTS%ELEMENTS(ne)%adjacentElements(nic),ERR,ERROR,*999)
              CALL LIST_REMOVE_DUPLICATES(ADJACENT_ELEMENTS_LIST(nic)%PTR,ERR,ERROR,*999)
              CALL LIST_DETACH_AND_DESTROY(ADJACENT_ELEMENTS_LIST(nic)%PTR,TOPOLOGY%ELEMENTS%ELEMENTS(ne)% &
                & adjacentElements(nic)%numberOfAdjacentElements,adjacentElements,ERR,ERROR,*999)
              ALLOCATE(TOPOLOGY%ELEMENTS%ELEMENTS(ne)%adjacentElements(nic)%adjacentElements(TOPOLOGY%ELEMENTS%ELEMENTS(ne)% &
                adjacentElements(nic)%numberOfAdjacentElements),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate element adjacent elements.",ERR,ERROR,*999)
              TOPOLOGY%ELEMENTS%ELEMENTS(ne)%adjacentElements(nic)%adjacentElements(1:TOPOLOGY%ELEMENTS%ELEMENTS(ne)% &
                adjacentElements(nic)%numberOfAdjacentElements) = adjacentElements(1:TOPOLOGY%ELEMENTS% &
                & ELEMENTS(ne)%adjacentElements(nic)%numberOfAdjacentElements)
              IF(ALLOCATED(adjacentElements)) DEALLOCATE(adjacentElements)
            ENDDO !nic
          ENDDO !ne           
        ELSE
          CALL FlagError("Mesh topology elements is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Mesh topology nodes is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh topology is not allocated.",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Number of elements = ",TOPOLOGY%ELEMENTS%numberOfElements,ERR,ERROR,*999)
      DO ne=1,TOPOLOGY%ELEMENTS%numberOfElements
        BASIS=>TOPOLOGY%ELEMENTS%ELEMENTS(ne)%BASIS
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Global element number : ",ne,ERR,ERROR,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of xi coordinates = ",BASIS%numberOfXiCoordinates, &
          & ERR,ERROR,*999)
        DO nic=-BASIS%numberOfXiCoordinates,BASIS%numberOfXiCoordinates
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Xi coordinate : ",nic,ERR,ERROR,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Number of adjacent elements = ", &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ne)%adjacentElements(nic)%numberOfAdjacentElements,ERR,ERROR,*999)
          IF(TOPOLOGY%ELEMENTS%ELEMENTS(ne)%adjacentElements(nic)%numberOfAdjacentElements>0) THEN
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,TOPOLOGY%ELEMENTS%ELEMENTS(ne)% &
              & adjacentElements(nic)%numberOfAdjacentElements,8,8,TOPOLOGY%ELEMENTS%ELEMENTS(ne)%adjacentElements(nic)% &
              & adjacentElements,'("        Adjacent elements :",8(X,I8))','(30x,8(X,I8))',ERR,ERROR,*999)
          ENDIF
        ENDDO !nic
      ENDDO !ne
    ENDIF
    
    EXITS("MeshTopology_ElementsAdjacentElementsCalculate")
    RETURN
999 IF(ALLOCATED(NODE_MATCHES)) DEALLOCATE(NODE_MATCHES)
    IF(ALLOCATED(adjacentElements)) DEALLOCATE(adjacentElements)
    IF(ASSOCIATED(NODE_MATCH_LIST)) CALL LIST_DESTROY(NODE_MATCH_LIST,DUMMY_ERR,DUMMY_ERROR,*998)
998 DO nic=-4,4
      IF(ASSOCIATED(ADJACENT_ELEMENTS_LIST(nic)%PTR)) CALL LIST_DESTROY(ADJACENT_ELEMENTS_LIST(nic)%PTR,DUMMY_ERR,DUMMY_ERROR,*997)
    ENDDO !ni
997 ERRORS("MeshTopology_ElementsAdjacentElementsCalculate",ERR,ERROR)
    EXITS("MeshTopology_ElementsAdjacentElementsCalculate")
    RETURN 1
    
  END SUBROUTINE MeshTopology_ElementsAdjacentElementsCalculate
  
  !
  !================================================================================================================================
  !

  !>Finalises the elements data structures for a mesh topology and deallocates any memory. \todo pass in elements
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_FINALISE(ELEMENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the mesh topology to finalise the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ne

    ENTERS("MESH_TOPOLOGY_ELEMENTS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      DO ne=1,ELEMENTS%numberOfElements
        CALL MESH_TOPOLOGY_ELEMENT_FINALISE(ELEMENTS%ELEMENTS(ne),ERR,ERROR,*999)
      ENDDO !ne
      DEALLOCATE(ELEMENTS%ELEMENTS)
      IF(ASSOCIATED(ELEMENTS%elementsTree)) CALL TREE_DESTROY(ELEMENTS%elementsTree,ERR,ERROR,*999)
      DEALLOCATE(ELEMENTS)
    ENDIF
 
    EXITS("MESH_TOPOLOGY_ELEMENTS_FINALISE")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the elements in a given mesh topology. \todo finalise on error
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: TOPOLOGY !<A pointer to the mesh topology to initialise the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("MESH_TOPOLOGY_ELEMENTS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
        CALL FlagError("Mesh already has topology elements associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%ELEMENTS,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology elements",ERR,ERROR,*999)
        TOPOLOGY%ELEMENTS%numberOfElements=0
        TOPOLOGY%ELEMENTS%meshComponentTopology=>TOPOLOGY
        NULLIFY(TOPOLOGY%ELEMENTS%ELEMENTS)
        NULLIFY(TOPOLOGY%ELEMENTS%elementsTree)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_TOPOLOGY_ELEMENTS_INITIALISE")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Initialises the elements in a given mesh topology. \todo finalise on error
  SUBROUTINE MESH_TOPOLOGY_DATA_POINTS_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: TOPOLOGY !<A pointer to the mesh topology to initialise the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("MESH_TOPOLOGY_DATA_POINTS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%dataPoints)) THEN
        CALL FlagError("Mesh already has topology data points associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%dataPoints,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate topology data points",ERR,ERROR,*999)
        TOPOLOGY%dataPoints%totalNumberOfProjectedData=0
        TOPOLOGY%dataPoints%meshComponentTopology=>topology
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_TOPOLOGY_DATA_POINTS_INITIALISE")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_DATA_POINTS_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESH_TOPOLOGY_DATA_POINTS_INITIALISE
  
  !
  !================================================================================================================================
  !

!!MERGE: ditto.
  
  !>Gets the user number for a global element identified by a given global number. \todo Check that the user number doesn't already exist. \see OpenCMISS::Iron::cmfe_MeshElementsUserNumberGet
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_NUMBER_GET(globalNumber,USER_NUMBER,ELEMENTS,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: globalNumber !<The global number of the elements to get.
    INTEGER(INTG), INTENT(OUT) :: USER_NUMBER !<The user number of the element to get
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<On return, a pointer to the elements to get the user number for \todo This should be the first parameter.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESH_TOPOLOGY_ELEMENTS_NUMBER_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(ELEMENTS%elementsFinished) THEN
        CALL FlagError("Elements have been finished",ERR,ERROR,*999)
      ELSE
        IF(globalNumber>=1.AND.globalNumber<=ELEMENTS%numberOfElements) THEN
          USER_NUMBER=ELEMENTS%ELEMENTS(globalNumber)%userNumber
        ELSE
          LOCAL_ERROR="Global element number "//TRIM(NUMBER_TO_VSTRING(globalNumber,"*",ERR,ERROR))// &
            & " is invalid. The limits are 1 to "//TRIM(NUMBER_TO_VSTRING(ELEMENTS%numberOfElements,"*",ERR,ERROR))
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Elements is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESH_TOPOLOGY_ELEMENTS_NUMBER_GET")
    RETURN
999 ERRORSEXITS("MESH_TOPOLOGY_ELEMENTS_NUMBER_GET",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_NUMBER_GET

  !
  !================================================================================================================================
  !

  !>Returns the user number for a global element identified by a given global number. \see OpenCMISS::Iron::cmfe_MeshElementsUserNumberGet
  SUBROUTINE MeshElements_ElementUserNumberGet(globalNumber,USER_NUMBER,ELEMENTS,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: globalNumber !<The global number of the elements to get.
    INTEGER(INTG), INTENT(OUT) :: USER_NUMBER !<The user number of the element to get
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the elements to set the user number for \todo This should be the first parameter.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MeshElements_ElementUserNumberGet",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(ELEMENTS%elementsFinished) THEN
        IF(globalNumber>=1.AND.globalNumber<=ELEMENTS%numberOfElements) THEN          
          USER_NUMBER=ELEMENTS%ELEMENTS(globalNumber)%userNumber
        ELSE
          LOCAL_ERROR="Global element number "//TRIM(NUMBER_TO_VSTRING(globalNumber,"*",ERR,ERROR))// &
            & " is invalid. The limits are 1 to "//TRIM(NUMBER_TO_VSTRING(ELEMENTS%numberOfElements,"*",ERR,ERROR))//"."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Elements have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Elements is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MeshElements_ElementUserNumberGet")
    RETURN
999 ERRORSEXITS("MeshElements_ElementUserNumberGet",ERR,ERROR)
    RETURN 1

   
  END SUBROUTINE MeshElements_ElementUserNumberGet

  !
  !================================================================================================================================
  !

  !>Changes/sets the user number for a global element identified by a given global number. \see OpenCMISS::Iron::cmfe_MeshElementsUserNumberSet
  SUBROUTINE MeshElements_ElementUserNumberSet(globalNumber,USER_NUMBER,ELEMENTS,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: globalNumber !<The global number of the elements to set.
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the element to set
    TYPE(MeshElementsType), POINTER :: ELEMENTS !<A pointer to the elements to set the user number for \todo This should be the first parameter.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code

    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: GLOBAL_ELEMENT_NUMBER,INSERT_STATUS
    TYPE(TREE_NODE_TYPE), POINTER :: TREE_NODE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MeshElements_ElementUserNumberSet",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(ELEMENTS%elementsFinished) THEN
        CALL FlagError("Elements have been finished.",ERR,ERROR,*999)
      ELSE
        IF(globalNumber>=1.AND.globalNumber<=ELEMENTS%numberOfElements) THEN
          NULLIFY(TREE_NODE)
          CALL TREE_SEARCH(ELEMENTS%elementsTree,USER_NUMBER,TREE_NODE,ERR,ERROR,*999)
          IF(ASSOCIATED(TREE_NODE)) THEN
            CALL TREE_NODE_VALUE_GET(ELEMENTS%elementsTree,TREE_NODE,GLOBAL_ELEMENT_NUMBER,ERR,ERROR,*999)
            LOCAL_ERROR="Element user number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
              & " is already used by global element number "// &
              & TRIM(NUMBER_TO_VSTRING(GLOBAL_ELEMENT_NUMBER,"*",ERR,ERROR))//"."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            CALL TREE_ITEM_DELETE(ELEMENTS%elementsTree,ELEMENTS%ELEMENTS(globalNumber)%userNumber,ERR,ERROR,*999)
            CALL TREE_ITEM_INSERT(ELEMENTS%elementsTree,USER_NUMBER,globalNumber,INSERT_STATUS,ERR,ERROR,*999)
            ELEMENTS%ELEMENTS(globalNumber)%userNumber=USER_NUMBER
          ENDIF
        ELSE
          LOCAL_ERROR="Global element number "//TRIM(NUMBER_TO_VSTRING(globalNumber,"*",ERR,ERROR))// &
            & " is invalid. The limits are 1 to "//TRIM(NUMBER_TO_VSTRING(ELEMENTS%numberOfElements,"*",ERR,ERROR))//"."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Elements is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MeshElements_ElementUserNumberSet")
    RETURN
999 ERRORSEXITS("MeshElements_ElementUserNumberSet",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE MeshElements_ElementUserNumberSet
  
  !
  !================================================================================================================================
  !

  !>Changes/sets the user numbers for all elements.
  SUBROUTINE MeshTopology_ElementsUserNumbersAllSet(elements,userNumbers,err,error,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: elements !<A pointer to the elements to set all the user numbers for 
    INTEGER(INTG), INTENT(IN) :: userNumbers(:) !<The user numbers to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementIdx,insertStatus
    TYPE(TREE_TYPE), POINTER :: newElementsTree
    TYPE(VARYING_STRING) :: localError

    NULLIFY(newElementsTree)

    ENTERS("MeshTopology_ElementsUserNumbersAllSet",err,error,*999)

    IF(ASSOCIATED(elements)) THEN
      IF(elements%elementsFinished) THEN
        CALL FlagError("Elements have been finished.",err,error,*999)
      ELSE
        IF(elements%numberOfElements==SIZE(userNumbers,1)) THEN
          !Check the users numbers to ensure that there are no duplicates          
          CALL TREE_CREATE_START(newElementsTree,err,error,*999)
          CALL TREE_INSERT_TYPE_SET(newElementsTree,TREE_NO_DUPLICATES_ALLOWED,err,error,*999)
          CALL TREE_CREATE_FINISH(newElementsTree,err,error,*999)
          DO elementIdx=1,elements%numberOfElements
            CALL TREE_ITEM_INSERT(newElementsTree,userNumbers(elementIdx),elementIdx,insertStatus,err,error,*999)
            IF(insertStatus/=TREE_NODE_INSERT_SUCESSFUL) THEN
              localError="The specified user number of "//TRIM(NumberToVstring(userNumbers(elementIdx),"*",err,error))// &
                & " for global element number "//TRIM(NUMBER_TO_VSTRING(elementIdx,"*",err,error))// &
                & " is a duplicate. The user element numbers must be unique."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ENDDO !elementIdx
          CALL TREE_DESTROY(elements%elementsTree,err,error,*999)
          elements%elementsTree=>newElementsTree
          NULLIFY(newElementsTree)
          DO elementIdx=1,elements%numberOfElements
            elements%ELEMENTS(elementIdx)%globalNumber=elementIdx
            elements%ELEMENTS(elementIdx)%userNumber=userNumbers(elementIdx)
          ENDDO !elementIdx
        ELSE
          localError="The number of specified element user numbers ("// &
            TRIM(NumberToVstring(SIZE(userNumbers,1),"*",err,error))// &
            ") does not match number of elements ("// &
            TRIM(NumberToVstring(elements%numberOfElements,"*",err,error))//")."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Elements is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopology_ElementsUserNumbersAllSet")
    RETURN
999 IF(ASSOCIATED(newElementsTree)) CALL TREE_DESTROY(newElementsTree,err,error,*998)
998 ERRORSEXITS("MeshTopology_ElementsUserNumbersAllSet",err,error)    
    RETURN 1
   
  END SUBROUTINE MeshTopology_ElementsUserNumbersAllSet
  
  !
  !================================================================================================================================
  !

  !>Calculates the data points in the given mesh topology.
  SUBROUTINE MeshTopology_DataPointsCalculateProjection(mesh,dataProjection,err,error,*)
  
    !Argument variables
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh topology to calcualte the data projection for
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points
    TYPE(MeshDataPointsType), POINTER :: dataPointsTopology
    TYPE(DataProjectionResultType), POINTER :: dataProjectionResult
    TYPE(MeshElementsType), POINTER :: elements
    INTEGER(INTG) :: dataPointIdx,elementIdx,exitTag,countIdx,projectionNumber,globalCountIdx,elementNumber

    ENTERS("MeshTopology_DataPointsCalculateProjection",ERR,ERROR,*999)

    IF(ASSOCIATED(mesh)) THEN
      IF(dataProjection%dataProjectionFinished) THEN 
        dataPoints=>dataProjection%dataPoints
        !Default the first mesh component topology to contain data points ! \TODO: need to be changed once the data points topology is moved under meshTopologyType.
        dataPointsTopology=>mesh%topology(1)%PTR%dataPoints
        !Extract the global number of the data projection 
        projectionNumber=dataProjection%globalNumber
        !Hard code the first mesh component since element topology is the same for all mesh components
        !\TODO: need to be changed once the elements topology is moved under meshTopologyType.
        elements=>mesh%TOPOLOGY(1)%PTR%ELEMENTS
        ALLOCATE(dataPointsTopology%elementDataPoint(elements%numberOfElements),STAT=ERR)     
        IF(ERR/=0) CALL FlagError("Could not allocate data points topology element.",ERR,ERROR,*999)
        DO elementIdx=1,elements%numberOfElements
          dataPointsTopology%elementDataPoint(elementIdx)%elementNumber=elements%ELEMENTS(elementIdx)%globalNumber
          dataPointsTopology%elementDataPoint(elementIdx)%numberOfProjectedData=0
        ENDDO        
        !Calculate number of projected data points on an element
        DO dataPointIdx=1,dataPoints%numberOfDataPoints
          dataProjectionResult=>dataProjection%dataProjectionResults(dataPointIdx)
          elementNumber=dataProjectionResult%elementNumber
          exitTag=dataProjectionResult%exitTag
          IF(exitTag/=DATA_PROJECTION_CANCELLED) THEN
            DO elementIdx=1,elements%numberOfElements
              IF(dataPointsTopology%elementDataPoint(elementIdx)%elementNumber==elementNumber) THEN
                dataPointsTopology%elementDataPoint(elementIdx)%numberOfProjectedData= &
                  & dataPointsTopology%elementDataPoint(elementIdx)%numberOfProjectedData+1;
              ENDIF
            ENDDO !elementIdx
          ENDIF
        ENDDO !dataPointIdx      
        !Allocate memory to store data indices and initialise them to be zero   
        DO elementIdx=1,elements%numberOfElements
          ALLOCATE(dataPointsTopology%elementDataPoint(elementIdx)%dataIndices(dataPointsTopology% &
            & elementDataPoint(elementIdx)%numberOfProjectedData),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate data points topology element data points.",ERR,ERROR,*999)
          DO countIdx=1,dataPointsTopology%elementDataPoint(elementIdx)%numberOfProjectedData
            dataPointsTopology%elementDataPoint(elementIdx)%dataIndices(countIdx)%userNumber=0
            dataPointsTopology%elementDataPoint(elementIdx)%dataIndices(countIdx)%globalNumber=0
          ENDDO
        ENDDO     
        !Record the indices of the data that projected on the elements 
        globalCountIdx=0
        dataPointsTopology%totalNumberOfProjectedData=0
        DO dataPointIdx=1,dataPoints%numberOfDataPoints 
          dataProjectionResult=>dataProjection%dataProjectionResults(dataPointIdx)
          elementNumber=dataProjectionResult%elementNumber
          exitTag=dataProjectionResult%exitTag
          IF(exitTag/=DATA_PROJECTION_CANCELLED) THEN
            DO elementIdx=1,elements%numberOfElements
              countIdx=1         
              IF(dataPointsTopology%elementDataPoint(elementIdx)%elementNumber==elementNumber) THEN
                globalCountIdx=globalCountIdx+1
                !Find the next data point index in this element
                DO WHILE(dataPointsTopology%elementDataPoint(elementIdx)%dataIndices(countIdx)%globalNumber/=0)
                  countIdx=countIdx+1
                ENDDO
                dataPointsTopology%elementDataPoint(elementIdx)%dataIndices(countIdx)%userNumber=dataPointIdx
                dataPointsTopology%elementDataPoint(elementIdx)%dataIndices(countIdx)%globalNumber=dataPointIdx!globalCountIdx (used this if only projected data are taken into account)
                dataPointsTopology%totalNumberOfProjectedData=dataPointsTopology%totalNumberOfProjectedData+1
              ENDIF
            ENDDO !elementIdx
          ENDIF
        ENDDO !dataPointIdx
        !Allocate memory to store total data indices in ascending order and element map
        ALLOCATE(dataPointsTopology%dataPoints(dataPointsTopology%totalNumberOfProjectedData),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate data points topology data points.",ERR,ERROR,*999)
        !The global number for the data points will be looping through elements.
        countIdx=1  
        DO elementIdx=1,elements%numberOfElements
          DO dataPointIdx=1,dataPointsTopology%elementDataPoint(elementIdx)%numberOfProjectedData
            dataPointsTopology%dataPoints(countIdx)%userNumber=dataPointsTopology%elementDataPoint(elementIdx)% &
              & dataIndices(dataPointIdx)%userNumber
             dataPointsTopology%dataPoints(countIdx)%globalNumber=dataPointsTopology%elementDataPoint(elementIdx)% &
              & dataIndices(dataPointIdx)%globalNumber
             dataPointsTopology%dataPoints(countIdx)%elementNumber=dataPointsTopology%elementDataPoint(elementIdx)% &
               & elementNumber
             countIdx=countIdx+1
          ENDDO !dataPointIdx
        ENDDO !elementIdx                      
      ELSE
        CALL FlagError("Data projection is not finished.",err,error,*999)
      ENDIF     
    ELSE
      CALL FlagError("Mesh is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopology_DataPointsCalculateProjection")
    RETURN
999 ERRORSEXITS("MeshTopology_DataPointsCalculateProjection",err,error)
    RETURN 1
  END SUBROUTINE MeshTopology_DataPointsCalculateProjection

  !
  !================================================================================================================================
  !

  !>Finalises the topology in the given mesh. \todo pass in the mesh topology
  SUBROUTINE MeshTopology_Finalise(mesh,err,error,*)

    !Argument variables
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh to finalise the topology for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx

    ENTERS("MeshTopology_Finalise",err,error,*999)

    IF(ASSOCIATED(mesh)) THEN
      DO componentIdx=1,mesh%NUMBER_OF_COMPONENTS
        CALL MeshTopology_ComponentFinalise(mesh%topology(componentIdx)%ptr,err,error,*999)
      ENDDO !componentIdx
      DEALLOCATE(mesh%topology)
    ELSE
      CALL FlagError("Mesh is not associated.",err,error,*999)
    ENDIF
 
    EXITS("MeshTopology_Finalise")
    RETURN
999 ERRORSEXITS("MeshTopology_Finalise",err,error)
    RETURN 1
   
  END SUBROUTINE MeshTopology_Finalise

  !
  !================================================================================================================================
  !

  !>Finalises the topology in the given mesh. 
  SUBROUTINE MeshTopology_ComponentFinalise(meshComponent,err,error,*)

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: meshComponent !<A pointer to the mesh component to finalise the topology for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshTopology_ComponentFinalise",err,error,*999)

    IF(ASSOCIATED(meshComponent)) THEN
      CALL MeshTopology_NodesFinalise(meshComponent%nodes,err,error,*999)
      CALL MESH_TOPOLOGY_ELEMENTS_FINALISE(meshComponent%elements,err,error,*999)
      CALL MeshTopology_DofsFinalise(meshComponent%dofs,err,error,*999)
      DEALLOCATE(meshComponent)
    ENDIF
 
    EXITS("MeshTopology_ComponentFinalise")
    RETURN
999 ERRORSEXITS("MeshTopology_ComponentFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE MeshTopology_ComponentFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the topology for a given mesh. \todo finalise on error
  SUBROUTINE MeshTopology_Initialise(mesh,err,error,*)

    !Argument variables
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh to initialise the mesh topology for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx
    
    ENTERS("MeshTopology_Initialise",err,error,*999)

    IF(ASSOCIATED(mesh)) THEN
      IF(ASSOCIATED(mesh%topology)) THEN
        CALL FlagError("Mesh already has topology associated.",err,error,*999)
      ELSE
        !Allocate mesh topology
        ALLOCATE(mesh%topology(mesh%NUMBER_OF_COMPONENTS),STAT=err)
        IF(err/=0) CALL FlagError("Mesh topology could not be allocated.",err,error,*999)
        DO componentIdx=1,mesh%NUMBER_OF_COMPONENTS
          ALLOCATE(mesh%topology(componentIdx)%ptr,STAT=err)
          IF(err/=0) CALL FlagError("Mesh topology component could not be allocated.",err,error,*999)
          mesh%topology(componentIdx)%ptr%mesh=>mesh
          NULLIFY(mesh%topology(componentIdx)%ptr%elements)
          NULLIFY(mesh%topology(componentIdx)%ptr%nodes)
          NULLIFY(mesh%topology(componentIdx)%ptr%dofs)
          NULLIFY(mesh%topology(componentIdx)%ptr%dataPoints)
          !Initialise the topology components
          CALL MESH_TOPOLOGY_ELEMENTS_INITIALISE(mesh%topology(componentIdx)%ptr,err,error,*999)
          CALL MeshTopology_NodesInitialise(mesh%topology(componentIdx)%ptr,err,error,*999)
          CALL MeshTopology_DofsInitialise(mesh%topology(componentIdx)%ptr,err,error,*999)
          CALL MESH_TOPOLOGY_DATA_POINTS_INITIALISE(mesh%topology(componentIdx)%ptr,err,error,*999)
        ENDDO !componentIdx
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopology_Initialise")
    RETURN
999 ERRORSEXITS("MeshTopology_Initialise",err,error)
    RETURN 1
  END SUBROUTINE MeshTopology_Initialise
  
  !
  !================================================================================================================================
  !

  !>Checks that a user element number exists in a mesh component. 
  SUBROUTINE MeshTopology_ElementCheckExistsMesh(mesh,meshComponentNumber,userElementNumber,elementExists,globalElementNumber, &
    & err,error,*)

    !Argument variables
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh to check the element exists on
    INTEGER(INTG), INTENT(IN) :: meshComponentNumber !<The mesh component to check the element exits on
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number to check if it exists
    LOGICAL, INTENT(OUT) :: elementExists !<On exit, is .TRUE. if the element user number exists in the mesh component, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: globalElementNumber !<On exit, if the element exists the global number corresponding to the user element number. If the element does not exist then global number will be 0.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(MeshElementsType), POINTER :: elements

    NULLIFY(elements)
    
    ENTERS("MeshTopology_ElementCheckExistsMesh",err,error,*999)

    IF(ASSOCIATED(mesh)) THEN
      IF(mesh%meshFinished) THEN
        CALL MESH_TOPOLOGY_ELEMENTS_GET(mesh,meshComponentNumber,elements,err,error,*999)
        CALL MeshTopology_ElementCheckExistsMeshElements(elements,userElementNumber,elementExists,globalElementNumber, &
          & err,error,*999)
      ELSE
        CALL FlagError("Mesh has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopology_ElementCheckExistsMesh")
    RETURN
999 ERRORSEXITS("MeshTopology_ElementCheckExistsMesh",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_ElementCheckExistsMesh
  
  !
  !================================================================================================================================
  !

  !>Checks that a user element number exists in a mesh elements. 
  SUBROUTINE MeshTopology_ElementCheckExistsMeshElements(meshElements,userElementNumber,elementExists,globalElementNumber, &
    & err,error,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: meshElements !<A pointer to the mesh elements to check the element exists on
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number to check if it exists
    LOGICAL, INTENT(OUT) :: elementExists !<On exit, is .TRUE. if the element user number exists in the mesh elements, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: globalElementNumber !<On exit, if the element exists the global number corresponding to the user element number. If the element does not exist then global number will be 0.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(TREE_NODE_TYPE), POINTER :: treeNode
    
    ENTERS("MeshTopology_ElementCheckExistsMeshElements",err,error,*999)

    elementExists=.FALSE.
    globalElementNumber=0
    IF(ASSOCIATED(meshElements)) THEN
      NULLIFY(treeNode)
      CALL TREE_SEARCH(meshElements%elementsTree,userElementNumber,treeNode,err,error,*999)
      IF(ASSOCIATED(treeNode)) THEN
        CALL TREE_NODE_VALUE_GET(meshElements%elementsTree,treeNode,globalElementNumber,err,error,*999)
        elementExists=.TRUE.
      ENDIF
    ELSE
      CALL FlagError("Mesh elements is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopology_ElementCheckExistsMeshElements")
    RETURN
999 ERRORSEXITS("MeshTopology_ElementCheckExistsMeshElements",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_ElementCheckExistsMeshElements
  
  !
  !================================================================================================================================
  !

  !>Gets a mesh element number that corresponds to a user element number from mesh element. An error will be raised if the user element number does not exist.
  SUBROUTINE MeshTopology_ElementGetMeshElements(meshElements,userElementNumber,globalElementNumber,err,error,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: meshElements !<A pointer to the mesh elements to get the element for
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number to get
    INTEGER(INTG), INTENT(OUT) :: globalElementNumber !<On exit, the global number corresponding to the user element number. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    LOGICAL :: elementExists
    TYPE(MeshType), POINTER :: mesh
    TYPE(MeshComponentTopologyType), POINTER :: meshComponentTopology
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("MeshTopology_ElementGetMeshElements",err,error,*999)

    CALL MeshTopology_ElementCheckExistsMeshElements(meshElements,userElementNumber,elementExists,globalElementNumber, &
      & err,error,*999)
    IF(.NOT.elementExists) THEN
      meshComponentTopology=>meshElements%meshComponentTopology
      IF(ASSOCIATED(meshComponentTopology)) THEN
        mesh=>meshComponentTopology%mesh
        IF(ASSOCIATED(mesh)) THEN
          localError="The user element number "//TRIM(NumberToVString(userElementNumber,"*",err,error))// &
            & " does not exist in mesh number "//TRIM(NumberToVString(mesh%userNumber,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ELSE
          localError="The user element number "//TRIM(NumberToVString(userElementNumber,"*",err,error))//" does not exist."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        localError="The user element number "//TRIM(NumberToVString(userElementNumber,"*",err,error))//" does not exist."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDIF
   
    EXITS("MeshTopology_ElementGetMeshElements")
    RETURN
999 globalElementNumber=0
    ERRORSEXITS("MeshTopology_ElementGetMeshElements",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_ElementGetMeshElements
  
  !
  !================================================================================================================================
  !
  
  !>Returns if the element in a mesh is on the boundary or not
  SUBROUTINE MeshTopology_ElementOnBoundaryGet(meshElements,userNumber,onBoundary,err,error,*)

    !Argument variables
    TYPE(MeshElementsType), POINTER :: meshElements !<A pointer to the mesh element containing the element to get the boundary type for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user element number to get the boundary type for
    INTEGER(INTG), INTENT(OUT) :: onBoundary !<On return, the boundary type of the specified user element number.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: globalNumber

    ENTERS("MeshTopology_ElementOnBoundaryGet",err,error,*999)

    IF(.NOT.ASSOCIATED(meshElements)) CALL FlagError("Mesh elements is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(meshElements%elements)) CALL FlagError("Mesh elements elements is not associated.",err,error,*999)
    
    CALL MeshTopology_ElementGet(meshElements,userNumber,globalNumber,err,error,*999)
    IF(meshElements%elements(globalNumber)%boundaryElement) THEN
      onBoundary=MESH_ON_DOMAIN_BOUNDARY
    ELSE
      onBoundary=MESH_OFF_DOMAIN_BOUNDARY
    ENDIF
    
    EXITS("MeshTopology_ElementOnBoundaryGet")
    RETURN
999 onBoundary=MESH_OFF_DOMAIN_BOUNDARY
    ERRORSEXITS("MeshTopology_ElementOnBoundaryGet",err,error)    
    RETURN 1
   
  END SUBROUTINE MeshTopology_ElementOnBoundaryGet
  
  !
  !================================================================================================================================
  !

  !>Checks that a user node number exists in a mesh component. 
  SUBROUTINE MeshTopology_NodeCheckExistsMesh(mesh,meshComponentNumber,userNodeNumber,nodeExists,meshNodeNumber,err,error,*)

    !Argument variables
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh to check the node exists on
    INTEGER(INTG), INTENT(IN) :: meshComponentNumber !<The mesh component to check the node exits on
    INTEGER(INTG), INTENT(IN) :: userNodeNumber !<The user node number to check if it exists
    LOGICAL, INTENT(OUT) :: nodeExists !<On exit, is .TRUE. if the node user number exists in the mesh component, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: meshNodeNumber !<On exit, if the node exists the mesh number corresponding to the user node number. If the node does not exist then mesh number will be 0.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: globalNodeNumber
    TYPE(MeshNodesType), POINTER :: meshNodes
    TYPE(NodesType), POINTER :: nodes
    TYPE(RegionType), POINTER :: region
    TYPE(TREE_NODE_TYPE), POINTER :: treeNode

    NULLIFY(meshNodes)
    NULLIFY(nodes)
    NULLIFY(region)
    
    ENTERS("MeshTopology_NodeCheckExistsMesh",err,error,*999)

    nodeExists=.FALSE.
    meshNodeNumber=0
    IF(ASSOCIATED(mesh)) THEN
      IF(mesh%meshFinished) THEN
        CALL MeshTopology_NodesGet(mesh,meshComponentNumber,meshNodes,err,error,*999)
        CALL Mesh_RegionGet(mesh,region,err,error,*999)
        nodes=>region%nodes
        IF(ASSOCIATED(nodes)) THEN
          CALL NODE_CHECK_EXISTS(nodes,userNodeNumber,nodeExists,globalNodeNumber,err,error,*999)
          NULLIFY(treeNode)
          CALL TREE_SEARCH(meshNodes%nodesTree,globalNodeNumber,treeNode,err,error,*999)
          IF(ASSOCIATED(treeNode)) THEN
            CALL TREE_NODE_VALUE_GET(meshNodes%nodesTree,treeNode,meshNodeNumber,err,error,*999)
            nodeExists=.TRUE.
          ENDIF
        ELSE
          CALL FlagError("Region nodes is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Mesh has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopology_NodeCheckExistsMesh")
    RETURN
999 ERRORSEXITS("MeshTopology_NodeCheckExistsMesh",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_NodeCheckExistsMesh
  
  !
  !================================================================================================================================
  !

  !>Checks that a user node number exists in a mesh nodes. 
  SUBROUTINE MeshTopology_NodeCheckExistsMeshNodes(meshNodes,userNodeNumber,nodeExists,meshNodeNumber,err,error,*)

    !Argument variables
    TYPE(MeshNodesType), POINTER :: meshNodes !<A pointer to the mesh nodes to check the node exists on
    INTEGER(INTG), INTENT(IN) :: userNodeNumber !<The user node number to check if it exists
    LOGICAL, INTENT(OUT) :: nodeExists !<On exit, is .TRUE. if the node user number exists in the mesh nodes, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: meshNodeNumber !<On exit, if the node exists the mesh number corresponding to the user node number. If the node does not exist then mesh number will be 0.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: globalNodeNumber
    TYPE(MeshType), POINTER :: mesh
    TYPE(MeshComponentTopologyType), POINTER :: meshComponentTopology
    TYPE(NodesType), POINTER :: nodes
    TYPE(RegionType), POINTER :: region
    TYPE(TREE_NODE_TYPE), POINTER :: treeNode

    NULLIFY(nodes)
    NULLIFY(region)
    
    ENTERS("MeshTopology_NodeCheckExistsMeshNodes",err,error,*999)

    nodeExists=.FALSE.
    meshNodeNumber=0
    IF(ASSOCIATED(meshNodes)) THEN
      meshComponentTopology=>meshNodes%meshComponentTopology
      IF(ASSOCIATED(meshComponentTopology)) THEN
        mesh=>meshComponentTopology%mesh
        IF(ASSOCIATED(mesh)) THEN
          IF(mesh%meshFinished) THEN
            CALL Mesh_RegionGet(mesh,region,err,error,*999)
            CALL Region_NodesGet(region,nodes,err,error,*999)
            CALL NODE_CHECK_EXISTS(nodes,userNodeNumber,nodeExists,globalNodeNumber,err,error,*999)
            NULLIFY(treeNode)
            CALL TREE_SEARCH(meshNodes%nodesTree,globalNodeNumber,treeNode,err,error,*999)
            IF(ASSOCIATED(treeNode)) THEN
              CALL TREE_NODE_VALUE_GET(meshNodes%nodesTree,treeNode,meshNodeNumber,err,error,*999)
              nodeExists=.TRUE.
            ENDIF
         ELSE
            CALL FlagError("Mesh has not been finished.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Mesh component topology mesh is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Mesh nodes mesh component topology is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh nodes is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopology_NodeCheckExistsMeshNodes")
    RETURN
999 ERRORSEXITS("MeshTopology_NodeCheckExistsMeshNodes",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_NodeCheckExistsMeshNodes
  
  !
  !================================================================================================================================
  !

  !>Gets a mesh node number that corresponds to a user node number from mesh nodes. An error will be raised if the user node number does not exist.
  SUBROUTINE MeshTopology_NodeGetMeshNodes(meshNodes,userNodeNumber,meshNodeNumber,err,error,*)

    !Argument variables
    TYPE(MeshNodesType), POINTER :: meshNodes !<A pointer to the mesh nodes to get the node from
    INTEGER(INTG), INTENT(IN) :: userNodeNumber !<The user node number to get
    INTEGER(INTG), INTENT(OUT) :: meshNodeNumber !<On exit, the mesh number corresponding to the user node number. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: meshComponentNumber
    LOGICAL :: nodeExists
    TYPE(MeshType), POINTER :: mesh
    TYPE(MeshComponentTopologyType), POINTER :: meshComponentTopology
    TYPE(VARYING_STRING) :: localError

    ENTERS("MeshTopology_NodeGetMeshNodes",err,error,*999)

    CALL MeshTopology_NodeCheckExistsMeshNodes(meshNodes,userNodeNumber,nodeExists,meshNodeNumber,err,error,*999)
    IF(.NOT.nodeExists) THEN
      meshComponentTopology=>meshNodes%meshComponentTopology
      IF(ASSOCIATED(meshComponentTopology)) THEN
        mesh=>meshComponentTopology%mesh
        IF(ASSOCIATED(mesh)) THEN
          meshComponentNumber=meshComponentTopology%meshComponentNumber
          localError="The user node number "//TRIM(NumberToVString(userNodeNumber,"*",err,error))// &
            & " does not exist in mesh component number "//TRIM(NumberToVString(meshComponentNumber,"*",err,error))// &
            & " of mesh number "//TRIM(NumberToVString(mesh%userNumber,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ELSE
          localError="The user node number "//TRIM(NumberToVString(userNodeNumber,"*",err,error))// &
            & " does not exist in mesh component number "//TRIM(NumberToVString(meshComponentNumber,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        localError="The user node number "//TRIM(NumberToVString(userNodeNumber,"*",err,error))//" does not exist."
        CALL FlagError(localError,err,error,*999)        
      ENDIF
    ENDIF
    
    EXITS("MeshTopology_NodeGetMeshNodes")
    RETURN
999 meshNodeNumber=0
    ERRORSEXITS("MeshTopology_NodeGetMeshNodes",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_NodeGetMeshNodes
  
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
    IF(ASSOCIATED(node%surroundingElements)) DEALLOCATE(node%surroundingElements)
  
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
    NULLIFY(node%surroundingElements)
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
  SUBROUTINE MeshTopology_NodesCalculate(topology,err,error,*)

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: topology !<A pointer to the mesh topology
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,elementIdx,insertStatus,localNodeIdx,globalNode,meshNodeIdx,meshNode,numberOfNodes 
    INTEGER(INTG), POINTER :: globalNodeNumbers(:)
    TYPE(BasisType), POINTER :: basis
    TYPE(MeshType), POINTER :: mesh
    TYPE(MeshElementsType), POINTER :: elements
    TYPE(MeshNodesType), POINTER :: meshNodes
    TYPE(NodesType), POINTER :: nodes
    TYPE(TREE_TYPE), POINTER :: globalNodesTree
    TYPE(TREE_NODE_TYPE), POINTER :: treeNode
    TYPE(VARYING_STRING) :: dummyError,localError

    NULLIFY(globalNodeNumbers)
    NULLIFY(globalNodesTree)
    
    ENTERS("MeshTopology_NodesCalculate",err,error,*998)

    IF(ASSOCIATED(topology)) THEN
      elements=>topology%elements
      IF(ASSOCIATED(elements)) THEN
        meshNodes=>topology%nodes
        IF(ASSOCIATED(meshNodes)) THEN
          mesh=>topology%mesh
          IF(ASSOCIATED(mesh)) THEN
            NULLIFY(nodes)
            CALL MeshGlobalNodesGet(mesh,nodes,err,error,*999)
            IF(ALLOCATED(meshNodes%nodes)) THEN
              localError="Mesh number "//TRIM(NumberToVString(mesh%userNumber,"*",err,error))// &
                & " already has allocated mesh topology nodes."
              CALL FlagError(localError,err,error,*998)
            ELSE
              !Work out what nodes are in the mesh
              CALL TREE_CREATE_START(globalNodesTree,err,error,*999)
              CALL TREE_INSERT_TYPE_SET(globalNodesTree,TREE_NO_DUPLICATES_ALLOWED,err,error,*999)
              CALL TREE_CREATE_FINISH(globalNodesTree,err,error,*999)
              DO elementIdx=1,elements%numberOfElements
                basis=>elements%elements(elementIdx)%basis
                DO localNodeIdx=1,basis%numberOfNodes
                  globalNode=elements%elements(elementIdx)%globalElementNodes(localNodeIdx)
                  CALL TREE_ITEM_INSERT(globalNodesTree,globalNode,globalNode,insertStatus,err,error,*999)
                ENDDO !localNodeIdx
              ENDDO !elementIdx
              CALL TREE_DETACH_AND_DESTROY(globalNodesTree,numberOfNodes,globalNodeNumbers,err,error,*999)
              !Set up the mesh nodes.
              ALLOCATE(meshNodes%nodes(numberOfNodes),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate mesh topology nodes nodes.",err,error,*999)
              CALL TREE_CREATE_START(meshNodes%nodesTree,err,error,*999)
              CALL TREE_INSERT_TYPE_SET(meshNodes%nodesTree,TREE_NO_DUPLICATES_ALLOWED,err,error,*999)
              CALL TREE_CREATE_FINISH(meshNodes%nodesTree,err,error,*999) 
              DO meshNodeIdx=1,numberOfNodes
                CALL MeshTopology_NodeInitialise(meshNodes%nodes(meshNodeIdx),err,error,*999)
                meshNodes%nodes(meshNodeIdx)%meshNumber=meshNodeIdx
                meshNodes%nodes(meshNodeIdx)%globalNumber=globalNodeNumbers(meshNodeIdx)
                meshNodes%nodes(meshNodeIdx)%userNumber=nodes%nodes(globalNodeNumbers(meshNodeIdx))%userNumber
                CALL TREE_ITEM_INSERT(meshNodes%nodesTree,globalNodeNumbers(meshNodeIdx),meshNodeIdx,insertStatus,err,error,*999)
              ENDDO !nodeIdx
              meshNodes%numberOfNodes=numberOfNodes
              IF(ASSOCIATED(globalNodeNumbers)) DEALLOCATE(globalNodeNumbers)
              !Now recalculate the mesh element nodes
              DO elementIdx=1,elements%numberOfElements
                basis=>elements%elements(elementIdx)%basis
                ALLOCATE(elements%elements(elementIdx)%meshElementNodes(basis%numberOfNodes),STAT=err)
                IF(err/=0) CALL FlagError("Could not allocate mesh topology elements mesh element nodes.",err,error,*999)
                DO localNodeIdx=1,basis%numberOfNodes
                  globalNode=elements%elements(elementIdx)%globalElementNodes(localNodeIdx)
                  NULLIFY(treeNode)
                  CALL TREE_SEARCH(meshNodes%nodesTree,globalNode,treeNode,err,error,*999)
                  IF(ASSOCIATED(treeNode)) THEN
                    CALL TREE_NODE_VALUE_GET(meshNodes%nodesTree,treeNode,meshNode,err,error,*999)
                    elements%elements(elementIdx)%meshElementNodes(localNodeIdx)=meshNode
                  ELSE
                    localError="Could not find global node "//TRIM(NumberToVString(globalNode,"*",err,error))//" (user node "// &
                      & TRIM(NumberToVString(nodes%nodes(globalNode)%userNumber,"*",err,error))//") in the mesh nodes."
                    CALL FlagError(localError,err,error,*999)
                  ENDIF
                ENDDO !localNodeIdx
              ENDDO !elementIdx
            ENDIF
          ELSE
            CALL FlagError("Mesh topology mesh is not associated.",err,error,*998)
          ENDIF
        ELSE
          CALL FlagError("Mesh topology nodes is not associated.",err,error,*998)
        ENDIF
      ELSE
        CALL FlagError("Mesh topology elements is not associated.",err,error,*998)
      ENDIF
    ELSE
      CALL FlagError("Mesh topology is not associated.",err,error,*998)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
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
    IF(ASSOCIATED(globalNodesTree)) CALL TREE_DESTROY(globalNodesTree,dummyErr,dummyError,*998)
998 ERRORSEXITS("MeshTopology_NodesCalculate",err,error)
    RETURN 1
   
  END SUBROUTINE MeshTopology_NodesCalculate

  !
  !================================================================================================================================
  !

  !>Destroys the nodes in a mesh topology.
  SUBROUTINE MeshTopology_NodesDestroy(nodes,err,error,*)

    !Argument variables
    TYPE(MeshNodesType), POINTER :: nodes !<A pointer to the mesh nodes to destroy 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshTopology_NodesDestroy",err,error,*999)

    IF(ASSOCIATED(nodes)) THEN
      CALL MeshTopology_NodesFinalise(nodes,err,error,*999)
    ELSE
      CALL FlagError("Mesh topology nodes is not associated",err,error,*999)
    ENDIF

    EXITS("MeshTopology_NodesDestroy")
    RETURN
999 ERRORSEXITS("MeshTopology_NodesDestroy",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_NodesDestroy
  
  !
  !================================================================================================================================
  !

  !>Finalises the given mesh topology node. 
  SUBROUTINE MeshTopology_NodeDerivativeFinalise(nodeDerivative,err,error,*)

    !Argument variables
    TYPE(MeshNodeDerivativeType) :: nodeDerivative !<The mesh node derivative to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshTopology_NodeDerivativeFinalise",err,error,*999)

    IF(ALLOCATED(nodeDerivative%userVersionNumbers)) DEALLOCATE(nodeDerivative%userVersionNumbers)
    IF(ALLOCATED(nodeDerivative%dofIndex)) DEALLOCATE(nodeDerivative%dofIndex)

    EXITS("MeshTopology_NodeDerivativeFinalise")
    RETURN
999 ERRORSEXITS("MeshTopology_NodeDerivativeFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_NodeDerivativeFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the given mesh topology node.
  SUBROUTINE MeshTopology_NodeDerivativeInitialise(nodeDerivative,err,error,*)

    !Argument variables
    TYPE(MeshNodeDerivativeType) :: nodeDerivative !<The mesh node derivative to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshTopology_NodeDerivativeInitialise",err,error,*999)

    nodeDerivative%numberOfVersions=0
    nodeDerivative%globalDerivativeIndex=0
    nodederivative%partialDerivativeIndex=0

    EXITS("MeshTopology_NodeDerivativeInitialise")
    RETURN
999 ERRORSEXITS("MeshTopology_NodeDerivativeInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE MeshTopology_NodeDerivativeInitialise

  !
  !================================================================================================================================
  !

  !>Calculates the number of derivatives at each node in a topology.
  SUBROUTINE MeshTopology_NodesDerivativesCalculate(topology,err,error,*)

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: topology !<A pointer to the mesh topology to calculate the derivates at each node for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: derivativeIdx,element,elementIdx,globalDerivative,localNodeIdx,maxNumberOfDerivatives,nodeIdx, &
      & numberOfDerivatives
    INTEGER(INTG), ALLOCATABLE :: derivatives(:)
    LOGICAL :: found
    TYPE(LIST_TYPE), POINTER :: nodeDerivativeList
    TYPE(MeshElementsType), POINTER :: elements
    TYPE(MeshNodesType), POINTER :: nodes
    TYPE(BasisType), POINTER :: basis
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("MeshTopology_NodesDerivativesCalculate",err,ERROR,*999)

     IF(ASSOCIATED(topology)) THEN
       elements=>topology%elements
       IF(ASSOCIATED(elements)) THEN
         nodes=>topology%nodes
         IF(ASSOCIATED(nodes)) THEN
          !Loop over the mesh nodes
          DO nodeIdx=1,nodes%numberOfNodes
            !Calculate the number of derivatives and versions at each node. This needs to be calculated by looking at the
            !mesh elements as we may have an adjacent element in another domain with a higher order basis also with versions.
            NULLIFY(nodeDerivativeList)
            CALL LIST_CREATE_START(nodeDerivativeList,err,error,*999)
            CALL LIST_DATA_TYPE_SET(nodeDerivativeList,LIST_INTG_TYPE,err,error,*999)
            CALL LIST_INITIAL_SIZE_SET(nodeDerivativeList,8,err,error,*999)
            CALL LIST_CREATE_FINISH(nodeDerivativeList,err,error,*999)
            maxNumberOfDerivatives=-1
            DO elementIdx=1,nodes%nodes(nodeIdx)%numberOfSurroundingElements
              element=nodes%nodes(nodeIdx)%surroundingElements(elementIdx)
              basis=>elements%elements(element)%basis
              !Find the local node corresponding to this node
              found=.FALSE.
              DO localNodeIdx=1,basis%numberOfNodes
                IF(elements%elements(element)%meshElementNodes(localNodeIdx)==nodeIdx) THEN
                  found=.TRUE.
                  EXIT
                ENDIF
              ENDDO !nn
              IF(found) THEN
                DO derivativeIdx=1,basis%numberOfDerivatives(localNodeIdx)
                  CALL LIST_ITEM_ADD(nodeDerivativeList,basis%partialDerivativeIndex(derivativeIdx,localNodeIdx),err,error,*999)
                ENDDO !derivativeIdx
                IF(basis%numberOfDerivatives(localNodeIdx)>maxNumberOfDerivatives) &
                  & maxNumberOfDerivatives=basis%numberOfDerivatives(localNodeidx)
              ELSE
                CALL FlagError("Could not find local node.",err,error,*999)
              ENDIF
            ENDDO !elem_idx
            CALL LIST_REMOVE_DUPLICATES(nodeDerivativeList,err,error,*999)
            CALL LIST_DETACH_AND_DESTROY(nodeDerivativeList,numberOfDerivatives,derivatives,err,error,*999)
            IF(numberOfDerivatives==maxNumberOfDerivatives) THEN
              !Set up the node derivatives.
              ALLOCATE(nodes%nodes(nodeIdx)%derivatives(maxNumberOfDerivatives),STAT=err)
              nodes%nodes(nodeIdx)%numberOfDerivatives=maxNumberOfDerivatives
              DO derivativeIdx=1,numberOfDerivatives
                CALL MeshTopology_NodeDerivativeInitialise(nodes%nodes(nodeIdx)%derivatives(derivativeIdx),err,error,*999)
                nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%partialDerivativeIndex = derivatives(derivativeIdx)
                globalDerivative=PARTIAL_DERIVATIVE_GLOBAL_DERIVATIVE_MAP(derivatives(derivativeIdx))
                IF(globalDerivative/=0) THEN
                  nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%globalDerivativeIndex=globalDerivative
                ELSE
                  localError="The partial derivative index of "//TRIM(NumberToVstring(derivatives(derivativeIdx),"*", &
                    & err,error))//" for derivative number "//TRIM(NumberToVstring(derivativeIdx,"*",err,error))// &
                    & " does not have a corresponding global derivative."
                  CALL FlagError(localError,err,error,*999)
                ENDIF
              ENDDO !derivativeIdx
              DEALLOCATE(derivatives)
            ELSE
              localError="Invalid mesh configuration. User node "// &
                & TRIM(NumberToVstring(nodes%nodes(nodeIdx)%userNumber,"*",err,error))// &
                & " has inconsistent derivative directions."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ENDDO !nodeIdx
        ELSE
          CALL FlagError("Mesh topology nodes is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Mesh topology elements is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh topology is not associated.",err,error,*999)
    ENDIF
   
    IF(diagnostics1) THEN
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Number of mesh global nodes = ",nodes%numberOfNodes,err,error,*999)
      DO nodeIdx=1,nodes%numberOfNodes
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Mesh global node number = ",nodeIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of derivatives = ",nodes%nodes(nodeIdx)%numberOfDerivatives, &
          & err,error,*999)
        DO derivativeIdx=1,nodes%nodes(nodeIdx)%numberOfDerivatives
          !TODO: change output string below so that it writes out derivativeIdx index as well
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Global derivative index(derivativeIdx) = ", &
            & nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%globalDerivativeIndex,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Partial derivative index(derivativeIdx) = ", &
            & nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%partialDerivativeIndex,err,error,*999)
        ENDDO !derivativeIdx
      ENDDO !node_idx
    ENDIF
    
    EXITS("MeshTopology_NodesDerivativesCalculate")
    RETURN
999 IF(ALLOCATED(derivatives)) DEALLOCATE(derivatives)
    IF(ASSOCIATED(nodeDerivativeList)) CALL LIST_DESTROY(nodeDerivativeList,err,error,*998)
998 ERRORSEXITS("MeshTopology_NodesDerivativesCalculate",err,error)
    RETURN 1
   
  END SUBROUTINE MeshTopology_NodesDerivativesCalculate

  !
  !================================================================================================================================
  !
  
  !>Returns the number of derivatives for a node in a mesh
  SUBROUTINE MeshTopology_NodeNumberOfDerivativesGet(meshNodes,userNumber,numberOfDerivatives,err,error,*)

    !Argument variables
    TYPE(MeshNodesType), POINTER :: meshNodes !<A pointer to the mesh nodes containing the nodes to get the number of derivatives for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user node number to get the number of derivatives for
    INTEGER(INTG), INTENT(OUT) :: numberOfDerivatives !<On return, the number of global derivatives at the specified user node number.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: meshComponentNumber,meshNumber
    LOGICAL :: nodeExists
    TYPE(MeshType), POINTER :: mesh
    TYPE(MeshComponentTopologyType), POINTER :: meshComponentTopology
    TYPE(VARYING_STRING) :: localError

    ENTERS("MeshTopology_NodeNumberOfDerivativesGet",err,error,*999)

    IF(ASSOCIATED(meshNodes)) THEN
      CALL MeshTopology_NodeCheckExists(meshNodes,userNumber,nodeExists,meshNumber,err,error,*999)
      IF(nodeExists) THEN
        numberOfDerivatives=meshNodes%nodes(meshNumber)%numberOfDerivatives
      ELSE
        meshComponentTopology=>meshNodes%meshComponentTopology
        IF(ASSOCIATED(meshComponentTopology)) THEN
          mesh=>meshComponentTopology%mesh
          IF(ASSOCIATED(mesh)) THEN
            meshComponentNumber=meshComponentTopology%meshComponentNumber
            localError="The user node number "//TRIM(NumberToVString(userNumber,"*",err,error))// &
              & " does not exist in mesh component number "//TRIM(NumberToVString(meshComponentNumber,"*",err,error))// &
              & " of mesh number "//TRIM(NumberToVString(mesh%userNumber,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ELSE
            CALL FlagError("Mesh component topology mesh is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Mesh nodes mesh component topology is not associated.",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Mesh nodes is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopology_NodeNumberOfDerivativesGet")
    RETURN
999 ERRORSEXITS("MeshTopology_NodeNumberOfDerivativesGet",err,error)    
    RETURN 1
   
  END SUBROUTINE MeshTopology_NodeNumberOfDerivativesGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the global derivative numbers for a node in mesh nodes
  SUBROUTINE MeshTopology_NodeDerivativesGet(meshNodes,userNumber,derivatives,err,error,*)

    !Argument variables
    TYPE(MeshNodesType), POINTER :: meshNodes !<A pointer to the mesh nodes containing the node to get the derivatives for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the node to get the derivatives for
    INTEGER(INTG), INTENT(OUT) :: derivatives(:) !<On return, the global derivatives at the specified node
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: derivativeIdx,meshComponentNumber,meshNumber,numberOfDerivatives
    LOGICAL :: nodeExists
    TYPE(MeshType), POINTER :: mesh
    TYPE(MeshComponentTopologyType), POINTER :: meshComponentTopology
    TYPE(VARYING_STRING) :: localError

    ENTERS("MeshTopology_NodeDerivativesGet",err,error,*999)

    IF(ASSOCIATED(meshNodes)) THEN
      CALL MeshTopology_NodeCheckExists(meshNodes,userNumber,nodeExists,meshNumber,err,error,*999)
      IF(nodeExists) THEN         
        numberOfDerivatives=meshNodes%nodes(meshNumber)%numberOfDerivatives
        IF(SIZE(derivatives,1)>=numberOfDerivatives) THEN
          DO derivativeIdx=1,numberOfDerivatives
            derivatives(derivativeIdx)=meshNodes%nodes(meshNumber)%derivatives(derivativeIdx)%globalDerivativeIndex
          ENDDO !derivativeIdx
        ELSE
          localError="The size of the supplied derivatives array of "// &
            & TRIM(NumberToVString(SIZE(derivatives,1),"*",err,error))// &
            & " is too small. The size should be >= "// &
            & TRIM(NumberToVString(numberOfDerivatives,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        meshComponentTopology=>meshNodes%meshComponentTopology
        IF(ASSOCIATED(meshComponentTopology)) THEN
          mesh=>meshComponentTopology%mesh
          IF(ASSOCIATED(mesh)) THEN
            meshComponentNumber=meshComponentTopology%meshComponentNumber
            localError="The user node number "//TRIM(NumberToVString(userNumber,"*",err,error))// &
              & " does not exist in mesh component number "//TRIM(NumberToVString(meshComponentNumber,"*",err,error))// &
              & " of mesh number "//TRIM(NumberToVString(mesh%userNumber,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ELSE
            CALL FlagError("Mesh component topology mesh is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Mesh nodes mesh component topology is not associated.",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Mesh nodes is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopology_NodeDerivativesGet")
    RETURN
999 ERRORSEXITS("MeshTopology_NodeDerivativesGet",err,error)    
    RETURN 1
   
  END SUBROUTINE MeshTopology_NodeDerivativesGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of versions for a derivative of a node in mesh nodes
  SUBROUTINE MeshTopology_NodeNumberOfVersionsGet(meshNodes,derivativeNumber,userNumber,numberOfVersions,err,error,*)

    !Argument variables
    TYPE(MeshNodesType), POINTER :: meshNodes !<A pointer to the mesh nodes containing the node to get the number of versions for
    INTEGER(INTG), INTENT(IN) :: derivativeNumber !<The derivative number of the node to get the number of versions for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the node to get the number of versions for
    INTEGER(INTG), INTENT(OUT) :: numberOfVersions !<On return, the number of versions for the specified derivative of the specified node
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: meshComponentNumber,meshNumber
    LOGICAL :: nodeExists
    TYPE(MeshType), POINTER :: mesh
    TYPE(MeshComponentTopologyType), POINTER :: meshComponentTopology
    TYPE(VARYING_STRING) :: localError

    ENTERS("MeshTopology_NodeNumberOfVersionsGet",err,error,*999)

    IF(ASSOCIATED(meshNodes)) THEN
      CALL MeshTopology_NodeCheckExists(meshNodes,userNumber,nodeExists,meshNumber,err,error,*999)
      IF(nodeExists) THEN
        IF(derivativeNumber>=1.AND.derivativeNumber<=meshNodes%nodes(meshNumber)%numberOfDerivatives) THEN
          numberOfVersions=meshNodes%nodes(meshNumber)%derivatives(derivativeNumber)%numberOfVersions
        ELSE
          localError="The specified derivative index of "// &
            & TRIM(NumberToVString(derivativeNumber,"*",err,error))// &
            & " is invalid. The derivative index must be >= 1 and <= "// &
            & TRIM(NumberToVString(meshNodes%nodes(meshNumber)%numberOfDerivatives,"*",err,error))// &
            & " for user node number "//TRIM(NumberToVString(userNumber,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        meshComponentTopology=>meshNodes%meshComponentTopology
        IF(ASSOCIATED(meshComponentTopology)) THEN
          mesh=>meshComponentTopology%mesh
          IF(ASSOCIATED(mesh)) THEN
            meshComponentNumber=meshComponentTopology%meshComponentNumber
            localError="The user node number "//TRIM(NumberToVString(userNumber,"*",err,error))// &
              & " does not exist in mesh component number "//TRIM(NumberToVString(meshComponentNumber,"*",err,error))// &
              & " of mesh number "//TRIM(NumberToVString(mesh%userNumber,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ELSE
            CALL FlagError("Mesh component topology mesh is not associated.",err,error,*999)
          ENDIF
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Mesh nodes is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopology_NodeNumberOfVersionsGet")
    RETURN
999 ERRORSEXITS("MeshTopology_NodeNumberOfVersionsGet",err,error)    
    RETURN 1
   
  END SUBROUTINE MeshTopology_NodeNumberOfVersionsGet
  
  !
  !================================================================================================================================
  !

  !>Returns the number of nodes for a node in a mesh
  SUBROUTINE MeshTopology_NodesNumberOfNodesGet(meshNodes,numberOfNodes,err,error,*)

    !Argument variables
    TYPE(MeshNodesType), POINTER :: meshNodes !<A pointer to the mesh nodes containing the nodes to get the number of nodes for
    INTEGER(INTG), INTENT(OUT) :: numberOfNodes !<On return, the number of nodes in the mesh.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
  
    ENTERS("MeshTopology_NodesNumberOfNodesGet",err,error,*999)

    IF(ASSOCIATED(meshNodes)) THEN
      numberOfNodes=meshNodes%numberOfNodes
    ELSE
      CALL FlagError("Mesh nodes is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopology_NodesNumberOfNodesGet")
    RETURN
999 ERRORSEXITS("MeshTopology_NodesNumberOfNodesGet",err,error)    
    RETURN 1
   
  END SUBROUTINE MeshTopology_NodesNumberOfNodesGet
  
  !
  !================================================================================================================================
  !

  !>Calculates the number of versions at each node in a topology.
  SUBROUTINE MeshTopology_NodesVersionCalculate(topology,err,error,*)

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: topology !<A pointer to the mesh topology to calculate the versions at each node for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: element,localNodeIdx,derivativeIdx,nodeIdx,numberOfVersions,versionIdx
    INTEGER(INTG), ALLOCATABLE :: versions(:)
    TYPE(LIST_PTR_TYPE), POINTER :: nodeVersionList(:,:)
    TYPE(MeshElementsType), POINTER :: elements
    TYPE(MeshNodesType), POINTER :: nodes
    TYPE(BasisType), POINTER :: basis
    
    ENTERS("MeshTopology_NodesVersionCalculate",err,error,*999)

    IF(ASSOCIATED(topology)) THEN
      elements=>topology%elements
      IF(ASSOCIATED(elements)) THEN
        nodes=>topology%nodes
        IF(ASSOCIATED(nodes)) THEN
          !Loop over the mesh elements
          !Calculate the number of versions at each node. This needs to be calculated by looking at all the mesh elements
          !as we may have an adjacent elements in another domain with a higher order basis along with different versions
          !being assigned to its derivatives.
          !\todo : See if there are any constraints that can be applied to restrict the amount of memory being allocated here
          ALLOCATE(nodeVersionList(MAXIMUM_GLOBAL_DERIV_NUMBER,nodes%numberOfNodes),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate node version list.",err,error,*999)
          DO nodeIdx=1,nodes%numberOfNodes
            DO derivativeIdx=1,nodes%nodes(nodeIdx)%numberOfDerivatives
              NULLIFY(nodeVersionList(derivativeIdx,nodeIdx)%ptr)
              CALL LIST_CREATE_START(nodeVersionList(derivativeIdx,nodeIdx)%ptr,err,error,*999)
              CALL LIST_DATA_TYPE_SET(nodeVersionList(derivativeIdx,nodeIdx)%ptr,LIST_INTG_TYPE,err,error,*999)
              CALL LIST_INITIAL_SIZE_SET(nodeVersionList(derivativeIdx,nodeIdx)%ptr,8,err,error,*999)
              CALL LIST_CREATE_FINISH(nodeVersionList(derivativeIdx,nodeIdx)%ptr,err,error,*999)
            ENDDO!derivativeIdx
          ENDDO!nodeIdx
          DO element=1,elements%numberOfElements
            basis=>elements%elements(element)%basis
            DO localNodeIdx=1,BASIS%numberOfNodes
              DO derivativeIdx=1,basis%numberOfDerivatives(localNodeIdx)
                CALL LIST_ITEM_ADD(nodeVersionList(derivativeIdx,elements%elements(element)% &
                  & meshElementNodes(localNodeIdx))%ptr,elements%elements(element)%userElementNodeVersions( &
                  & derivativeIdx,localNodeIdx),err,error,*999)
              ENDDO!derivativeIdx
            ENDDO!localNodeIdx
          ENDDO!element
          DO nodeIdx=1,nodes%numberOfNodes
            DO derivativeIdx=1,nodes%nodes(nodeIdx)%numberOfDerivatives
              CALL LIST_REMOVE_DUPLICATES(nodeVersionList(derivativeIdx,nodeIdx)%ptr,err,error,*999)
              CALL LIST_DETACH_AND_DESTROY(nodeVersionList(derivativeIdx,nodeIdx)%ptr,numberOfVersions,versions, &
                & err,error,*999)
              nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%numberOfVersions = MAXVAL(versions(1:numberOfVersions))
              ALLOCATE(nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%userVersionNumbers(nodes%nodes(nodeIdx)% &
                & derivatives(derivativeIdx)%numberOfVersions),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate node global derivative index.",err,error,*999)
              DO versionIdx=1,nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%numberOfVersions 
                nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%userVersionNumbers(versionIdx) = versionIdx
              ENDDO !versionIdx
              DEALLOCATE(versions)
            ENDDO!derivativeIdx
          ENDDO!nodeIdx
          DEALLOCATE(nodeVersionList)
          NULLIFY(nodeVersionList)
        ELSE
          CALL FlagError("Mesh topology nodes is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Mesh topology elements is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh topology is not associated.",err,error,*999)
    ENDIF
   
    IF(diagnostics1) THEN
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Number of mesh global nodes = ",nodes%numberOfNodes,err,error,*999)
      DO nodeIdx=1,nodes%numberOfNodes
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Mesh global node number = ",nodeIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of derivatives = ", &
          & nodes%nodes(nodeIdx)%numberOfDerivatives,err,error,*999)
        DO derivativeIdx=1,nodes%nodes(nodeIdx)%numberOfDerivatives
          !\todo : change output string below so that it writes out derivativeIdx index as well
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Global derivative index(derivativeIdx) = ", &
            & nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%globalDerivativeIndex,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Partial derivative index(derivativeIdx) = ", &
            & nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%partialDerivativeIndex,err,error,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1, &
            & nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%numberOfVersions,8,8, &
            & nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%userVersionNumbers, & 
            & '("    User Version index(derivativeIdx,:) :",8(X,I2))','(36X,8(X,I2))',err,error,*999)
        ENDDO!derivativeIdx
      ENDDO !nodeIdx
    ENDIF
    
    EXITS("MeshTopology_NodesVersionCalculate")
    RETURN
999 IF(ALLOCATED(versions)) DEALLOCATE(versions)
    IF(ASSOCIATED(nodeVersionList)) THEN
      DO nodeIdx=1,nodes%numberOfNodes
        DO derivativeIdx=1,nodes%nodes(nodeIdx)%numberOfDerivatives
          CALL LIST_DESTROY(nodeVersionList(derivativeIdx,nodeIdx)%ptr,err,error,*998)
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
  
  !>Returns if the node in a mesh is on the boundary or not
  SUBROUTINE MeshTopology_NodeOnBoundaryGet(meshNodes,userNumber,onBoundary,err,error,*)

    !Argument variables
    TYPE(MeshNodesType), POINTER :: meshNodes !<A pointer to the mesh nodes containing the nodes to get the boundary type for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user node number to get the boundary type for
    INTEGER(INTG), INTENT(OUT) :: onBoundary !<On return, the boundary type of the specified user node number.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: meshNumber
 
    ENTERS("MeshTopology_NodeOnBoundaryGet",err,error,*999)

    IF(.NOT.ASSOCIATED(meshNodes)) CALL FlagError("Mesh nodes is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(meshNodes%nodes)) CALL FlagError("Mesh nodes nodes is not allocated.",err,error,*999)
    
    CALL MeshTopology_NodeGet(meshNodes,userNumber,meshNumber,err,error,*999)
    IF(meshNodes%nodes(meshNumber)%boundaryNode) THEN
      onBoundary=MESH_ON_DOMAIN_BOUNDARY
    ELSE
      onBoundary=MESH_OFF_DOMAIN_BOUNDARY
    ENDIF
    
    EXITS("MeshTopology_NodeOnBoundaryGet")
    RETURN
999 ERRORSEXITS("MeshTopology_NodeOnBoundaryGet",err,error)    
    RETURN 1
   
  END SUBROUTINE MeshTopology_NodeOnBoundaryGet
  
  !
  !================================================================================================================================
  !

  !>Calculates the element numbers surrounding a node for a mesh.
  SUBROUTINE MeshTopology_SurroundingElementsCalculate(topology,err,error,*)

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: topology !<A pointer to the mesh topology to calculate the elements surrounding each node for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: element,elementIdx,insertPosition,localNodeIdx,node,surroundingElementNumber
    INTEGER(INTG), POINTER :: newSurroundingElements(:)
    LOGICAL :: foundElement
    TYPE(BasisType), POINTER :: basis
    TYPE(MeshElementsType), POINTER :: elements
    TYPE(MeshNodesType), POINTER :: nodes

    NULLIFY(newSurroundingElements)

    ENTERS("MeshTopology_SurroundingElementsCalculate",err,error,*999)
    
    IF(ASSOCIATED(topology)) THEN
      elements=>topology%elements      
      IF(ASSOCIATED(elements)) THEN
        nodes=>topology%nodes       
        IF(ASSOCIATED(nodes)) THEN
          IF(ALLOCATED(nodes%nodes)) THEN
            DO elementIdx=1,elements%numberOfElements
              basis=>elements%elements(elementIdx)%basis
              DO localNodeIdx=1,basis%numberOfNodes
                node=elements%elements(elementIdx)%meshElementNodes(localNodeIdx)
                foundElement=.FALSE.
                element=1
                insertPosition=1
                DO WHILE(element<=nodes%nodes(node)%numberOfSurroundingElements.AND..NOT.foundElement)
                  surroundingElementNumber=nodes%nodes(node)%surroundingElements(element)
                  IF(surroundingElementNumber==elementIdx) THEN
                    foundElement=.TRUE.
                  ENDIF
                  element=element+1
                  IF(elementIdx>=surroundingElementNumber) THEN
                    insertPosition=element
                  ENDIF
                ENDDO
                IF(.NOT.foundElement) THEN
                  !Insert element into surrounding elements
                  ALLOCATE(newSurroundingElements(nodes%nodes(node)%numberOfSurroundingElements+1),STAT=err)
                  IF(err/=0) CALL FlagError("Could not allocate new surrounding elements.",err,error,*999)
                  IF(ASSOCIATED(nodes%nodes(node)%surroundingElements)) THEN
                    newSurroundingElements(1:insertPosition-1)=nodes%nodes(node)%surroundingElements(1:insertPosition-1)
                    newSurroundingElements(insertPosition)=elementIdx
                    newSurroundingElements(insertPosition+1:nodes%nodes(node)%numberOfSurroundingElements+1)= &
                      & nodes%nodes(node)%surroundingElements(insertPosition:nodes%nodes(node)%numberOfSurroundingElements)
                    DEALLOCATE(nodes%nodes(node)%surroundingElements)
                  ELSE
                    newSurroundingElements(1)=elementIdx
                  ENDIF
                  nodes%nodes(node)%surroundingElements=>newSurroundingElements
                  nodes%nodes(node)%numberOfSurroundingElements=nodes%nodes(node)%numberOfSurroundingElements+1
                ENDIF
              ENDDO !localNodeIdx
            ENDDO !elementIdx
          ELSE
            CALL FlagError("Mesh topology nodes nodes have not been allocated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Mesh topology nodes are not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Mesh topology elements is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Mesh topology not associated.",err,error,*999)
    ENDIF

    EXITS("MeshTopology_SurroundingElementsCalculate")
    RETURN
999 IF(ASSOCIATED(newSurroundingElements)) DEALLOCATE(newSurroundingElements)
    ERRORSEXITS("MeshTopology_SurroundingElementsCalculate",err,error)
    RETURN 1   
  END SUBROUTINE MeshTopology_SurroundingElementsCalculate
  
  !
  !===============================================================================================================================
  !

  !>Finalises the nodes data structures for a mesh topology and deallocates any memory.
  SUBROUTINE MeshTopology_NodesFinalise(nodes,err,error,*)

    !Argument variables
    TYPE(MeshNodesType), POINTER :: nodes !<A pointer to the mesh topology nodes to finalise the nodes for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: nodeIdx

    ENTERS("MeshTopology_NodesFinalise",err,error,*999)

    IF(ASSOCIATED(nodes)) THEN
      IF(ALLOCATED(nodes%nodes)) THEN
        DO nodeIdx=1,SIZE(nodes%nodes,1)
          CALL MeshTopology_NodeFinalise(nodes%nodes(nodeIdx),err,error,*999)
        ENDDO !nodesIdx
        DEALLOCATE(nodes%nodes)
      ENDIF
      IF(ASSOCIATED(nodes%nodesTree)) CALL TREE_DESTROY(nodes%nodesTree,err,error,*999)
      DEALLOCATE(nodes)
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
  SUBROUTINE MeshTopology_NodesInitialise(topology,err,error,*)

    !Argument variables
    TYPE(MeshComponentTopologyType), POINTER :: topology !<A pointer to the mesh topology to initialise the nodes for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MeshTopology_NodesInitialise",err,error,*999)

    IF(ASSOCIATED(topology)) THEN
      IF(ASSOCIATED(topology%nodes)) THEN
        CALL FlagError("Mesh already has topology nodes associated.",err,error,*999)
      ELSE
        ALLOCATE(topology%nodes,STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate topology nodes.",err,error,*999)
        topology%nodes%numberOfNodes=0
        topology%nodes%meshComponentTopology=>topology
        NULLIFY(topology%nodes%nodesTree)
      ENDIF
    ELSE
      CALL FlagError("Topology is not associated.",err,error,*999)
    ENDIF
    
    EXITS("MeshTopology_NodesInitialise")
    RETURN
999 ERRORSEXITS("MeshTopology_NodesInitialise",err,error)
    RETURN 1
  END SUBROUTINE MeshTopology_NodesInitialise

  !
  !================================================================================================================================
  !

  !>Finalises the meshes and deallocates all memory
  SUBROUTINE MESHES_FINALISE(MESHES,ERR,ERROR,*)

   !Argument variables
    TYPE(MeshesType), POINTER :: MESHES !<A pointer to the meshes to finalise the meshes for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(MeshType), POINTER :: MESH
 
    ENTERS("MESHES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(MESHES)) THEN
      DO WHILE(MESHES%numberOfMeshes>0)
        MESH=>MESHES%MESHES(1)%PTR
        CALL MESH_DESTROY(MESH,ERR,ERROR,*999)
      ENDDO !mesh_idx
      DEALLOCATE(MESHES)
    ELSE
      CALL FlagError("Meshes is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("MESHES_FINALISE")
    RETURN
999 ERRORSEXITS("MESHES_FINALISE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE MESHES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the generic meshes.
  SUBROUTINE MESHES_INITIALISE_GENERIC(MESHES,ERR,ERROR,*)

    !Argument variables
    TYPE(MeshesType), POINTER :: MESHES !<A pointer to the meshes to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    ENTERS("MESHES_INITIALISE_GENERIC",ERR,ERROR,*998)

    IF(ASSOCIATED(MESHES)) THEN
      CALL FlagError("Meshes is already associated.",ERR,ERROR,*998)
    ELSE
      ALLOCATE(MESHES,STAT=ERR)
      IF(ERR/=0) CALL FlagError("Meshes could not be allocated",ERR,ERROR,*999)
      NULLIFY(MESHES%REGION)
      NULLIFY(MESHES%INTERFACE)
      MESHES%numberOfMeshes=0
      NULLIFY(MESHES%MESHES)
    ENDIF
    
    EXITS("MESHES_INITIALISE_GENERIC")
    RETURN
999 CALL MESHES_FINALISE(MESHES,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("MESHES_INITIALISE_GENERIC",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESHES_INITIALISE_GENERIC

  !
  !================================================================================================================================
  !

  !>Initialises the meshes for the given interface.
  SUBROUTINE MESHES_INITIALISE_INTERFACE(INTERFACE,ERR,ERROR,*)

    !Argument variables
    TYPE(InterfaceType), POINTER :: INTERFACE !<A pointer to the interface to initialise the meshes for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESHES_INITIALISE_INTERFACE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(ASSOCIATED(INTERFACE%MESHES)) THEN
        LOCAL_ERROR="Interface number "//TRIM(NumberToVString(INTERFACE%userNumber,"*",ERR,ERROR))// &
          & " already has a mesh associated"
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ELSE
        CALL MESHES_INITIALISE_GENERIC(INTERFACE%MESHES,ERR,ERROR,*999)
        INTERFACE%MESHES%INTERFACE=>INTERFACE
      ENDIF
    ELSE
      CALL FlagError("Interface is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESHES_INITIALISE_INTERFACE")
    RETURN
999 ERRORSEXITS("MESHES_INITIALISE_INTERFACE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESHES_INITIALISE_INTERFACE

  !
  !================================================================================================================================
  !

  !>Initialises the meshes for the given region.
  SUBROUTINE MESHES_INITIALISE_REGION(REGION,ERR,ERROR,*)

    !Argument variables
    TYPE(RegionType), POINTER :: REGION !<A pointer to the region to initialise the meshes for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MESHES_INITIALISE_REGION",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%MESHES)) THEN
        LOCAL_ERROR="Region number "//TRIM(NumberToVString(REGION%userNumber,"*",ERR,ERROR))// &
          & " already has a mesh associated"
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ELSE
        CALL MESHES_INITIALISE_GENERIC(REGION%MESHES,ERR,ERROR,*999)
        REGION%MESHES%REGION=>REGION
      ENDIF
    ELSE
      CALL FlagError("Region is not associated",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MESHES_INITIALISE_REGION")
    RETURN
999 ERRORSEXITS("MESHES_INITIALISE_REGION",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESHES_INITIALISE_REGION

  !
  !================================================================================================================================
  !

  !>Initialises the embedded meshes.
  SUBROUTINE EMBEDDED_MESH_INITIALISE(MESH_EMBEDDING,ERR,ERROR,*)

    !Argument variables
    !TYPE(MESH_EMBEDDING_TYPE), INTENT(INOUT) :: MESH_EMBEDDING !<Mesh embedding to initialise
    TYPE(MESH_EMBEDDING_TYPE), POINTER :: MESH_EMBEDDING !<Mesh embedding to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("EMBEDDED_MESH_INITIALISE",ERR,ERROR,*998)
    
    ALLOCATE(MESH_EMBEDDING,STAT=ERR)
    NULLIFY(MESH_EMBEDDING%PARENT_MESH)
    NULLIFY(MESH_EMBEDDING%CHILD_MESH)
    
    EXITS("EMBEDDED_MESH_INITIALISE")
    RETURN
!999 CALL EMBEDDED_MESH_FINALISE(MESH_EMBEDDING,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("EMBEDDED_MESH_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE EMBEDDED_MESH_INITIALISE

  !
  !================================================================================================================================
  !

  !>Creates an embedding of one mesh in another
  SUBROUTINE MESH_EMBEDDING_CREATE(MESH_EMBEDDING, PARENT_MESH, CHILD_MESH,ERR,ERROR,*)
!    TYPE(MESH_EMBEDDING_TYPE), INTENT(INOUT) :: MESH_EMBEDDING !<Mesh embedding to create.
    TYPE(MESH_EMBEDDING_TYPE), POINTER :: MESH_EMBEDDING !<Mesh embedding to create.
    TYPE(MeshType), POINTER, INTENT(IN) :: PARENT_MESH !<The parent mesh
    TYPE(MeshType), POINTER, INTENT(IN) :: CHILD_MESH  !<The child mesh
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: NGP = 0, ne

    ENTERS("MESH_EMBEDDING_CREATE",ERR,ERROR,*999) 
    
    WRITE(*,*) 'parent mesh', PARENT_MESH%numberOfElements
    WRITE(*,*) 'child mesh', child_MESH%numberOfElements
    CALL EMBEDDED_MESH_INITIALISE(MESH_EMBEDDING,ERR,ERROR,*999)

    DO ne=1,PARENT_MESH%numberOfElements
      NGP = MAX(NGP,PARENT_MESH%TOPOLOGY(1)%PTR%ELEMENTS%ELEMENTS(ne)%BASIS%QUADRATURE%&
        & QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR%NUMBER_OF_GAUSS)
    ENDDO !ne

    MESH_EMBEDDING%PARENT_MESH => PARENT_MESH
    MESH_EMBEDDING%CHILD_MESH  => CHILD_MESH
    ALLOCATE(MESH_EMBEDDING%CHILD_NODE_XI_POSITION(PARENT_MESH%numberOfElements),STAT=ERR)
    IF(ERR/=0) CALL FlagError("Could not allocate child node positions.",ERR,ERROR,*999)
    ALLOCATE(MESH_EMBEDDING%GAUSS_POINT_XI_POSITION(NGP,PARENT_MESH%numberOfElements),STAT=ERR)
    IF(ERR/=0) CALL FlagError("Could not allocate gauss point positions.",ERR,ERROR,*999)
    
    EXITS("MESH_EMBEDDING_CREATE")
    RETURN 

999 ERRORSEXITS("MESH_EMBEDDING_CREATE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESH_EMBEDDING_CREATE

  !
  !================================================================================================================================
  !

  !>Sets the positions of nodes in the child mesh for one element in the parent mesh
  SUBROUTINE MESH_EMBEDDING_SET_CHILD_NODE_POSITION(MESH_EMBEDDING, ELEMENT_NUMBER, NODE_NUMBERS, XI_COORDS,ERR,ERROR,*)
    TYPE(MESH_EMBEDDING_TYPE), INTENT(INOUT) :: MESH_EMBEDDING !<The mesh embedding object
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER  !<Element number in the parent mesh
    INTEGER(INTG), INTENT(IN) :: NODE_NUMBERS(:) !<NODE_NUMBERS(node_idx) Node numbers in child mesh for the node_idx'th embedded node in the ELEMENT_NUMBER'th element of the parent mesh
    REAL(DP), INTENT(IN)      :: XI_COORDS(:,:)  !<XI_COORDS(:,node_idx) Xi coordinates of the node_idx'th embedded node in the ELEMENT_NUMBER'th

    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    ENTERS("MESH_EMBEDDING_SET_CHILD_NODE_POSITION",ERR,ERROR,*999)

    IF(ELEMENT_NUMBER<1 .OR. ELEMENT_NUMBER > MESH_EMBEDDING%PARENT_MESH%numberOfElements) THEN
      CALL FlagError("Element number out of range",ERR,ERROR,*999)
    ENDIF

    MESH_EMBEDDING%CHILD_NODE_XI_POSITION(ELEMENT_NUMBER)%NUMBER_OF_NODES = SIZE(NODE_NUMBERS)

    ALLOCATE(MESH_EMBEDDING%CHILD_NODE_XI_POSITION(ELEMENT_NUMBER)%NODE_NUMBERS(SIZE(NODE_NUMBERS)))
    MESH_EMBEDDING%CHILD_NODE_XI_POSITION(ELEMENT_NUMBER)%NODE_NUMBERS(1:SIZE(NODE_NUMBERS)) = NODE_NUMBERS(1:SIZE(NODE_NUMBERS))

    ALLOCATE(MESH_EMBEDDING%CHILD_NODE_XI_POSITION(ELEMENT_NUMBER)%XI_COORDS(SIZE(XI_COORDS,1),SIZE(XI_COORDS,2)))
    MESH_EMBEDDING%CHILD_NODE_XI_POSITION(ELEMENT_NUMBER)%XI_COORDS(1:SIZE(XI_COORDS,1),1:SIZE(XI_COORDS,2)) = &
      & XI_COORDS(1:SIZE(XI_COORDS,1),1:SIZE(XI_COORDS,2))
    
    EXITS("MESH_EMBEDDING_SET_CHILD_NODE_POSITION")
    RETURN
999 ERRORSEXITS("MESH_EMBEDDING_SET_CHILD_NODE_POSITION",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESH_EMBEDDING_SET_CHILD_NODE_POSITION

  !
  !================================================================================================================================
  !

  !>Sets the positions of a Gauss point of the parent mesh in terms of element/xi coordinate in the child mesh
  SUBROUTINE MESH_EMBEDDING_SET_GAUSS_POINT_DATA(MESH_EMBEDDING, PARENT_ELEMENT_NUMBER, GAUSSPT_NUMBER,&
    & PARENT_XI_COORD, CHILD_ELEMENT_NUMBER, CHILD_XI_COORD,ERR,ERROR,*)
    TYPE(MESH_EMBEDDING_TYPE), INTENT(INOUT) :: MESH_EMBEDDING   !<The mesh embedding object
    INTEGER(INTG), INTENT(IN) :: PARENT_ELEMENT_NUMBER           !<Element number in the parent mesh
    INTEGER(INTG), INTENT(IN) :: GAUSSPT_NUMBER                  !<Gauss point number in this element
    REAL(DP), INTENT(IN) :: PARENT_XI_COORD(:)              !<Xi coordinate in parent element

    INTEGER(INTG), INTENT(IN) :: CHILD_ELEMENT_NUMBER !<Element number in the child mesh
    REAL(DP), INTENT(IN) :: CHILD_XI_COORD(:)    !<Xi coordinate in child element

    INTEGER(INTG), INTENT(OUT) :: ERR           !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR  !<The error string

    ENTERS("MESH_EMBEDDING_SET_GAUSS_POINT_DATA",ERR,ERROR,*999)

    IF(PARENT_ELEMENT_NUMBER<1 .OR. PARENT_ELEMENT_NUMBER > MESH_EMBEDDING%PARENT_MESH%numberOfElements) THEN
      CALL FlagError("Parent element number out of range",ERR,ERROR,*999)
    ENDIF
    IF(CHILD_ELEMENT_NUMBER<1 .OR. CHILD_ELEMENT_NUMBER > MESH_EMBEDDING%CHILD_MESH%numberOfElements) THEN
      CALL FlagError("Child element number out of range",ERR,ERROR,*999)
    ENDIF
    IF(GAUSSPT_NUMBER<1 .OR. GAUSSPT_NUMBER > SIZE(MESH_EMBEDDING%GAUSS_POINT_XI_POSITION,1)) THEN
      CALL FlagError("Gauss point number out of range",ERR,ERROR,*999)
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
999 ERRORSEXITS("MESH_EMBEDDING_SET_GAUSS_POINT_DATA",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MESH_EMBEDDING_SET_GAUSS_POINT_DATA


  !
  !================================================================================================================================
  !

END MODULE MeshRoutines


