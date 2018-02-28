!> \file
!> \author Chris Bradley
!> \brief This module contains all interface routines.
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
!> Contributor(s): David Nordsletten, Thiranja Prasad Babarenda Gamage
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delte
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!> This module contains all interface routines.
MODULE INTERFACE_ROUTINES

  USE BaseRoutines
  USE BasisRoutines
  USE DataPointRoutines
  USE DataProjectionRoutines
  USE FIELD_ROUTINES
  USE GENERATED_MESH_ROUTINES
  USE INPUT_OUTPUT
  USE InterfaceAccessRoutines
  USE INTERFACE_CONDITIONS_ROUTINES
  USE ISO_VARYING_STRING
  USE Kinds
  USE Lists
  USE MESH_ROUTINES
  USE NODE_ROUTINES
  USE Strings
  USE Types

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  INTERFACE INTERFACE_LABEL_GET
    MODULE PROCEDURE INTERFACE_LABEL_GET_C
    MODULE PROCEDURE INTERFACE_LABEL_GET_VS
  END INTERFACE INTERFACE_LABEL_GET
  
  INTERFACE INTERFACE_LABEL_SET
    MODULE PROCEDURE INTERFACE_LABEL_SET_C
    MODULE PROCEDURE INTERFACE_LABEL_SET_VS
  END INTERFACE INTERFACE_LABEL_SET

  PUBLIC INTERFACE_MESH_ADD

  PUBLIC INTERFACE_CREATE_START, INTERFACE_CREATE_FINISH

  PUBLIC INTERFACE_COORDINATE_SYSTEM_SET

  PUBLIC INTERFACE_DESTROY, INTERFACE_MESH_CONNECTIVITY_DESTROY, InterfacePointsConnectivity_Destroy

  PUBLIC INTERFACE_LABEL_GET,INTERFACE_LABEL_SET

  PUBLIC INTERFACE_MESH_CONNECTIVITY_CREATE_START, INTERFACE_MESH_CONNECTIVITY_CREATE_FINISH

  PUBLIC INTERFACES_FINALISE,INTERFACES_INITIALISE

  PUBLIC INTERFACE_MESH_CONNECTIVITY_ELEMENT_XI_SET,InterfaceMeshConnectivity_ElementNumberSet
  
  PUBLIC INTERFACE_MESH_CONNECTIVITY_NODE_NUMBER_SET

  PUBLIC INTERFACE_MESH_CONNECTIVITY_BASIS_SET
  
  PUBLIC InterfacePointsConnectivity_CreateStart,InterfacePointsConnectivity_CreateFinish
  
  PUBLIC InterfacePointsConnectivity_DataReprojection
  
  PUBLIC InterfacePointsConnectivity_ElementNumberGet,InterfacePointsConnectivity_ElementNumberSet
  
  PUBLIC InterfacePointsConnectivity_PointXiGet,InterfacePointsConnectivity_PointXiSet
  
  PUBLIC InterfacePointsConnectivity_UpdateFromProjection
  
CONTAINS

  !
  !================================================================================================================================
  !

  SUBROUTINE INTERFACE_MESH_ADD(INTERFACE,MESH,MESH_INDEX,err,error,*)   

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to add a mesh to
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to add to the interface
    INTEGER(INTG), INTENT(OUT) :: MESH_INDEX !<On return, the index of the added mesh in the list of meshes in the interface
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: mesh_idx
    LOGICAL :: MESH_ALREADY_COUPLED
    TYPE(MESH_TYPE), POINTER :: COUPLED_MESH
    TYPE(MESH_PTR_TYPE), POINTER :: NEW_COUPLED_MESHES(:)
    TYPE(REGION_TYPE), POINTER :: COUPLED_MESH_REGION,MESH_REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    NULLIFY(NEW_COUPLED_MESHES)
    
    ENTERS("INTERFACE_MESH_ADD",err,error,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(INTERFACE%INTERFACE_FINISHED) THEN
        CALL FlagError("Interface has already been finished.",err,error,*999)
      ELSE
        IF(ASSOCIATED(MESH)) THEN
          IF(MESH%MESH_FINISHED) THEN
            MESH_REGION=>MESH%REGION
            IF(ASSOCIATED(MESH_REGION)) THEN
              ALLOCATE(NEW_COUPLED_MESHES(INTERFACE%NUMBER_OF_COUPLED_MESHES+1),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate new coupled meshes.",err,error,*999)
              !Check that the mesh is not already in the list of meshes for the interface.
              IF(INTERFACE%NUMBER_OF_COUPLED_MESHES>0) THEN
                IF(ASSOCIATED(INTERFACE%COUPLED_MESHES)) THEN
                  MESH_ALREADY_COUPLED=.FALSE.
                  DO mesh_idx=1,INTERFACE%NUMBER_OF_COUPLED_MESHES
                    COUPLED_MESH=>INTERFACE%COUPLED_MESHES(mesh_idx)%PTR
                    IF(ASSOCIATED(COUPLED_MESH)) THEN
                      COUPLED_MESH_REGION=>COUPLED_MESH%REGION
                      IF(ASSOCIATED(COUPLED_MESH_REGION)) THEN
                        IF(MESH_REGION%USER_NUMBER==COUPLED_MESH_REGION%USER_NUMBER) THEN
                          IF(MESH%USER_NUMBER==COUPLED_MESH%USER_NUMBER) THEN
                            MESH_ALREADY_COUPLED=.TRUE.
                            EXIT
                          ENDIF
                        ENDIF
                      ELSE
                        LOCAL_ERROR="Coupled interface mesh region for mesh index "// &
                          & TRIM(NUMBER_TO_VSTRING(mesh_idx,"*",err,error))//" is not associated."
                        CALL FlagError(LOCAL_ERROR,err,error,*999)
                      ENDIF
                    ELSE
                      LOCAL_ERROR="Coupled interface mesh for mesh index "//TRIM(NUMBER_TO_VSTRING(mesh_idx,"*",err,error))// &
                        & " is not associated."
                      CALL FlagError(LOCAL_ERROR,err,error,*999)
                    ENDIF
                    NEW_COUPLED_MESHES(mesh_idx)%PTR=>INTERFACE%COUPLED_MESHES(mesh_idx)%PTR
                  ENDDO !mesh_idx
                  IF(MESH_ALREADY_COUPLED) THEN
                    LOCAL_ERROR="The supplied mesh has already been added to the list of coupled meshes at mesh index "// &
                      & TRIM(NUMBER_TO_VSTRING(mesh_idx,"*",err,error))//"."
                    CALL FlagError(LOCAL_ERROR,err,error,*999)
                  ENDIF
                  DEALLOCATE(INTERFACE%COUPLED_MESHES)
                ELSE
                  CALL FlagError("Interface coupled meshes is not associated.",err,error,*999)
                ENDIF
              ENDIF
              !Add the mesh to the list of coupled meshes
              NEW_COUPLED_MESHES(INTERFACE%NUMBER_OF_COUPLED_MESHES+1)%PTR=>MESH
              INTERFACE%COUPLED_MESHES=>NEW_COUPLED_MESHES
              !Increment the number of coupled meshes and return the index
              INTERFACE%NUMBER_OF_COUPLED_MESHES=INTERFACE%NUMBER_OF_COUPLED_MESHES+1
              MESH_INDEX=INTERFACE%NUMBER_OF_COUPLED_MESHES
            ELSE
              CALL FlagError("Mesh region is not associated.",err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Mesh has not been finished.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Mesh is not associated.",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Interface is not associated.",err,error,*999)
    ENDIF
    
    EXITS("INTERFACE_MESH_ADD")
    RETURN
999 IF(ASSOCIATED(NEW_COUPLED_MESHES)) DEALLOCATE(NEW_COUPLED_MESHES)
    ERRORSEXITS("INTERFACE_MESH_ADD",err,error)
    RETURN 1
    
  END SUBROUTINE INTERFACE_MESH_ADD

  !
  !================================================================================================================================
  !

  !>Finishes the creation of an interface. \see OPENCMISS::Iron::cmfe_InterfaceCreateFinish
  SUBROUTINE INTERFACE_CREATE_FINISH(INTERFACE,err,error,*) 

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to finish the creation of
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
     
    ENTERS("INTERFACE_CREATE_FINISH",err,error,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(INTERFACE%INTERFACE_FINISHED) THEN
        CALL FlagError("Interface has already been finished.",err,error,*999)
      ELSE
        IF(INTERFACE%NUMBER_OF_COUPLED_MESHES<2) THEN
          LOCAL_ERROR="Invalid mesh coupling. Only "//TRIM(NUMBER_TO_VSTRING(INTERFACE%NUMBER_OF_COUPLED_MESHES,"*",err,error))// &
            & " have been coupled. The number of coupled meshes must be >= 2."
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        ENDIF
        INTERFACE%INTERFACE_FINISHED=.TRUE.
      ENDIF
    ELSE
      CALL FlagError("Interface is not associated.",err,error,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Interface :",err,error,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  User number = ",INTERFACE%USER_NUMBER,err,error,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Global number = ",INTERFACE%GLOBAL_NUMBER,err,error,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Label = ",INTERFACE%LABEL,err,error,*999)
      IF(ASSOCIATED(INTERFACE%INTERFACES)) THEN
        IF(ASSOCIATED(INTERFACE%INTERFACES%PARENT_REGION)) THEN
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Parent region user number = ",INTERFACE%INTERFACES% &
            & PARENT_REGION%USER_NUMBER,err,error,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Parent region label = ",INTERFACE%INTERFACES% &
            & PARENT_REGION%LABEL,err,error,*999)        
        ELSE
          CALL FlagError("Interfaces parent region is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Interface interfaces is not associated.",err,error,*999)
      ENDIF
    ENDIF
    
    EXITS("INTERFACE_CREATE_FINISH")
    RETURN
999 ERRORSEXITS("INTERFACE_CREATE_FINISH",err,error)
    RETURN 1
  END SUBROUTINE INTERFACE_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Starts the creation of an interface on a parent region. \see OPENCMISS::Iron::cmfe_InterfaceCreateStart
  SUBROUTINE INTERFACE_CREATE_START(USER_NUMBER,PARENT_REGION,INTERFACE,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the interface to create
    TYPE(REGION_TYPE), POINTER :: PARENT_REGION !<A pointer to the parent region to create the interface on.
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<On exit, a pointer to the created interface. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: interface_idx
    TYPE(INTERFACE_TYPE), POINTER :: NEW_INTERFACE
    TYPE(INTERFACE_PTR_TYPE), POINTER :: NEW_INTERFACES(:)
    TYPE(VARYING_STRING) :: LOCAL_ERROR,LOCAL_STRING

    NULLIFY(NEW_INTERFACE)
    NULLIFY(NEW_INTERFACES)
    
    ENTERS("INTERFACE_CREATE_START",err,error,*998)

    IF(ASSOCIATED(PARENT_REGION)) THEN
      IF(ASSOCIATED(INTERFACE)) THEN
        CALL FlagError("Interface is already associated.",err,error,*998)
      ELSE
        NULLIFY(INTERFACE)
        CALL INTERFACE_USER_NUMBER_FIND(USER_NUMBER,PARENT_REGION,INTERFACE,err,error,*998)
        IF(ASSOCIATED(INTERFACE)) THEN
          LOCAL_ERROR="Interface number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",err,error))// &
            & " has already been created on region number "//TRIM(NUMBER_TO_VSTRING(PARENT_REGION%USER_NUMBER,"*",err,error))//"."
          CALL FlagError(LOCAL_ERROR,err,error,*998)
        ELSE        
          NULLIFY(INTERFACE)
          !Allocate and set default interface properties.
          CALL INTERFACE_INITIALISE(NEW_INTERFACE,err,error,*999)
          NEW_INTERFACE%USER_NUMBER=USER_NUMBER
          NEW_INTERFACE%GLOBAL_NUMBER=PARENT_REGION%INTERFACES%NUMBER_OF_INTERFACES+1
          LOCAL_STRING="Interface_"//NUMBER_TO_VSTRING(USER_NUMBER,"*",err,error)
          NEW_INTERFACE%LABEL=CHAR(LOCAL_STRING)
          IF(ERR/=0) GOTO 999
          NEW_INTERFACE%INTERFACES=>PARENT_REGION%INTERFACES
          NEW_INTERFACE%PARENT_REGION=>PARENT_REGION
          !Add new initerface into list of interfaces in the parent region
          ALLOCATE(NEW_INTERFACES(PARENT_REGION%INTERFACES%NUMBER_OF_INTERFACES+1),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate new interfaces.",err,error,*999)
          DO interface_idx=1,PARENT_REGION%INTERFACES%NUMBER_OF_INTERFACES
            NEW_INTERFACES(interface_idx)%PTR=>PARENT_REGION%INTERFACES%INTERFACES(interface_idx)%PTR
          ENDDO !interface_idx
          NEW_INTERFACES(PARENT_REGION%INTERFACES%NUMBER_OF_INTERFACES+1)%PTR=>NEW_INTERFACE
          IF(ASSOCIATED(PARENT_REGION%INTERFACES%INTERFACES)) DEALLOCATE(PARENT_REGION%INTERFACES%INTERFACES)
          PARENT_REGION%INTERFACES%INTERFACES=>NEW_INTERFACES
          PARENT_REGION%INTERFACES%NUMBER_OF_INTERFACES=PARENT_REGION%INTERFACES%NUMBER_OF_INTERFACES+1
          INTERFACE=>NEW_INTERFACE
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Parent region is not associated.",err,error,*998)
    ENDIF
    
    EXITS("INTERFACE_CREATE_START")
    RETURN
999 IF(ASSOCIATED(NEW_INTERFACES)) DEALLOCATE(NEW_INTERFACES)
    CALL INTERFACE_FINALISE(INTERFACE,err,error,*998)
998 ERRORSEXITS("INTERFACE_CREATE_START",err,error)
    RETURN 1
  END SUBROUTINE INTERFACE_CREATE_START

  !
  !================================================================================================================================
  !

  !>Sets the coordinate system of an interface.  \see OPENCMISS::Iron::cmfe_Interface_CoordinateSystemSet
  SUBROUTINE INTERFACE_COORDINATE_SYSTEM_SET(INTERFACE,COORDINATE_SYSTEM,err,error,*)

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to set the coordinate system for
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM !<The coordinate system to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("INTERFACE_COORDINATE_SYSTEM_SET",err,error,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(INTERFACE%INTERFACE_FINISHED) THEN
        CALL FlagError("Interface has been finished.",err,error,*999)
      ELSE
        IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
          IF(COORDINATE_SYSTEM%COORDINATE_SYSTEM_FINISHED) THEN
            INTERFACE%COORDINATE_SYSTEM=>COORDINATE_SYSTEM
          ELSE
            CALL FlagError("Coordinate system has not been finished.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Coordinate system is not associated.",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Interface is not associated.",err,error,*999)
    ENDIF
    
    EXITS("INTERFACE_COORDINATE_SYSTEM_SET")
    RETURN
999 ERRORSEXITS("INTERFACE_COORDINATE_SYSTEM_SET",err,error)
    RETURN 1
  END SUBROUTINE INTERFACE_COORDINATE_SYSTEM_SET
  
!
  !================================================================================================================================
  !

  !>Destroys an interface. \see OPENCMISS::Iron::cmfe_InterfaceDestroy
  SUBROUTINE INTERFACE_DESTROY(INTERFACE,err,error,*) 

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: interface_idx,interface_position
    TYPE(INTERFACE_PTR_TYPE), POINTER :: NEW_INTERFACES(:)
    TYPE(INTERFACES_TYPE), POINTER :: INTERFACES
     
    NULLIFY(NEW_INTERFACES)

    ENTERS("INTERFACE_DESTROY",err,error,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      INTERFACES=>INTERFACE%INTERFACES
      IF(ASSOCIATED(INTERFACES)) THEN
        interface_position=INTERFACE%GLOBAL_NUMBER

        !Destroy all the interface condition components
        CALL INTERFACE_FINALISE(INTERFACE,err,error,*999)
        
        !Remove the interface from the list of interfaces
        IF(INTERFACES%NUMBER_OF_INTERFACES>1) THEN
          ALLOCATE(NEW_INTERFACES(INTERFACES%NUMBER_OF_INTERFACES-1),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate new interface conditions.",err,error,*999)
          DO interface_idx=1,INTERFACES%NUMBER_OF_INTERFACES
            IF(interface_idx<interface_position) THEN
              NEW_INTERFACES(interface_idx)%PTR=>INTERFACES%INTERFACES(interface_idx)%PTR
            ELSE IF(interface_idx>interface_position) THEN
              INTERFACES%INTERFACES(interface_idx)%PTR%GLOBAL_NUMBER=INTERFACES%INTERFACES(interface_idx)%PTR%GLOBAL_NUMBER-1
              NEW_INTERFACES(interface_idx-1)%PTR=>INTERFACES%INTERFACES(interface_idx)%PTR
            ENDIF
          ENDDO !interface_idx
          IF(ASSOCIATED(INTERFACES%INTERFACES)) DEALLOCATE(INTERFACES%INTERFACES)
          INTERFACES%INTERFACES=>NEW_INTERFACES
          INTERFACES%NUMBER_OF_INTERFACES=INTERFACES%NUMBER_OF_INTERFACES-1
        ELSE
          DEALLOCATE(INTERFACES%INTERFACES)
          INTERFACES%NUMBER_OF_INTERFACES=0
        ENDIF
        
      ELSE
        CALL FlagError("Interface interfaces is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Interface is not associated.",err,error,*999)
    ENDIF    
    
    EXITS("INTERFACE_DESTROY")
    RETURN
999 ERRORSEXITS("INTERFACE_DESTROY",err,error)
    RETURN 1
  END SUBROUTINE INTERFACE_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalises an interface and deallocates all memory.
  SUBROUTINE INTERFACE_FINALISE(INTERFACE,err,error,*) 

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
     
    ENTERS("INTERFACE_FINALISE",err,error,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(ASSOCIATED(INTERFACE%COUPLED_MESHES)) DEALLOCATE(INTERFACE%COUPLED_MESHES)
      IF(ASSOCIATED(INTERFACE%MESH_CONNECTIVITY)) &
         & CALL INTERFACE_MESH_CONNECTIVITY_FINALISE(INTERFACE%MESH_CONNECTIVITY,err,error,*999)
      IF(ASSOCIATED(INTERFACE%pointsConnectivity)) &
        & CALL InterfacePointsConnectivity_Finalise(INTERFACE%pointsConnectivity,err,error,*999)
      IF(ASSOCIATED(INTERFACE%NODES)) CALL NODES_DESTROY(INTERFACE%NODES,err,error,*999)
      CALL GENERATED_MESHES_FINALISE(INTERFACE%generatedMeshes,err,error,*999)
      CALL MESHES_FINALISE(INTERFACE%MESHES,err,error,*999)
      CALL FIELDS_FINALISE(INTERFACE%FIELDS,err,error,*999)
      CALL INTERFACE_CONDITIONS_FINALISE(INTERFACE%INTERFACE_CONDITIONS,err,error,*999)
      CALL DataPointSets_Finalise(interface%dataPointSets,err,error,*999)
      DEALLOCATE(INTERFACE)
    ENDIF
    
    EXITS("INTERFACE_FINALISE")
    RETURN
999 ERRORSEXITS("INTERFACE_FINALISE",err,error)
    RETURN 1
  END SUBROUTINE INTERFACE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises an interface.
  SUBROUTINE INTERFACE_INITIALISE(INTERFACE,err,error,*) 

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to initialise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
     
    ENTERS("INTERFACE_INITIALISE",err,error,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      CALL FlagError("Interface is already associated.",err,error,*999)
    ELSE
      ALLOCATE(INTERFACE,STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate interface.",err,error,*999)
      INTERFACE%USER_NUMBER=0
      INTERFACE%GLOBAL_NUMBER=0
      INTERFACE%INTERFACE_FINISHED=.FALSE.
      INTERFACE%LABEL=""
      NULLIFY(INTERFACE%INTERFACES)
      NULLIFY(INTERFACE%PARENT_REGION)
      INTERFACE%NUMBER_OF_COUPLED_MESHES=0
      NULLIFY(INTERFACE%COUPLED_MESHES)
      NULLIFY(INTERFACE%MESH_CONNECTIVITY)
      NULLIFY(INTERFACE%pointsConnectivity)
      NULLIFY(INTERFACE%NODES)
      NULLIFY(INTERFACE%MESHES)
      NULLIFY(INTERFACE%generatedMeshes)
      NULLIFY(INTERFACE%FIELDS)
      NULLIFY(INTERFACE%INTERFACE_CONDITIONS)
      NULLIFY(INTERFACE%COORDINATE_SYSTEM)
      NULLIFY(INTERFACE%dataPointSets)
      CALL DataPointSets_Initialise(interface,err,error,*999)
      CALL MESHES_INITIALISE(INTERFACE,err,error,*999)
      CALL GENERATED_MESHES_INITIALISE(INTERFACE,err,error,*999)
      CALL FIELDS_INITIALISE(INTERFACE,err,error,*999)
      CALL INTERFACE_CONDITIONS_INITIALISE(INTERFACE,err,error,*999)
    ENDIF
    
    EXITS("INTERFACE_INITIALISE")
    RETURN
999 ERRORSEXITS("INTERFACE_INITIALISE",err,error)
    RETURN 1
  END SUBROUTINE INTERFACE_INITIALISE

  !
  !================================================================================================================================
  !

  !>Returns the label of an interface for a character label. \see OPENCMISS::Iron::cmfe_InterfaceLabelGet
  SUBROUTINE INTERFACE_LABEL_GET_C(INTERFACE,LABEL,err,error,*)

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to get the label for
    CHARACTER(LEN=*), INTENT(OUT) :: LABEL !<On return the interface label.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: C_LENGTH,VS_LENGTH

    ENTERS("INTERFACE_LABEL_GET_C",err,error,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      C_LENGTH=LEN(LABEL)
      VS_LENGTH=LEN_TRIM(INTERFACE%LABEL)
      IF(C_LENGTH>VS_LENGTH) THEN
        LABEL=CHAR(INTERFACE%LABEL,VS_LENGTH)
      ELSE
        LABEL=CHAR(INTERFACE%LABEL,C_LENGTH)
      ENDIF
    ELSE
      CALL FlagError("Interface is not associated.",err,error,*999)
    ENDIF
    
    EXITS("INTERFACE_LABEL_GET_C")
    RETURN
999 ERRORSEXITS("INTERFACE_LABEL_GET_C",err,error)
    RETURN 1
    
  END SUBROUTINE INTERFACE_LABEL_GET_C

   !
  !================================================================================================================================
  !

  !>Returns the label of an interface for a varying string label. \see OPENCMISS::Iron::cmfe_InterfaceLabelGet
  SUBROUTINE INTERFACE_LABEL_GET_VS(INTERFACE,LABEL,err,error,*)

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to get the label for
    TYPE(VARYING_STRING), INTENT(OUT) :: LABEL !<On return the interface label.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("INTERFACE_LABEL_GET_VS",err,error,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      !\todo The following line crashes the AIX compiler unless it has a VAR_STR(CHAR()) around it
      LABEL=VAR_STR(CHAR(INTERFACE%LABEL))
    ELSE
      CALL FlagError("Interface is not associated.",err,error,*999)
    ENDIF
    
    EXITS("INTERFACE_LABEL_GET_VS")
    RETURN
999 ERRORSEXITS("INTERFACE_LABEL_GET_VS",err,error)
    RETURN 1
    
  END SUBROUTINE INTERFACE_LABEL_GET_VS

  !
  !================================================================================================================================
  !

  !>Sets the label of an interface for a character label. \see OPENCMISS::Iron::cmfe_InterfaceLabelSet
  SUBROUTINE INTERFACE_LABEL_SET_C(INTERFACE,LABEL,err,error,*)

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to set the label for 
    CHARACTER(LEN=*), INTENT(IN) :: LABEL !<The label to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("INTERFACE_LABEL_SET_C",err,error,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(INTERFACE%INTERFACE_FINISHED) THEN
        CALL FlagError("Interface has been finished.",err,error,*999)
      ELSE
        INTERFACE%LABEL=LABEL
      ENDIF
    ELSE
      CALL FlagError("Interface is not associated.",err,error,*999)
    ENDIF
    
    EXITS("INTERFACE_LABEL_SET_C")
    RETURN
999 ERRORSEXITS("INTERFACE_LABEL_SET_C",err,error)
    RETURN 1
  END SUBROUTINE INTERFACE_LABEL_SET_C

  !
  !================================================================================================================================
  !

  !>Sets the label of an interface for a varying string label. \see OPENCMISS::Iron::cmfe_InterfaceLabelSet
  SUBROUTINE INTERFACE_LABEL_SET_VS(INTERFACE,LABEL,err,error,*)

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to set the label for 
    TYPE(VARYING_STRING), INTENT(IN) :: LABEL !<The label to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("INTERFACE_LABEL_SET_VS",err,error,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(INTERFACE%INTERFACE_FINISHED) THEN
        CALL FlagError("Interface has been finished.",err,error,*999)
      ELSE
        INTERFACE%LABEL=LABEL
      ENDIF
    ELSE
      CALL FlagError("Interface is not associated.",err,error,*999)
    ENDIF
    
    EXITS("INTERFACE_LABEL_SET_VS")
    RETURN
999 ERRORSEXITS("INTERFACE_LABEL_SET_VS",err,error)
    RETURN 1
  END SUBROUTINE INTERFACE_LABEL_SET_VS

  !
  !================================================================================================================================
  !

  !>Finalises a meshes connectivity for an interface.
  SUBROUTINE INTERFACE_MESH_CONNECTIVITY_CREATE_FINISH(INTERFACE_MESH_CONNECTIVITY,err,error,*) 

    !Argument variables
    TYPE(INTERFACE_MESH_CONNECTIVITY_TYPE), POINTER :: INTERFACE_MESH_CONNECTIVITY !<A pointer to the interface meshes connectivity to finish creating
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: InterfaceElementIdx,CoupledMeshIdx
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("INTERFACE_MESH_CONNECTIVITY_CREATE_FINISH",err,error,*999)

    IF(ASSOCIATED(INTERFACE_MESH_CONNECTIVITY)) THEN
      IF(INTERFACE_MESH_CONNECTIVITY%MESH_CONNECTIVITY_FINISHED) THEN
        CALL FlagError("Interface meshes connectivity has already been finished.",err,error,*999)
      ELSE
        !Check if connectivity from each interface element to an appropriate element in the coupled meshes has been setup
        DO InterfaceElementIdx=1,INTERFACE_MESH_CONNECTIVITY%NUMBER_OF_INTERFACE_ELEMENTS
          DO CoupledMeshIdx=1,INTERFACE_MESH_CONNECTIVITY%NUMBER_OF_COUPLED_MESHES
            IF (INTERFACE_MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY(InterfaceElementIdx,CoupledMeshIdx)% &
              & COUPLED_MESH_ELEMENT_NUMBER==0) THEN
              LOCAL_ERROR="The connectivity from interface element " &
                //TRIM(NUMBER_TO_VSTRING(InterfaceElementIdx,"*",err,error))//" to an element in coupled mesh " & 
                //TRIM(NUMBER_TO_VSTRING(CoupledMeshIdx,"*",err,error))//" has not been defined."
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            ENDIF
          ENDDO !CoupledMeshIdx
        ENDDO !InterfaceElementIdx
        !Calculate line or face numbers for coupled mesh elements that are connected to the interface mesh
        CALL InterfaceMeshConnectivity_ConnectedLinesFacesCalculate(INTERFACE_MESH_CONNECTIVITY,err,error,*999)
        INTERFACE_MESH_CONNECTIVITY%MESH_CONNECTIVITY_FINISHED=.TRUE.
      ENDIF
    ELSE
      CALL FlagError("Interface meshes connectivity is not associated.",err,error,*999)
    ENDIF

    EXITS("INTERFACE_MESH_CONNECTIVITY_CREATE_FINISH")
    RETURN
999 ERRORSEXITS("INTERFACE_MESH_CONNECTIVITY_CREATE_FINISH",err,error)
    RETURN 1
  END SUBROUTINE INTERFACE_MESH_CONNECTIVITY_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Initialises a meshes connectivity for an interface.
  SUBROUTINE INTERFACE_MESH_CONNECTIVITY_CREATE_START(INTERFACE,MESH,INTERFACE_MESH_CONNECTIVITY,err,error,*) 

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to create the meshes connectivity for
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(INTERFACE_MESH_CONNECTIVITY_TYPE), POINTER :: INTERFACE_MESH_CONNECTIVITY !<On return, a pointer to the created meshes connectivity
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("INTERFACE_MESH_CONNECTIVITY_CREATE_START",err,error,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(INTERFACE%INTERFACE_FINISHED) THEN
        IF(ASSOCIATED(INTERFACE%MESH_CONNECTIVITY)) THEN
          CALL FlagError("The interface already has a meshes connectivity associated.",err,error,*999)
        ELSE
          !Initialise the meshes connectivity
          CALL INTERFACE_MESH_CONNECTIVITY_INITIALISE(INTERFACE,MESH,err,error,*999)
          !Return the pointer
          INTERFACE_MESH_CONNECTIVITY=>INTERFACE%MESH_CONNECTIVITY
        ENDIF
      ELSE
        CALL FlagError("Interface has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Interface is not associated.",err,error,*999)
    ENDIF
    
    EXITS("INTERFACE_MESH_CONNECTIVITY_CREATE_START")
    RETURN
999 ERRORSEXITS("INTERFACE_MESH_CONNECTIVITY_CREATE_START",err,error)
    RETURN 1
    
  END SUBROUTINE INTERFACE_MESH_CONNECTIVITY_CREATE_START

  !
  !================================================================================================================================
  !

  !>Finalises a meshes connectivity and deallocates all memory
  SUBROUTINE INTERFACE_MESH_CONNECTIVITY_DESTROY(INTERFACE_MESH_CONNECTIVITY,err,error,*) 

    !Argument variables
    TYPE(INTERFACE_MESH_CONNECTIVITY_TYPE), POINTER :: INTERFACE_MESH_CONNECTIVITY !<A pointer to the interface meshes connectivity to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("INTERFACE_MESH_CONNECTIVITY_DESTROY",err,error,*999)

    IF(ASSOCIATED(INTERFACE_MESH_CONNECTIVITY)) THEN
      CALL INTERFACE_MESH_CONNECTIVITY_FINALISE(INTERFACE_MESH_CONNECTIVITY,err,error,*999)
    ELSE
      CALL FlagError("Interface meshes connectivity is not associated.",err,error,*999)
    ENDIF
       
    EXITS("INTERFACE_MESH_CONNECTIVITY_DESTROY")
    RETURN
999 ERRORSEXITS("INTERFACE_MESH_CONNECTIVITY_DESTROY",err,error)
    RETURN 1
  END SUBROUTINE INTERFACE_MESH_CONNECTIVITY_DESTROY

  !
  !================================================================================================================================
  !

  !>Sets the interface mesh connectivity basis
  SUBROUTINE INTERFACE_MESH_CONNECTIVITY_BASIS_SET(INTERFACE_MESH_CONNECTIVITY,BASIS,err,error,*)

    !Argument variables
    TYPE(INTERFACE_MESH_CONNECTIVITY_TYPE), POINTER :: INTERFACE_MESH_CONNECTIVITY !<A pointer to interface mesh connectivity to set the element number of elements for.
    TYPE(BASIS_TYPE), POINTER :: BASIS
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(INTERFACE_ELEMENT_CONNECTIVITY_TYPE), POINTER :: ELEMENT_CONNECTIVITY
    INTEGER(INTG) :: InterfaceElementIdx,CoupledMeshIdx,NumberOfInterfaceElementNodes,NumberOfCoupledMeshXiDirections

    ENTERS("INTERFACE_MESH_CONNECTIVITY_BASIS_SET",err,error,*999)

    IF(ASSOCIATED(INTERFACE_MESH_CONNECTIVITY)) THEN
      IF(INTERFACE_MESH_CONNECTIVITY%MESH_CONNECTIVITY_FINISHED) THEN
        CALL FlagError("Interface mesh connectivity already been finished.",err,error,*999)
      ELSE
        IF(ASSOCIATED(INTERFACE_MESH_CONNECTIVITY%BASIS)) THEN
          CALL FlagError("Mesh connectivity basis already associated.",err,error,*999)
        ELSE
          IF(ASSOCIATED(BASIS)) THEN
            INTERFACE_MESH_CONNECTIVITY%BASIS=>BASIS
            !Now that the mesh connectivity basis is set the number of interface element nodes can be determined and now MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY(InterfaceElementIdx,CoupledMeshIdx)%XI can be allocated
            !\todo NumberOfCoupledMeshXiDirections currently set to the number of interface mesh xi directions + 1. Restructure ELEMENT_CONNECTIVITY type see below
            NumberOfCoupledMeshXiDirections=INTERFACE_MESH_CONNECTIVITY%INTERFACE_MESH%NUMBER_OF_DIMENSIONS+1
            NumberOfInterfaceElementNodes=INTERFACE_MESH_CONNECTIVITY%BASIS%NUMBER_OF_NODES
            DO InterfaceElementIdx=1,INTERFACE_MESH_CONNECTIVITY%NUMBER_OF_INTERFACE_ELEMENTS
              DO CoupledMeshIdx=1,INTERFACE_MESH_CONNECTIVITY%NUMBER_OF_COUPLED_MESHES
                ELEMENT_CONNECTIVITY=>INTERFACE_MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY(InterfaceElementIdx,CoupledMeshIdx)
                !\todo Update mesh component index to look at the number of mesh components in each element. 
                !\todo Currently this defaults to the first mesh component ie %XI(NumberOfInterfaceMeshXi,1,NumberOfInterfaceElementNodes)). 
                !\todo The interface mesh types will also need to be restructured.
                !eg. INTERFACE%MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY(InterfaceElementIdx,CoupledMeshIdx)%MESH_COMPONENT(MeshComponentIdx)%XI(NumberOfCoupledMeshXiDirections,NumberOfInterfaceElementNodes) and adding appropriate initialize/finialize routines
                ALLOCATE(ELEMENT_CONNECTIVITY%XI(NumberOfCoupledMeshXiDirections,1,NumberOfInterfaceElementNodes),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate interface element connectivity.",err,error,*999)
                ELEMENT_CONNECTIVITY%XI=0.0_DP
              ENDDO !CoupledMeshIdx
            ENDDO !InterfaceElementIdx
          ELSE
            CALL FlagError("Basis to set mesh connectivity not associated.",err,error,*999)
          ENDIF
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Interface mesh connectivity is not associated.",err,error,*999)
    ENDIF

    EXITS("INTERFACE_MESH_CONNECTIVITY_BASIS_SET")
    RETURN
999 ERRORSEXITS("INTERFACE_MESH_CONNECTIVITY_BASIS_SET",err,error)
    RETURN 1
    
  END SUBROUTINE INTERFACE_MESH_CONNECTIVITY_BASIS_SET

  !
  !================================================================================================================================
  !
    
  !>Sets the mapping from an xi position of a coupled mesh element to a node of an interface mesh element
  SUBROUTINE INTERFACE_MESH_CONNECTIVITY_ELEMENT_XI_SET(INTERFACE_MESH_CONNECTIVITY,INTERFACE_MESH_ELEMENT_NUMBER, &
    & COUPLED_MESH_INDEX,COUPLED_MESH_ELEMENT_NUMBER,INTERFACE_MESH_LOCAL_NODE_NUMBER,INTERFACE_MESH_COMPONENT_NUMBER,XI, &
    & err,error,*)

    !Argument variables
    TYPE(INTERFACE_MESH_CONNECTIVITY_TYPE), POINTER :: INTERFACE_MESH_CONNECTIVITY !<A pointer to the interface mesh connectivity for the interface mesh
    INTEGER(INTG), INTENT(IN) :: INTERFACE_MESH_ELEMENT_NUMBER !<The interface mesh element number to which the specified coupled mesh element would be connected
    INTEGER(INTG), INTENT(IN) :: COUPLED_MESH_INDEX !<The index of the coupled mesh at the interface to set the element connectivity for
    INTEGER(INTG), INTENT(IN) :: COUPLED_MESH_ELEMENT_NUMBER !<The coupled mesh element to define the element xi connectivity from
    INTEGER(INTG), INTENT(IN) :: INTERFACE_MESH_LOCAL_NODE_NUMBER !<The interface mesh node to assign the coupled mesh element xi to
    INTEGER(INTG), INTENT(IN) :: INTERFACE_MESH_COMPONENT_NUMBER !<The interface mesh node's component to assign the coupled mesh element xi to
    REAL(DP), INTENT(IN) :: XI(:) !<XI(xi_idx). The xi value for the xi_idx'th xi direction in the coupled mesh element.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: coupled_mesh_number_of_xi
    TYPE(INTERFACE_ELEMENT_CONNECTIVITY_TYPE), POINTER :: ELEMENT_CONNECTIVITY
    
    ENTERS("INTERFACE_MESH_CONNECTIVITY_ELEMENT_XI_SET",err,error,*999)

    IF(ASSOCIATED(INTERFACE_MESH_CONNECTIVITY)) THEN
      IF(INTERFACE_MESH_CONNECTIVITY%MESH_CONNECTIVITY_FINISHED) THEN
        CALL FlagError("Interface mesh connectivity already been finished.",err,error,*999)
      ELSE
        IF(ALLOCATED(INTERFACE_MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY)) THEN
          IF((INTERFACE_MESH_ELEMENT_NUMBER>0).AND. &
            & (INTERFACE_MESH_ELEMENT_NUMBER<=INTERFACE_MESH_CONNECTIVITY%NUMBER_OF_INTERFACE_ELEMENTS)) THEN
            IF((COUPLED_MESH_INDEX>0).AND.(COUPLED_MESH_INDEX<=INTERFACE_MESH_CONNECTIVITY%NUMBER_OF_COUPLED_MESHES)) THEN
              IF((COUPLED_MESH_ELEMENT_NUMBER>0).AND.(COUPLED_MESH_ELEMENT_NUMBER<= &
                & INTERFACE_MESH_CONNECTIVITY%INTERFACE%COUPLED_MESHES(COUPLED_MESH_INDEX)%PTR%NUMBER_OF_ELEMENTS))THEN
                IF((INTERFACE_MESH_COMPONENT_NUMBER>0).AND. &
                  & (INTERFACE_MESH_COMPONENT_NUMBER<=INTERFACE_MESH_CONNECTIVITY%INTERFACE_MESH%NUMBER_OF_COMPONENTS)) THEN
                  IF((INTERFACE_MESH_LOCAL_NODE_NUMBER>0).AND.(INTERFACE_MESH_LOCAL_NODE_NUMBER<= &
                    & INTERFACE_MESH_CONNECTIVITY%BASIS%NUMBER_OF_NODES))THEN
                    ELEMENT_CONNECTIVITY=>INTERFACE_MESH_CONNECTIVITY% &
                      & ELEMENT_CONNECTIVITY(INTERFACE_MESH_ELEMENT_NUMBER,COUPLED_MESH_INDEX)
                    IF(ELEMENT_CONNECTIVITY%COUPLED_MESH_ELEMENT_NUMBER==COUPLED_MESH_ELEMENT_NUMBER)THEN
                      coupled_mesh_number_of_xi = INTERFACE_MESH_CONNECTIVITY%INTERFACE%COUPLED_MESHES(COUPLED_MESH_INDEX)%PTR% & 
                        & TOPOLOGY(1)%PTR%ELEMENTS%ELEMENTS(COUPLED_MESH_ELEMENT_NUMBER)%BASIS%NUMBER_OF_XI
                      IF(SIZE(XI)==coupled_mesh_number_of_xi) THEN
                        !\todo the ELEMENT_CONNECTIVITY%XI array needs to be restructured to efficiently utilize memory when coupling bodies with 2xi directions to bodies with 3xi directions using an interface.
                        ELEMENT_CONNECTIVITY%XI(1:coupled_mesh_number_of_xi,INTERFACE_MESH_COMPONENT_NUMBER, &
                          & INTERFACE_MESH_LOCAL_NODE_NUMBER)=XI(1:coupled_mesh_number_of_xi)
                      ELSE
                        CALL FlagError("The size of the xi array provided does not match the coupled mesh element's' number"// &
                          & " of xi.",err,error,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Coupled mesh element number doesn't match that set to the interface.",err,error,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Interface local node number is out of range.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Interface component number is out of range.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Coupled mesh element number out of range.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Interface coupled mesh index number out of range.",err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Interface mesh element number out of range.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Interface elements connectivity array not allocated.",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Interface mesh connectivity is not associated.",err,error,*999)
    ENDIF

    EXITS("INTERFACE_MESH_CONNECTIVITY_ELEMENT_XI_SET")
    RETURN
999 ERRORSEXITS("INTERFACE_MESH_CONNECTIVITY_ELEMENT_XI_SET",err,error)
    RETURN 1
    
  END SUBROUTINE INTERFACE_MESH_CONNECTIVITY_ELEMENT_XI_SET

  !
  !================================================================================================================================
  !
  
  !>Sets the connectivity between an element in a coupled mesh to an element in the interface mesh
  SUBROUTINE InterfaceMeshConnectivity_ElementNumberSet(INTERFACE_MESH_CONNECTIVITY,INTERFACE_MESH_ELEMENT_NUMBER, &
      & COUPLED_MESH_INDEX,COUPLED_MESH_ELEMENT_NUMBER,err,error,*)

    !Argument variables
    TYPE(INTERFACE_MESH_CONNECTIVITY_TYPE), POINTER :: INTERFACE_MESH_CONNECTIVITY !<A pointer to the interface mesh connectivity for the interface mesh
    INTEGER(INTG), INTENT(IN) :: INTERFACE_MESH_ELEMENT_NUMBER !<The interface mesh element number to which the specified coupled mesh element would be connected
    INTEGER(INTG), INTENT(IN) :: COUPLED_MESH_INDEX !<The index of the coupled mesh at the interface to set the element connectivity for
    INTEGER(INTG), INTENT(IN) :: COUPLED_MESH_ELEMENT_NUMBER !<The coupled mesh element to be connected to the interface
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(INTERFACE_ELEMENT_CONNECTIVITY_TYPE), POINTER :: ELEMENT_CONNECTIVITY
    
    ENTERS("InterfaceMeshConnectivity_ElementNumberSet",err,error,*999)

    IF(ASSOCIATED(INTERFACE_MESH_CONNECTIVITY)) THEN
      IF(INTERFACE_MESH_CONNECTIVITY%MESH_CONNECTIVITY_FINISHED) THEN
        CALL FlagError("Interface mesh connectivity has already been finished.",err,error,*999)
      ELSE
        IF((INTERFACE_MESH_ELEMENT_NUMBER>0).AND.(INTERFACE_MESH_ELEMENT_NUMBER<= &
          & INTERFACE_MESH_CONNECTIVITY%NUMBER_OF_INTERFACE_ELEMENTS)) THEN
          IF((COUPLED_MESH_INDEX>0).AND.(COUPLED_MESH_INDEX<=INTERFACE_MESH_CONNECTIVITY%NUMBER_OF_COUPLED_MESHES)) THEN
            IF (ALLOCATED(INTERFACE_MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY)) THEN
              ELEMENT_CONNECTIVITY=>INTERFACE_MESH_CONNECTIVITY% &
                & ELEMENT_CONNECTIVITY(INTERFACE_MESH_ELEMENT_NUMBER,COUPLED_MESH_INDEX)
              ELEMENT_CONNECTIVITY%COUPLED_MESH_ELEMENT_NUMBER=COUPLED_MESH_ELEMENT_NUMBER
            ELSE
              CALL FlagError("Interface elements connectivity array not allocated.",err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Interface coupled mesh index number out of range.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Interface mesh element number out of range.",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Interface mesh connectivity is not associated.",err,error,*999)
    ENDIF
    
    EXITS("InterfaceMeshConnectivity_ElementNumberSet")
    RETURN
999 ERRORSEXITS("InterfaceMeshConnectivity_ElementNumberSet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMeshConnectivity_ElementNumberSet

  !
  !================================================================================================================================
  !
  
  !>Sets the connectivity between an element in a coupled mesh to an element in the interface mesh
  SUBROUTINE INTERFACE_MESH_CONNECTIVITY_NODE_NUMBER_SET(NODES,INTERFACE_MESH_NODE_NUMBERS, &
      & FIRST_COUPLED_MESH_INDEX,FIRST_COUPLED_MESH_NODE_NUMBERS,SECOND_COUPLED_MESH_INDEX,SECOND_COUPLED_MESH_NODE_NUMBERS, &
      & err,error,*)        
        
    !Argument variables
    TYPE(NODES_TYPE), POINTER :: NODES !<A pointer to the interface mesh connectivity for the interface mesh
    INTEGER(INTG), INTENT(IN) :: INTERFACE_MESH_NODE_NUMBERS(:) !<The interface mesh element number to which the specified coupled mesh element would be connected
    INTEGER(INTG), INTENT(IN) :: FIRST_COUPLED_MESH_INDEX,SECOND_COUPLED_MESH_INDEX !<The index of the coupled mesh at the interface to set the element connectivity for
    INTEGER(INTG), INTENT(IN) :: FIRST_COUPLED_MESH_NODE_NUMBERS(:),SECOND_COUPLED_MESH_NODE_NUMBERS(:) !<The coupled mesh element to be connected to the interface
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: nodeIndex
    
    ENTERS("INTERFACE_MESH_CONNECTIVITY_NODE_NUMBER_SET",err,error,*999)

    IF(ASSOCIATED(NODES)) THEN
      IF(NODES%INTERFACE%MESH_CONNECTIVITY%MESH_CONNECTIVITY_FINISHED) THEN
        PRINT *, 'CHECK how to circumvent! interface_routines.f90:1133'
        CALL FlagError("Interface mesh connectivity has already been finished.",err,error,*999)
      ELSE
        !Default to two coupled meshes
        ALLOCATE(NODES%COUPLED_NODES(2,SIZE(INTERFACE_MESH_NODE_NUMBERS(:))))
        DO nodeIndex=1,SIZE(INTERFACE_MESH_NODE_NUMBERS(:))
          NODES%COUPLED_NODES(FIRST_COUPLED_MESH_INDEX,nodeIndex)=FIRST_COUPLED_MESH_NODE_NUMBERS(nodeIndex)
          NODES%COUPLED_NODES(SECOND_COUPLED_MESH_INDEX,nodeIndex)=SECOND_COUPLED_MESH_NODE_NUMBERS(nodeIndex)
        ENDDO
      ENDIF
    ELSE
      CALL FlagError("Nodes are not associated.",err,error,*999)
    ENDIF
    
    EXITS("INTERFACE_MESH_CONNECTIVITY_NODE_NUMBER_SET")
    RETURN
999 ERRORSEXITS("INTERFACE_MESH_CONNECTIVITY_NODE_NUMBER_SET",err,error)
    RETURN 1
    
  END SUBROUTINE INTERFACE_MESH_CONNECTIVITY_NODE_NUMBER_SET

  !
  !================================================================================================================================
  !
  
  !>Calculate line or face numbers for coupled mesh elements that are connected to the interface mesh
  SUBROUTINE InterfaceMeshConnectivity_ConnectedLinesFacesCalculate(interfaceMeshConnectivity,err,error,*)

    !Argument variables
    TYPE(INTERFACE_MESH_CONNECTIVITY_TYPE), POINTER :: interfaceMeshConnectivity !<A pointer to interface mesh connectivity to calculate line or face numbers for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    REAL(DP) :: areaCoordinates(4,4),xi,xiCoordinates(4,4),xiDifference
    INTEGER(INTG) :: constantXiIdx,coupledMeshIdx,coupledMeshXiIdx,interfaceElementIdx,localNodeIdx,numberOfInterfaceElementNodes, &
      & numberOfInterfaceMeshXi,numberOfMeshXi,numberOfMeshXiCoordinates,xiCoordIdx,zeroXiCoordIdx
    LOGICAL :: differentXi
    TYPE(BASIS_TYPE), POINTER :: interfaceBasis
    TYPE(INTERFACE_ELEMENT_CONNECTIVITY_TYPE), POINTER :: elementConnectivity
    TYPE(MESH_TYPE), POINTER :: coupledMesh, interfaceMesh
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("InterfaceMeshConnectivity_ConnectedLinesFacesCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceMeshConnectivity)) CALL FlagError("Interface mesh connectivity is not associated.",err,error,*999)

    interfaceBasis=>interfaceMeshConnectivity%basis
    IF(.NOT.ASSOCIATED(interfaceBasis)) CALL FlagError("Interface basis is not associated.",err,error,*999)
    interfaceMesh=>interfaceMeshConnectivity%INTERFACE_MESH
    numberOfInterfaceElementNodes=interfaceBasis%NUMBER_OF_NODES
    numberOfInterfaceMeshXi=interfaceMesh%NUMBER_OF_DIMENSIONS     
    
!!\TODO: We need to check the basis in order to interpret the xi and decided to use area coordinates of normal coordinates. This really should be via the coupled mesh element being considered. For know assume either a LHTP basis or a simplex basis based on the interface basis. 
    SELECT CASE(interfaceBasis%TYPE)
    CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
      DO coupledMeshIdx=1,interfaceMeshConnectivity%NUMBER_OF_COUPLED_MESHES
        coupledMesh=>interfaceMeshConnectivity%INTERFACE%COUPLED_MESHES(coupledMeshIdx)%ptr
        IF(.NOT.ASSOCIATED(coupledMesh)) THEN
          localError="The interface coupled mesh for coupled mesh index "// &
            & TRIM(NumberToVString(coupledMeshIdx,"*",err,error))//" is not associated."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        numberOfMeshXi=coupledMesh%NUMBER_OF_DIMENSIONS
        DO interfaceElementIdx=1,interfaceMeshConnectivity%NUMBER_OF_INTERFACE_ELEMENTS
          elementConnectivity=>interfaceMeshConnectivity%ELEMENT_CONNECTIVITY(interfaceElementIdx,coupledMeshIdx)
          SELECTCASE(numberOfInterfaceMeshXi)
          CASE(1) !Lines
            DO coupledMeshXiIdx=1,numberOfMeshXi
              ! Calculate difference between first node and last node of an element
              xiDifference=elementConnectivity%xi(coupledMeshXiIdx,1,1)- &
                & elementConnectivity%xi(coupledMeshXiIdx,1,numberOfInterfaceElementNodes)              
              IF(ABS(xiDifference)<ZERO_TOLERANCE) THEN
                IF(ABS(elementConnectivity%xi(coupledMeshXiIdx,1,numberOfInterfaceElementNodes))<ZERO_TOLERANCE) THEN
                  elementConnectivity%CONNECTED_LINE=3-(coupledMeshXiIdx-1)*2
                ELSE
                  elementConnectivity%CONNECTED_LINE=4-(coupledMeshXiIdx-1)*2
                ENDIF
              ENDIF
            ENDDO
          CASE(2) !Faces
            constantXiIdx=0
            DO coupledMeshXiIdx=1,numberOfMeshXi
              xi=elementConnectivity%xi(coupledMeshXiIdx,1,1)
              differentXi=.FALSE.
              DO localNodeIdx=2,numberOfInterfaceElementNodes
                IF(ABS(elementConnectivity%xi(coupledMeshXiIdx,1,localNodeIdx)-xi)>ZERO_TOLERANCE) THEN 
                  differentXi=.TRUE.
                  EXIT
                ENDIF
              ENDDO !localNodeIdx
              IF(.NOT.differentXi) THEN
                IF(constantXiIdx==0) THEN
                  constantXiIdx=coupledMeshXiIdx
                ELSE
                  localError="Invalid xi connectivity. Xi directions "//TRIM(NumberToVString(coupledMeshXiIdx,"*",err,error))// &
                    & " and "//TRIM(NumberToVString(constantXiIdx,"*",err,error))// &
                    & " are both constant acrosss all interface nodes for coupled mesh number "// &
                    & TRIM(NumberToVString(coupledMeshIdx,"*",err,error))//"."
                  CALL FlagError(localError,err,error,*999)
                ENDIF
              ENDIF
            ENDDO !coupledMeshXiIdx
            IF(constantXiIdx==0) THEN
              localError="Invalid xi connectivity. There is no xi direction that is constant across all interface nodes "// &
                & "for coupled mesh number "//TRIM(NumberToVString(coupledMeshIdx,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ELSE
              IF(ABS(elementConnectivity%xi(constantXiIdx,1,1))<ZERO_TOLERANCE) THEN
                elementConnectivity%CONNECTED_FACE=constantXiIdx
              ELSE
                elementConnectivity%CONNECTED_FACE=constantXiIdx+3
              ENDIF
            ENDIF
          CASE DEFAULT 
            localError="The number of interface mesh dimension of "// &
              & TRIM(NUMBER_TO_VSTRING(numberOfInterfaceMeshXi,"*",err,error))//" is invalid"
            CALL FlagError(localError,err,error,*999)
          ENDSELECT
        ENDDO
      ENDDO
    CASE(BASIS_SIMPLEX_TYPE)
      DO coupledMeshIdx=1,interfaceMeshConnectivity%NUMBER_OF_COUPLED_MESHES
        coupledMesh=>interfaceMeshConnectivity%INTERFACE%COUPLED_MESHES(coupledMeshIdx)%ptr
        IF(.NOT.ASSOCIATED(coupledMesh)) THEN
          localError="The interface coupled mesh for coupled mesh index "// &
            & TRIM(NumberToVString(coupledMeshIdx,"*",err,error))//" is not associated."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        numberOfMeshXi=coupledMesh%NUMBER_OF_DIMENSIONS
        numberOfMeshXiCoordinates=numberOfMeshXi+1
        DO interfaceElementIdx=1,interfaceMeshConnectivity%NUMBER_OF_INTERFACE_ELEMENTS
          elementConnectivity=>interfaceMeshConnectivity%ELEMENT_CONNECTIVITY(interfaceElementIdx,coupledMeshIdx)
          zeroXiCoordIdx=0
          DO localNodeIdx=1,numberOfInterfaceElementNodes
            xiCoordinates(1:numberOfMeshXi,localNodeIdx)=elementConnectivity%xi(1:numberOfMeshXi,1,localNodeIdx)
            CALL Basis_XiToAreaCoordinates(xiCoordinates(1:numberOfMeshXi,localNodeIdx), &
              & areaCoordinates(1:numberOfMeshXiCoordinates,localNodeIdx),err,error,*999)
          ENDDO !localNodeIdx
          DO xiCoordIdx=1,numberOfMeshXiCoordinates
            IF(ALL(areaCoordinates(xiCoordIdx,1:numberOfInterfaceElementNodes)<ZERO_TOLERANCE)) THEN
              IF(zeroXiCoordIdx==0) THEN
                zeroXiCoordIdx=xiCoordIdx
              ELSE
                localError="Invalid xi connectivity. Xi directions "//TRIM(NumberToVString(xiCoordIdx,"*",err,error))// &
                  & " and "//TRIM(NumberToVString(zeroXiCoordIdx,"*",err,error))// &
                  & " both have a zero area coordinate all interface nodes for interface element number "// &
                  & TRIM(NumberToVString(interfaceElementIdx,"*",err,error))//" and coupled mesh number "// &
                  & TRIM(NumberToVString(coupledMeshIdx,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              ENDIF
            ENDIF
          ENDDO !xiCoordIdx
          IF(zeroXiCoordIdx==0) THEN
            localError="Invalid xi connectivity. There are no xi directions that have a zero area coordinate "// &
              & "across all interface nodes for interface element number "// &
              & TRIM(NumberToVString(interfaceElementIdx,"*",err,error))//" and coupled mesh number "// &
              & TRIM(NumberToVString(coupledMeshIdx,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          SELECTCASE(numberOfInterfaceMeshXi)
          CASE(1) !Lines
            elementConnectivity%CONNECTED_LINE=zeroXiCoordIdx
          CASE(2) !Faces
            elementConnectivity%CONNECTED_FACE=zeroXiCoordIdx
          CASE DEFAULT
            localError="The number of interface mesh dimension of "// &
              & TRIM(NUMBER_TO_VSTRING(numberOfInterfaceMeshXi,"*",err,error))//" is invalid"
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ENDDO !interfaceElementIdx
      ENDDO !coupledMeshIdx
    CASE DEFAULT
      localError="The interface basis type of "//TRIM(NumberToVString(interfaceBasis%type,"*",err,error))// &
        & " is invalid or not implemented."
      CALL FlagError(localError,err,error,*999)
    END SELECT
        
    EXITS("InterfaceMeshConnectivity_ConnectedLinesFacesCalculate")
    RETURN
999 ERRORS("InterfaceMeshConnectivity_ConnectedLinesFacesCalculate",err,error)
    EXITS("InterfaceMeshConnectivity_ConnectedLinesFacesCalculate")
    RETURN 1
  
  END SUBROUTINE InterfaceMeshConnectivity_ConnectedLinesFacesCalculate

  !
  !================================================================================================================================
  !

  !>Finalises the meshes connectivity and deallocates all memory
  SUBROUTINE INTERFACE_MESH_CONNECTIVITY_FINALISE(INTERFACE_MESH_CONNECTIVITY,err,error,*)

    !Argument variables
    TYPE(INTERFACE_MESH_CONNECTIVITY_TYPE) :: INTERFACE_MESH_CONNECTIVITY !<The interface mesh connectivity to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
     
    ENTERS("INTERFACE_MESH_CONNECTIVITY_FINALISE",err,error,*999)

    CALL InterfaceMeshConnectivity_ElementFinalise(INTERFACE_MESH_CONNECTIVITY,err,error,*999)
    NULLIFY(INTERFACE_MESH_CONNECTIVITY%INTERFACE)
    NULLIFY(INTERFACE_MESH_CONNECTIVITY%INTERFACE_MESH)
    NULLIFY(INTERFACE_MESH_CONNECTIVITY%BASIS)
    INTERFACE_MESH_CONNECTIVITY%NUMBER_OF_INTERFACE_ELEMENTS=0
    INTERFACE_MESH_CONNECTIVITY%NUMBER_OF_COUPLED_MESHES=0
       
    EXITS("INTERFACE_MESH_CONNECTIVITY_FINALISE")
    RETURN
999 ERRORSEXITS("INTERFACE_MESH_CONNECTIVITY_FINALISE",err,error)
    RETURN 1
    
  END SUBROUTINE INTERFACE_MESH_CONNECTIVITY_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the interface mesh connectivity.
  SUBROUTINE INTERFACE_MESH_CONNECTIVITY_INITIALISE(INTERFACE,MESH,err,error,*)

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to initialise the mesh connectivity for
    TYPE(MESH_TYPE), POINTER :: MESH
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("INTERFACE_MESH_CONNECTIVITY_INITIALISE",err,error,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(ASSOCIATED(INTERFACE%MESH_CONNECTIVITY)) THEN
        CALL FlagError("Interface mesh connectivity is already associated.",err,error,*999)
      ELSE
        ALLOCATE(INTERFACE%MESH_CONNECTIVITY,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate interface mesh connectivity.",err,error,*999)
        INTERFACE%MESH_CONNECTIVITY%INTERFACE=>INTERFACE
        INTERFACE%MESH_CONNECTIVITY%MESH_CONNECTIVITY_FINISHED=.FALSE.
        INTERFACE%MESH_CONNECTIVITY%INTERFACE_MESH=>MESH
        NULLIFY(INTERFACE%MESH_CONNECTIVITY%BASIS)
        CALL InterfaceMeshConnectivity_ElementInitialise(INTERFACE,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Interface is not associated.",err,error,*999)
    ENDIF
    
    EXITS("INTERFACE_MESH_CONNECTIVITY_INITIALISE")
    RETURN
999 ERRORSEXITS("INTERFACE_MESH_CONNECTIVITY_INITIALISE",err,error)
    RETURN 1
  END SUBROUTINE INTERFACE_MESH_CONNECTIVITY_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Calculate the coupled mesh elements that are connected to each interface element
  SUBROUTINE InterfacePointsConnectivity_CoupledElementsCalculate(interfacePointsConnectivity,coupledMeshIdx,err,error,*) 

    !Argument variables
    TYPE(InterfacePointsConnectivityType), POINTER :: interfacePointsConnectivity !<A pointer to the interface points connectivity to calculate coupled elements for
    INTEGER(INTG), INTENT(IN) :: coupledMeshIdx !<The coupled mesh index
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(LIST_TYPE), POINTER :: elementNumbersList
    INTEGER(INTG) :: elementIdx,dataPointIdx,globalDataPointNumber,globalElementNumber,numberOfElementDataPoints, &
      & numberOfCoupledElements,coupledElementIdx
    INTEGER(INTG), ALLOCATABLE :: elementNumbers(:)
  
    ENTERS("InterfacePointsConnectivity_CoupledElementsCalculate",err,error,*999)

    IF(ASSOCIATED(interfacePointsConnectivity)) THEN
      IF(ALLOCATED(interfacePointsConnectivity%coupledElements)) THEN
        interfacePointsConnectivity%maxNumberOfCoupledElements(coupledMeshIdx)=0; !Initialise the number of coupled mesh elements 
        DO elementIdx=1,SIZE(interfacePointsConnectivity%coupledElements,1) !Number of interface elements
          numberOfElementDataPoints=interfacePointsConnectivity%interfaceMesh%TOPOLOGY(1)%PTR%dataPoints% &
            & elementDataPoint(elementIdx)%numberOfProjectedData !Get the number of data points in interface mesh element
          !Set up list
          NULLIFY(elementNumbersList)
          CALL LIST_CREATE_START(elementNumbersList,err,error,*999)
          CALL LIST_DATA_TYPE_SET(elementNumbersList,LIST_INTG_TYPE,err,error,*999)
          CALL LIST_INITIAL_SIZE_SET(elementNumbersList,numberOfElementDataPoints,err,error,*999)
          CALL LIST_CREATE_FINISH(elementNumbersList,err,error,*999)
          DO dataPointIdx=1,numberOfElementDataPoints
            globalDataPointNumber=interfacePointsConnectivity%interfaceMesh%TOPOLOGY(1)%PTR%dataPoints% &
              & elementDataPoint(elementIdx)%dataIndices(dataPointIdx)%globalNumber
            globalElementNumber=interfacePointsConnectivity%pointsConnectivity(globalDataPointNumber,coupledMeshIdx)% &
              & coupledMeshElementNumber
            CALL LIST_ITEM_ADD(elementNumbersList,globalElementNumber,err,error,*999)
          ENDDO !dataPointIdx
          CALL LIST_REMOVE_DUPLICATES(elementNumbersList,err,error,*999)
          CALL LIST_DETACH_AND_DESTROY(elementNumbersList,numberOfCoupledElements,elementNumbers, &
            & err,error,*999)
          IF(ALLOCATED(interfacePointsConnectivity%coupledElements(elementIdx,coupledMeshIdx)%elementNumbers)) & 
            & DEALLOCATE(interfacePointsConnectivity%coupledElements(elementIdx,coupledMeshIdx)%elementNumbers) !for updating coupledElements after projection
          ALLOCATE(interfacePointsConnectivity%coupledElements(elementIdx,coupledMeshIdx)% &
            & elementNumbers(numberOfCoupledElements),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate coupled mesh element numbers.",err,error,*999)
          DO coupledElementIdx=1,numberOfCoupledElements
            interfacePointsConnectivity%coupledElements(elementIdx,coupledMeshIdx)%elementNumbers(coupledElementIdx)= &
              & elementNumbers(coupledElementIdx)
          ENDDO
          interfacePointsConnectivity%coupledElements(elementIdx,coupledMeshIdx)%numberOfCoupledElements=numberOfCoupledElements
          IF(ALLOCATED(elementNumbers)) DEALLOCATE(elementNumbers)
          IF(interfacePointsConnectivity%maxNumberOfCoupledElements(coupledMeshIdx)<numberOfCoupledElements) THEN
            interfacePointsConnectivity%maxNumberOfCoupledElements(coupledMeshIdx)=numberOfCoupledElements
          ENDIF
        ENDDO !elementIdx
      ELSE
        CALL FlagError("Interface points connectivity coupled elements is not allocated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Interface points connectivity is not associated.",err,error,*999)
    ENDIF
    
    EXITS("InterfacePointsConnectivity_CoupledElementsCalculate")
    RETURN
999 ERRORS("InterfacePointsConnectivity_CoupledElementsCalculate",err,error)
    EXITS("InterfacePointsConnectivity_CoupledElementsCalculate")
    RETURN 1
    
  END SUBROUTINE InterfacePointsConnectivity_CoupledElementsCalculate
  
  !
  !================================================================================================================================
  !

  !>Finalise the points connectivity coupled mesh elements 
  SUBROUTINE InterfacePointsConnectivity_CoupledElementsFinalise(interfacePointsConnectivity,err,error,*) 

    !Argument variables
    TYPE(InterfacePointsConnectivityType), POINTER :: interfacePointsConnectivity !<A pointer to the interface points connectivity whose coupled mesh elements is to be finalised
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementIdx,coupledMeshIdx
  
    ENTERS("InterfacePointsConnectivity_CoupledElementsFinalise",err,error,*999)

    IF(ASSOCIATED(interfacePointsConnectivity)) THEN
      IF(ALLOCATED(interfacePointsConnectivity%coupledElements)) THEN
        DO coupledMeshIdx=1,SIZE(interfacePointsConnectivity%coupledElements,2)
          DO elementIdx=1,SIZE(interfacePointsConnectivity%coupledElements,1)
            interfacePointsConnectivity%coupledElements(elementIdx,coupledMeshIdx)%numberOfCoupledElements=0
            IF(ALLOCATED(interfacePointsConnectivity%coupledElements(elementIdx,coupledMeshIdx)%elementNumbers)) &
              & DEALLOCATE(interfacePointsConnectivity%coupledElements(elementIdx,coupledMeshIdx)%elementNumbers)
          ENDDO
        ENDDO
        DEALLOCATE(interfacePointsConnectivity%coupledElements)
      END IF
      IF(ALLOCATED(interfacePointsConnectivity%maxNumberOfCoupledElements)) &
        & DEALLOCATE(interfacePointsConnectivity%maxNumberOfCoupledElements)
    ELSE
      CALL FlagError("Interface points connectivity is not associated.",err,error,*999)
    ENDIF
    
    EXITS("InterfacePointsConnectivity_CoupledElementsFinalise")
    RETURN
999 ERRORS("InterfacePointsConnectivity_CoupledElementsFinalise",err,error)
    EXITS("InterfacePointsConnectivity_CoupledElementsFinalise")
    RETURN 1
    
  END SUBROUTINE InterfacePointsConnectivity_CoupledElementsFinalise
  
  !
  !================================================================================================================================
  !

  !>Initialise the coupled mesh elements for points connectivity
  SUBROUTINE InterfacePointsConnectivity_CoupledElementsInitialise(interfacePointsConnectivity,err,error,*) 

    !Argument variables
    TYPE(InterfacePointsConnectivityType), POINTER :: interfacePointsConnectivity !<A pointer to the interface points connectivity to initialise coupled elements for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfInterfaceElements,numberOfCoupledMeshes,coupledMeshIdx,elementIdx
    INTEGER(INTG) :: dummyErr !<The error code
    TYPE(VARYING_STRING)  :: dummyError !<The error string
  
    ENTERS("InterfacePointsConnectivity_CoupledElementsInitialise",err,error,*999)

     IF(ASSOCIATED(interfacePointsConnectivity)) THEN
       IF(ALLOCATED(interfacePointsConnectivity%coupledElements)) THEN
         CALL FlagError("Interface points connectivity coupled elements is already allocated.",err,error,*999)
       ELSE
         IF(ASSOCIATED(interfacePointsConnectivity%interface)) THEN
           IF(ASSOCIATED(interfacePointsConnectivity%interfaceMesh)) THEN
             numberOfInterfaceElements=interfacePointsConnectivity%interfaceMesh%NUMBER_OF_ELEMENTS
             numberOfCoupledMeshes=interfacePointsConnectivity%interface%NUMBER_OF_COUPLED_MESHES
             ALLOCATE(interfacePointsConnectivity%coupledElements(numberOfInterfaceElements,numberOfCoupledMeshes),STAT=ERR)
             IF(ERR/=0) CALL FlagError("Could not allocate points connectivity coupled element.",err,error,*999)
             DO coupledMeshIdx=1,numberOfCoupledMeshes
               DO elementIdx=1,numberOfInterfaceElements
                 interfacePointsConnectivity%coupledElements(elementIdx,coupledMeshIdx)%numberOfCoupledElements=0
               ENDDO !elementIdx
             ENDDO !coupledMeshIdx
           ELSE
             CALL FlagError("Interface mesh is not associated.",err,error,*999)
           ENDIF
         ELSE
            CALL FlagError("Interface is not associated.",err,error,*999)
         ENDIF
       ENDIF
     ELSE
       CALL FlagError("Interface points connectivity is not associated.",err,error,*999)
     ENDIF
    
    EXITS("InterfacePointsConnectivity_CoupledElementsInitialise")
    RETURN
999 CALL InterfacePointsConnectivity_CoupledElementsFinalise(interfacePointsConnectivity,dummyErr,dummyError,*998) 
998 ERRORS("InterfacePointsConnectivity_CoupledElementsInitialise",err,error)
    EXITS("InterfacePointsConnectivity_CoupledElementsInitialise")
    RETURN 1
    
  END SUBROUTINE InterfacePointsConnectivity_CoupledElementsInitialise
  
  !
  !================================================================================================================================
  !

  !>Finish create interface points connectivity
  SUBROUTINE InterfacePointsConnectivity_CreateFinish(interfacePointsConnectivity,err,error,*) 

    !Argument variables
    TYPE(InterfacePointsConnectivityType), POINTER :: interfacePointsConnectivity !<A pointer to the interface points connectivity to finish creating
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(INTERFACE_TYPE), POINTER :: interface
    INTEGER(INTG) :: coupledMeshIdx
  
    ENTERS("InterfacePointsConnectivity_CreateFinish",err,error,*999)

    IF(ASSOCIATED(interfacePointsConnectivity)) THEN
      IF(interfacePointsConnectivity%pointsConnectivityFinished) THEN
        CALL FlagError("Interface points connectivity has already been finished.",err,error,*999)
      ELSE
        CALL InterfacePointsConnectivity_ReducedXiCalculate(interfacePointsConnectivity,err,error,*999) 
        interface=>interfacePointsConnectivity%interface
        IF(ASSOCIATED(interface)) THEN
          CALL InterfacePointsConnectivity_CoupledElementsInitialise(interfacePointsConnectivity,err,error,*999) 
          DO coupledMeshIdx=1,interfacePointsConnectivity%interface%NUMBER_OF_COUPLED_MESHES
            CALL InterfacePointsConnectivity_CoupledElementsCalculate(interfacePointsConnectivity,coupledMeshIdx,err,error,*999) 
          ENDDO  
        ELSE
          CALL FlagError("Interface is not associated.",err,error,*999)
        ENDIF
        interfacePointsConnectivity%pointsConnectivityFinished=.TRUE.   
      ENDIF
    ELSE
      CALL FlagError("Interface points connectivity is not associated.",err,error,*999)
    ENDIF
    
    EXITS("InterfacePointsConnectivity_CreateFinish")
    RETURN
999 ERRORSEXITS("InterfacePointsConnectivity_CreateFinish",err,error)
    RETURN 1
  END SUBROUTINE InterfacePointsConnectivity_CreateFinish

  !
  !================================================================================================================================
  !

  !>Start create interface points connectivity
  SUBROUTINE InterfacePointsConnectivity_CreateStart(INTERFACE,interfaceMesh,dataPoints,interfacePointsConnectivity,err,error,*) 

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: interface !<A pointer to the interface to create the points connectivity for
    TYPE(MESH_TYPE), POINTER :: interfaceMesh !<A pointer to the interface mesh for which the points connectivity is created
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points in the interface 
    TYPE(InterfacePointsConnectivityType), POINTER :: interfacePointsConnectivity !<On return, a pointer to the created points connectivity
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("InterfacePointsConnectivity_CreateStart",err,error,*999)

    IF(ASSOCIATED(interface)) THEN
      IF(interface%INTERFACE_FINISHED) THEN
        IF(ASSOCIATED(interface%pointsConnectivity)) THEN
          CALL FlagError("The interface already has a points connectivity associated.",err,error,*999)
        ELSE
          IF(ASSOCIATED(interfaceMesh)) THEN
            IF(ASSOCIATED(dataPoints)) THEN
              IF(dataPoints%dataPointsFinished) THEN
                IF(ASSOCIATED(INTERFACE,dataPoints%INTERFACE)) THEN
                  CALL InterfacePointsConnectivity_Initialise(INTERFACE,interfaceMesh,dataPoints,err,error,*999)
                  !Return the pointer
                  interfacePointsConnectivity=>interface%pointsConnectivity
                ELSE
                  CALL FlagError("The specified data points are not defined on the specified interface.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("The specified data points have not been finished.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Data points are not associated.",err,error,*999)  
            ENDIF
          ELSE
            CALL FlagError("Interface mesh is not associated.",err,error,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FlagError("Interface has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Interface is not associated.",err,error,*999)
    ENDIF
    
    EXITS("InterfacePointsConnectivity_CreateStart")
    RETURN
999 ERRORSEXITS("InterfacePointsConnectivity_CreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE InterfacePointsConnectivity_CreateStart
  
  !
  !================================================================================================================================
  !

  !>Reproject data points for points connectivity
  SUBROUTINE InterfacePointsConnectivity_DataReprojection(interface,interfaceCondition,err,error,*) 

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: interface !<A pointer to the interface where data reprojection is performed
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: interfaceCondition !<A pointer to the interface where data reprojection is performed
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(InterfacePointsConnectivityType), POINTER :: interfacePointsConnectivity
    TYPE(FIELD_TYPE), POINTER :: dependentFieldFixed,dependentFieldProjection
    TYPE(DataPointsType), POINTER :: dataPoints
    TYPE(DataProjectionType), POINTER :: dataProjection
    TYPE(FIELD_INTERPOLATED_POINT_PTR_TYPE), POINTER :: interpolatedPoints(:)
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: interpolatedPoint
    TYPE(FIELD_INTERPOLATION_PARAMETERS_PTR_TYPE), POINTER :: interpolationParameters(:)
    INTEGER(INTG) :: fixedBodyIdx,projectionBodyIdx,dataPointIdx
    INTEGER(INTG) :: elementNumber,numberOfGeometricComponents
    INTEGER(INTG) :: coupledMeshFaceLineNumber,component
  
    ENTERS("InterfacePointsConnectivity_DataReprojection",err,error,*999)
    
    NULLIFY(interpolatedPoints)
    NULLIFY(interpolationParameters)
    fixedBodyIdx=2 !\todo: need to generalise
    projectionBodyIdx=1

    IF(ASSOCIATED(interface)) THEN
      IF(ASSOCIATED(interfaceCondition)) THEN
        interfacePointsConnectivity=>interface%pointsConnectivity
        IF(ASSOCIATED(interfacePointsConnectivity)) THEN
          dataPoints=>interfacePointsConnectivity%dataPoints
          IF(ASSOCIATED(dataPoints)) THEN
            IF(ASSOCIATED(dataPoints%dataProjections)) THEN

              !Evaluate data points positions
              dependentFieldFixed=>interfaceCondition%DEPENDENT%FIELD_VARIABLES(fixedBodyIdx)%PTR%FIELD
              IF(ASSOCIATED(dependentFieldFixed)) THEN
                numberOfGeometricComponents=dependentFieldFixed%GEOMETRIC_FIELD%VARIABLES(1)%NUMBER_OF_COMPONENTS
                CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(dependentFieldFixed,interpolationParameters,err,error,*999, &
                  & FIELD_GEOMETRIC_COMPONENTS_TYPE)
                CALL FIELD_INTERPOLATED_POINTS_INITIALISE(interpolationParameters,interpolatedPoints,err,error,*999, &
                  & FIELD_GEOMETRIC_COMPONENTS_TYPE)
                interpolatedPoint=>interpolatedPoints(FIELD_U_VARIABLE_TYPE)%PTR
                DO dataPointIdx=1,dataPoints%numberOfDataPoints
                  elementNumber=interfacePointsConnectivity%pointsConnectivity(dataPointIdx,fixedBodyIdx)% &
                    & coupledMeshElementNumber
                  coupledMeshFaceLineNumber=dependentFieldFixed%DECOMPOSITION%TOPOLOGY%ELEMENTS% &
                    & ELEMENTS(elementNumber)% &
                    & ELEMENT_FACES(interfacePointsConnectivity%pointsConnectivity(dataPointIdx,fixedBodyIdx)% &
                    & elementLineFaceNumber)
                  CALL FIELD_INTERPOLATION_PARAMETERS_FACE_GET(FIELD_VALUES_SET_TYPE,coupledMeshFaceLineNumber, &
                    & interpolationParameters(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
                  CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,interfacePointsConnectivity%pointsConnectivity(dataPointIdx, &
                    & fixedBodyIdx)%reducedXi(:),interpolatedPoint,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE) !Interpolate contact data points on each surface
                  DO component=1,numberOfGeometricComponents
                    dataPoints%dataPoints(dataPointIdx)%position(component) = interpolatedPoint%VALUES(component,NO_PART_DERIV)
                  ENDDO !component
                ENDDO !dataPointIdx
              ELSE
                CALL FlagError("Fixed dependent field is not associated.",err,error,*999)
              ENDIF
            
              !Data reprojection and update points connectivity information with the projection results
              dataProjection=>dataPoints%dataProjections%dataProjections(projectionBodyIdx+1)%PTR 
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"ProjectedBodyDataProjectionLabel",err,error,*999)
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,dataProjection%label,err,error,*999)
              IF(ASSOCIATED(dataProjection)) THEN
                dependentFieldProjection=>interfaceCondition%DEPENDENT%FIELD_VARIABLES(projectionBodyIdx)%PTR%FIELD
                IF(ASSOCIATED(dependentFieldProjection)) THEN
                  !Projection the data points (with know spatial positions) on the projection dependent field 
                  CALL DataProjection_DataPointsProjectionEvaluate(dataProjection,FIELD_VALUES_SET_TYPE,err,error,*999)
                  CALL InterfacePointsConnectivity_UpdateFromProjection(InterfacePointsConnectivity,dataProjection, &
                    & projectionBodyIdx,err,error,*999) 
                ELSE
                  CALL FlagError("Projection dependent field is not associated.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Interface data projection is not associated.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Data points data projections is not associated.",err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Interface data points is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Interface points connectivity is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Interface condition is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Interface is not associated.",err,error,*999)
    ENDIF
    
    EXITS("InterfacePointsConnectivity_DataReprojection")
    RETURN
999 ERRORS("InterfacePointsConnectivity_DataReprojection",err,error)
    EXITS("InterfacePointsConnectivity_DataReprojection")
    RETURN 1
    
  END SUBROUTINE InterfacePointsConnectivity_DataReprojection
  
  !
  !================================================================================================================================
  !

  !>Destroy interface points connectivity
  SUBROUTINE InterfacePointsConnectivity_Destroy(interfacePointsConnectivity,err,error,*) 

    !Argument variables
    TYPE(InterfacePointsConnectivityType), POINTER :: interfacePointsConnectivity !<A pointer to interface points connectivity to be destroyed
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("InterfacePointsConnectivity_Destroy",err,error,*999)

    IF(ASSOCIATED(interfacePointsConnectivity)) THEN
      CALL InterfacePointsConnectivity_Finalise(interfacePointsConnectivity,err,error,*999) 
    ELSE
      CALL FlagError("Interface points connectivity is not associated.",err,error,*999)
    ENDIF
    
    EXITS("InterfacePointsConnectivity_Destroy")
    RETURN
999 ERRORSEXITS("InterfacePointsConnectivity_Destroy",err,error)
    RETURN 1
    
  END SUBROUTINE InterfacePointsConnectivity_Destroy
  
  !
  !================================================================================================================================
  !
  
  !>Gets the number of coupled mesh elements which are linked to a specific interface element.
  SUBROUTINE InterfacePointsConnectivity_ElementNumberGet(pointsConnectivity,dataPointUserNumber,coupledMeshIndexNumber, &
      & meshComponentNumber,coupledMeshUserElementNumber,err,error,*)
      
    !Argument variables
    TYPE(InterfacePointsConnectivityType), POINTER :: pointsConnectivity !<A pointer to interface points connectivity to set the element number of elements for.
    INTEGER(INTG), INTENT(IN) :: dataPointUserNumber !<The index of the data point.
    INTEGER(INTG), INTENT(IN) :: coupledMeshIndexNumber !<The index of the coupled mesh in the interface to set element number for
    INTEGER(INTG), INTENT(OUT) :: coupledMeshUserElementNumber !<The coupled mesh element user number
    INTEGER(INTG), INTENT(IN) :: meshComponentNumber !<The mesh component number of the interface mesh that points connectivity is associated to
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(MESH_TYPE), POINTER :: interfaceMesh
    INTEGER(INTG) :: dataPointGlobalNumber
    LOGICAL :: dataPointExists

    ENTERS("InterfacePointsConnectivity_ElementNumberGet",err,error,*999)
    
    IF(ASSOCIATED(pointsConnectivity)) THEN 
      CALL DataPoint_CheckExists(pointsConnectivity%dataPoints,dataPointUserNumber,dataPointExists, &
        & dataPointGlobalNumber,err,error,*999)
      IF(dataPointExists) THEN
        interfaceMesh=>pointsConnectivity%interfaceMesh
        coupledMeshUserElementNumber=pointsConnectivity%pointsConnectivity(dataPointGlobalNumber,coupledMeshIndexNumber)% &
          & coupledMeshElementNumber
      ELSE
        CALL FlagError("Data point with user number ("//TRIM(NUMBER_TO_VSTRING &
          & (dataPointUserNumber,"*",err,error))//") does not exist.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Interface points connectivity is not associated.",err,error,*999)
    ENDIF
    
    EXITS("InterfacePointsConnectivity_ElementNumberGet")
    RETURN
999 ERRORS("InterfacePointsConnectivity_ElementNumberGet",err,error)
    EXITS("InterfacePointsConnectivity_ElementNumberGet")
    RETURN 1
    
  END SUBROUTINE InterfacePointsConnectivity_ElementNumberGet
  
  !
  !================================================================================================================================
  !
  
  !>Sets the number of coupled mesh elements which are linked to a specific interface element.
  SUBROUTINE InterfacePointsConnectivity_ElementNumberSet(pointsConnectivity,dataPointUserNumber,coupledMeshIndexNumber, &
      & coupledMeshUserElementNumber,meshComponentNumber,err,error,*)
      
    !Argument variables
    TYPE(InterfacePointsConnectivityType), POINTER :: pointsConnectivity !<A pointer to interface points connectivity to set the element number of elements for.
    INTEGER(INTG), INTENT(IN) :: dataPointUserNumber !<The index of the data point.
    INTEGER(INTG), INTENT(IN) :: coupledMeshIndexNumber !<The index of the coupled mesh in the interface to set element number for
    INTEGER(INTG), INTENT(IN) :: coupledMeshUserElementNumber !<The coupled mesh element user number
    INTEGER(INTG), INTENT(IN) :: meshComponentNumber !<The mesh component number of the interface mesh that points connectivity is associated to
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointGlobalNumber,elementMeshNumber
    LOGICAL :: dataPointExists,elementExists

    ENTERS("InterfacePointsConnectivity_ElementNumberSet",err,error,*999)
    
    IF(ASSOCIATED(pointsConnectivity)) THEN
      IF(pointsConnectivity%pointsConnectivityFinished) THEN
        CALL FlagError("Interface points connectivity has already been finished.",err,error,*999)
      ELSE
        CALL DataPoint_CheckExists(pointsConnectivity%dataPoints,dataPointUserNumber,dataPointExists, &
          & dataPointGlobalNumber,err,error,*999)
        IF(dataPointExists) THEN
          IF ((coupledMeshIndexNumber<=pointsConnectivity%dataPoints%numberOfDataPoints).OR. &
              & (coupledMeshIndexNumber>0)) THEN
            IF (ALLOCATED(pointsConnectivity%pointsConnectivity)) THEN
              CALL MeshTopology_ElementCheckExists(pointsConnectivity%INTERFACE%COUPLED_MESHES(coupledMeshIndexNumber)%PTR, &
                & meshComponentNumber,coupledMeshUserElementNumber,elementExists,elementMeshNumber,err,error,*999) !Make sure user element exists       
              IF(elementExists) THEN
                pointsConnectivity%pointsConnectivity(dataPointGlobalNumber,coupledMeshIndexNumber)%coupledMeshElementNumber= &
                  & elementMeshNumber
              ELSE
                CALL FlagError("Element with user number ("//TRIM(NUMBER_TO_VSTRING &
                  & (coupledMeshUserElementNumber,"*",err,error))//") does not exist.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Interface points connectivity array not allocated.",err,error,*999)
            END IF
          ELSE
            CALL FlagError("Interface coupled mesh index number out of range.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Data point with user number ("//TRIM(NUMBER_TO_VSTRING &
              & (dataPointUserNumber,"*",err,error))//") does not exist.",err,error,*999)
        ENDIF
      ENDIF  
    ELSE
      CALL FlagError("Interface points connectivity is not associated.",err,error,*999)
    ENDIF
    
    EXITS("InterfacePointsConnectivity_ElementNumberSet")
    RETURN
999 ERRORS("InterfacePointsConnectivity_ElementNumberSet",err,error)
    EXITS("InterfacePointsConnectivity_ElementNumberSet")
    RETURN 1
    
  END SUBROUTINE InterfacePointsConnectivity_ElementNumberSet
  
  !
  !================================================================================================================================
  !

  !>Finalise interface points connectivity
  SUBROUTINE InterfacePointsConnectivity_Finalise(interfacePointsConnectivity,err,error,*) 

    !Argument variables
    TYPE(InterfacePointsConnectivityType), POINTER :: interfacePointsConnectivity !<A pointer to interface points connectivity to be finalised
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: coupledMeshIdx,dataPointIdx
    
    ENTERS("InterfacePointsConnectivity_Finalise",err,error,*999)
    
    IF(ASSOCIATED(interfacePointsConnectivity)) THEN
      IF(ALLOCATED(interfacePointsConnectivity%PointsConnectivity)) THEN
        DO coupledMeshIdx=1,size(interfacePointsConnectivity%pointsConnectivity,2)
          DO dataPointIdx=1,size(interfacePointsConnectivity%pointsConnectivity,1) !Deallocate memory for each data point
            CALL InterfacePointsConnectivity_PointFinalise(interfacePointsConnectivity%pointsConnectivity(dataPointIdx, &
              & coupledMeshIdx),err,error,*999) 
          ENDDO
        ENDDO
        DEALLOCATE(interfacePointsConnectivity%pointsConnectivity)  
      ENDIF
      CALL InterfacePointsConnectivity_CoupledElementsFinalise(interfacePointsConnectivity,err,error,*999) 
      NULLIFY(interfacePointsConnectivity%interface)
      NULLIFY(interfacePointsConnectivity%interfaceMesh) 
      IF(ALLOCATED(interfacePointsConnectivity%coupledElements)) DEALLOCATE(interfacePointsConnectivity%coupledElements)
      IF(ALLOCATED(interfacePointsConnectivity%maxNumberOfCoupledElements))  &
        & DEALLOCATE(interfacePointsConnectivity%maxNumberOfCoupledElements)
      DEALLOCATE(interfacePointsConnectivity)
    ENDIF
    
    EXITS("InterfacePointsConnectivity_Finalise")
    RETURN
999 ERRORSEXITS("InterfacePointsConnectivity_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE InterfacePointsConnectivity_Finalise
  
  !
  !================================================================================================================================
  !

  !>Finalise interface point connectivity
  SUBROUTINE InterfacePointsConnectivity_PointFinalise(interfacePointConnectivity,err,error,*) 

    !Argument variables
    TYPE(InterfacePointConnectivityType) :: interfacePointConnectivity !< An interface point connectivity to be finalised
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("InterfacePointsConnectivity_PointFinalise",err,error,*999)
    
    interfacePointConnectivity%coupledMeshElementNumber=0
    interfacePointConnectivity%elementLineFaceNumber=0
    IF(ALLOCATED(interfacePointConnectivity%xi)) DEALLOCATE(interfacePointConnectivity%xi)
    IF(ALLOCATED(interfacePointConnectivity%reducedXi)) DEALLOCATE(interfacePointConnectivity%reducedXi)
    
    EXITS("InterfacePointsConnectivity_PointFinalise")
    RETURN
999 ERRORSEXITS("InterfacePointsConnectivity_PointFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE InterfacePointsConnectivity_PointFinalise
  
  !
  !================================================================================================================================
  !

  !>Calculate full xi locations in points connectivity from reduced xi
  SUBROUTINE InterfacePointsConnectivity_FullXiCalculate(InterfacePointsConnectivity,coupledMeshIdx,err,error,*) 

    !Argument variables
    TYPE(InterfacePointsConnectivityType), POINTER :: InterfacePointsConnectivity !<A pointer to the interface points connectivity to finish creating
    INTEGER(INTG), INTENT(IN) :: coupledMeshIdx !<Coupled mesh index to calculate the full xi location for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointIdx,interfaceMeshDimensions,coupledMeshDimensions
    TYPE(INTERFACE_TYPE), POINTER :: interface
    TYPE(InterfacePointConnectivityType), POINTER :: pointConnectivity
    TYPE(MESH_TYPE), POINTER :: interfaceMesh,coupledMesh
    TYPE(VARYING_STRING) :: localError
  
    ENTERS("InterfacePointsConnectivity_FullXiCalculate",err,error,*999)

     IF(ASSOCIATED(InterfacePointsConnectivity)) THEN
       interfaceMesh=>InterfacePointsConnectivity%interfaceMesh
       IF(ASSOCIATED(interfaceMesh)) THEN
         interfaceMeshDimensions=InterfacePointsConnectivity%interfaceMesh%NUMBER_OF_DIMENSIONS
         interface=>InterfacePointsConnectivity%interface
         IF(ASSOCIATED(interface)) THEN
           coupledMesh=>interface%COUPLED_MESHES(coupledMeshIdx)%PTR
           IF(ASSOCIATED(coupledMesh)) THEN
             coupledMeshDimensions=coupledMesh%NUMBER_OF_DIMENSIONS
             IF(interfaceMeshDimensions==coupledMeshDimensions) THEN !e.g. If 1D-2D, 2D-3D coupling, interface dimension is 1D and 2D respectively for 1st body
               DO dataPointIdx=1,InterfacePointsConnectivity%dataPoints%numberOfDataPoints
                 InterfacePointsConnectivity%pointsConnectivity(dataPointIdx,coupledMeshIdx)%xi(:)= &
                   & InterfacePointsConnectivity%pointsConnectivity(dataPointIdx,coupledMeshIdx)%reducedXi(:)
               ENDDO !dataPointIdx
             ELSE
               !Update full xi location from reduced xi and element face/line number
               SELECT CASE(coupledMeshDimensions)
               CASE(2)
                 DO dataPointIdx=1,InterfacePointsConnectivity%dataPoints%numberOfDataPoints
                   pointConnectivity=>InterfacePointsConnectivity%pointsConnectivity(dataPointIdx,coupledMeshIdx)
                   SELECT CASE(pointConnectivity%elementLineFaceNumber)
                   CASE(1)
                     pointConnectivity%xi(1)=pointConnectivity%reducedXi(1)
                     pointConnectivity%xi(2)=0.0_DP
                   CASE(2)
                     pointConnectivity%xi(1)=pointConnectivity%reducedXi(1)
                     pointConnectivity%xi(2)=1.0_DP
                   CASE(3)
                     pointConnectivity%xi(1)=0.0_DP
                     pointConnectivity%xi(2)=pointConnectivity%reducedXi(1)
                   CASE(4)
                     pointConnectivity%xi(1)=1.0_DP
                     pointConnectivity%xi(2)=pointConnectivity%reducedXi(1)
                   CASE DEFAULT
                     localError="Invalid local line number "// &
                       & TRIM(NUMBER_TO_VSTRING(pointConnectivity%elementLineFaceNumber,"*",err,error))//" ."
                     CALL FlagError(localError,err,error,*999)
                   END SELECT
                 ENDDO !dataPointIdx
               CASE(3)
                 DO dataPointIdx=1,InterfacePointsConnectivity%dataPoints%numberOfDataPoints
                   pointConnectivity=>InterfacePointsConnectivity%pointsConnectivity(dataPointIdx,coupledMeshIdx)
                   SELECT CASE(pointConnectivity%elementLineFaceNumber)
                   CASE(1)
                     pointConnectivity%xi(1)=1.0_DP
                     pointConnectivity%xi(2)=pointConnectivity%reducedXi(1)
                     pointConnectivity%xi(3)=pointConnectivity%reducedXi(2)
                   CASE(2)
                     pointConnectivity%xi(1)=0.0_DP
                     pointConnectivity%xi(2)=pointConnectivity%reducedXi(1)
                     pointConnectivity%xi(3)=pointConnectivity%reducedXi(2)
                   CASE(3)
                     pointConnectivity%xi(1)=pointConnectivity%reducedXi(1)
                     pointConnectivity%xi(2)=1.0_DP
                     pointConnectivity%xi(3)=pointConnectivity%reducedXi(2)
                   CASE(4)
                     pointConnectivity%xi(1)=pointConnectivity%reducedXi(1)
                     pointConnectivity%xi(2)=0.0_DP
                     pointConnectivity%xi(3)=pointConnectivity%reducedXi(2)
                   CASE(5)
                     pointConnectivity%xi(1)=pointConnectivity%reducedXi(1)
                     pointConnectivity%xi(2)=pointConnectivity%reducedXi(2)
                     pointConnectivity%xi(3)=1.0_DP
                   CASE(6)
                     pointConnectivity%xi(1)=pointConnectivity%reducedXi(1)
                     pointConnectivity%xi(2)=pointConnectivity%reducedXi(2)
                     pointConnectivity%xi(3)=0.0_DP
                   CASE DEFAULT
                     localError="Invalid local face number "// &
                       & TRIM(NUMBER_TO_VSTRING(pointConnectivity%elementLineFaceNumber,"*",err,error))//" ."
                     CALL FlagError(localError,err,error,*999)
                   END SELECT
                 ENDDO !dataPointIdx
               CASE DEFAULT
                 localError="Invalid coupled mesh dimension "// &
                   & TRIM(NUMBER_TO_VSTRING(coupledMeshDimensions,"*",err,error))//" ."
                 CALL FlagError(localError,err,error,*999)
               END SELECT
             ENDIF
           ELSE
             CALL FlagError("Coupled mesh is not associated.",err,error,*999)
           ENDIF
         ELSE
           CALL FlagError("Interface is not associated.",err,error,*999)
         ENDIF
       ELSE
         CALL FlagError("Interface mesh is not associated.",err,error,*999)
       ENDIF
     ELSE
       CALL FlagError("Interface points connectivity is not associated.",err,error,*999)
     ENDIF
    
    EXITS("InterfacePointsConnectivity_FullXiCalculate")
    RETURN
999 ERRORSEXITS("InterfacePointsConnectivity_FullXiCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE InterfacePointsConnectivity_FullXiCalculate
  
    !
  !================================================================================================================================
  !

  !>Initialises the interface mesh connectivity.
  SUBROUTINE InterfacePointsConnectivity_Initialise(INTERFACE,interfaceMesh,dataPoints,err,error,*)

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: interface !<A pointer to the interface to initialise points connectivity for
    TYPE(MESH_TYPE), POINTER :: interfaceMesh !<A pointer to the interface mesh to initialise points connectivity with
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points to initialise the points connectivity with
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointIdx,meshIdx,coupledMeshDimension,interfaceMeshDimension
    INTEGER(INTG) :: dummyErr !<The error code
    TYPE(VARYING_STRING)  :: dummyError !<The error string
     
    ENTERS("InterfacePointsConnectivity_Initialise",err,error,*999)

    IF(ASSOCIATED(interface)) THEN
      IF(ASSOCIATED(interface%pointsConnectivity)) THEN
        CALL FlagError("Interface has already got points connectivity associated.",err,error,*998)
      ELSE
        IF(ASSOCIATED(interfaceMesh)) THEN
          IF(ASSOCIATED(dataPoints)) THEN
            !Initialise the poins connectivity
            ALLOCATE(INTERFACE%pointsConnectivity,STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate interface points connectivity.",err,error,*998)
            INTERFACE%pointsConnectivity%INTERFACE=>interface
            INTERFACE%pointsConnectivity%pointsConnectivityFinished=.FALSE.
            INTERFACE%pointsConnectivity%interfaceMesh=>interfaceMesh
            INTERFACE%pointsConnectivity%dataPoints=>dataPoints
            interfaceMeshDimension=interfaceMesh%NUMBER_OF_DIMENSIONS
            IF(interface%NUMBER_OF_COUPLED_MESHES>0) THEN 
              IF(dataPoints%numberOfDataPoints>0) THEN
                ALLOCATE(interface%pointsConnectivity%pointsConnectivity(dataPoints%numberOfDataPoints, &
                  & interface%NUMBER_OF_COUPLED_MESHES),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate interface point connectivity.",err,error,*999)
                DO dataPointIdx=1,dataPoints%numberOfDataPoints
                  DO meshIdx=1,interface%NUMBER_OF_COUPLED_MESHES
                    coupledMeshDimension=interface%COUPLED_MESHES(meshIdx)%PTR%NUMBER_OF_DIMENSIONS
                    CALL InterfacePointsConnectivity_PointInitialise(interface%pointsConnectivity%pointsConnectivity(dataPointIdx, &
                      & meshIdx),coupledMeshDimension,interfaceMeshDimension,err,error,*999)
                  ENDDO!meshIdx
                ENDDO!dataPointIdx
              ELSE
                CALL FlagError("Number of interface data points must be > 0.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Number of coupled meshes in the interface must be > 0.",err,error,*999)
            ENDIF 
            ALLOCATE(interface%pointsConnectivity%maxNumberOfCoupledElements(interface%NUMBER_OF_COUPLED_MESHES),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate interface max number of coupled mesh elements.",err,error,*999)
            INTERFACE%pointsConnectivity%maxNumberOfCoupledElements=0
          ELSE
            CALL FlagError("Data points are not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Interface mesh is not associated.",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Interface is not associated.",err,error,*999)
    ENDIF

    EXITS("InterfacePointsConnectivity_Initialise")
    RETURN
999 CALL InterfacePointsConnectivity_Finalise(interface%pointsConnectivity,dummyErr,dummyError,*998) 
998 ERRORSEXITS("InterfacePointsConnectivity_Initialise",err,error)
    RETURN 1
  END SUBROUTINE InterfacePointsConnectivity_Initialise
  
  !
  !================================================================================================================================
  !

  !>Initialises the interface mesh connectivity.
  SUBROUTINE InterfacePointsConnectivity_PointInitialise(interfacePointConnectivity,coupledMeshDimension,interfaceMeshDimension, &
      & err,error,*)

    !Argument variables
    TYPE(InterfacePointConnectivityType) :: interfacePointConnectivity !<An interface point connectivity to be initliased for
    INTEGER(INTG), INTENT(IN) :: coupledMeshDimension,interfaceMeshDimension
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr !<The error code
    TYPE(VARYING_STRING)  :: dummyError !<The error string
     
    ENTERS("InterfacePointsConnectivity_PointInitialise",err,error,*999)

    interfacePointConnectivity%coupledMeshElementNumber=0
    interfacePointConnectivity%elementLineFaceNumber=0
    !Allocate memory for coupled mesh full and reduced xi location
    ALLOCATE(interfacePointConnectivity%xi(coupledMeshDimension),STAT=ERR)
    IF(ERR/=0) CALL FlagError("Could not allocate interface point connectivity full xi.",err,error,*999)
    interfacePointConnectivity%xi=0.0_DP
    ALLOCATE(interfacePointConnectivity%reducedXi(interfaceMeshDimension),STAT=ERR)
    IF(ERR/=0) CALL FlagError("Could not allocate interface point connectivity reduced xi.",err,error,*999)
    interfacePointConnectivity%reducedXi=0.0_DP
       
    EXITS("InterfacePointsConnectivity_PointInitialise")
    RETURN
999 CALL InterfacePointsConnectivity_PointFinalise(interfacePointConnectivity,dummyErr,dummyError,*998)
998 ERRORSEXITS("InterfacePointsConnectivity_PointInitialise",err,error)
    RETURN 1
  END SUBROUTINE InterfacePointsConnectivity_PointInitialise
  
  !
  !================================================================================================================================
  !
    
  !>Gets the xi coordinate mapping between the data points in interface and xi coordinates in a coupled region mesh
  SUBROUTINE InterfacePointsConnectivity_PointXiGet(pointsConnectivity,dataPointUserNumber,coupledMeshIndexNumber, &
    & xi,err,error,*)

    !Argument variables
    TYPE(InterfacePointsConnectivityType), POINTER :: pointsConnectivity !<A pointer to interface points connectivity to set the element number of elements for.
    INTEGER(INTG), INTENT(IN) :: dataPointUserNumber !<The index of the data point.
    INTEGER(INTG), INTENT(IN) :: coupledMeshIndexNumber !<The index of the coupled mesh in the interface to set the number of elements for.
    REAL(DP), INTENT(OUT) :: xi(:) !<xi(xi_idx). The full xi location in the coupled mesh element.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointGlobalNumber
    LOGICAL :: dataPointExists
    
    ENTERS("InterfacePointsConnectivity_PointXiGet",err,error,*999)
    
    ! Preliminary error checks to verify user input information
    IF(ASSOCIATED(pointsConnectivity)) THEN
      IF (ALLOCATED(pointsConnectivity%pointsConnectivity)) THEN
        CALL DataPoint_CheckExists(pointsConnectivity%dataPoints,dataPointUserNumber,dataPointExists, &
          & dataPointGlobalNumber,err,error,*999)
        IF(dataPointExists) THEN
          IF(SIZE(xi)>=SIZE(pointsConnectivity%pointsConnectivity(dataPointUserNumber,coupledMeshIndexNumber)%xi,1)) THEN
            xi=pointsConnectivity%pointsConnectivity(dataPointUserNumber,coupledMeshIndexNumber)%xi(:)
          ELSE
            CALL FlagError("Input xi array size is smaller than points connectivity xi array.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Data point with user number ("//TRIM(NUMBER_TO_VSTRING &
              & (dataPointUserNumber,"*",err,error))//") does not exist.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Interface points connectivity array not allocated.",err,error,*999)
      ENDIF 
    ELSE
      CALL FlagError("Interface points connectivity is not associated.",err,error,*999)
    ENDIF   

    EXITS("InterfacePointsConnectivity_PointXiGet")
    RETURN
999 ERRORSEXITS("InterfacePointsConnectivity_PointXiGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfacePointsConnectivity_PointXiGet
  
  !
  !================================================================================================================================
  !
    
  !>Sets the xi coordinate mapping between the data points in interface and xi coordinates in a coupled region mesh
  SUBROUTINE InterfacePointsConnectivity_PointXiSet(pointsConnectivity,dataPointUserNumber,coupledMeshIndexNumber, &
    & xi,err,error,*)

    !Argument variables
    TYPE(InterfacePointsConnectivityType), POINTER :: pointsConnectivity !<A pointer to interface points connectivity to set the element number of elements for.
    INTEGER(INTG), INTENT(IN) :: dataPointUserNumber !<The index of the data point.
    INTEGER(INTG), INTENT(IN) :: coupledMeshIndexNumber !<The index of the coupled mesh in the interface to set the number of elements for.
    REAL(DP), INTENT(IN) :: xi(:) !<xi(xi_idx). The full xi location in the coupled mesh element.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointGlobalNumber
    LOGICAL :: dataPointExists
    
    ENTERS("InterfacePointsConnectivity_PointXiSet",err,error,*999)
    
    ! Preliminary error checks to verify user input information
    IF(ASSOCIATED(pointsConnectivity)) THEN
      IF(pointsConnectivity%pointsConnectivityFinished) THEN
        CALL FlagError("Interface mesh connectivity already been finished.",err,error,*999)
      ELSE
        IF (ALLOCATED(pointsConnectivity%pointsConnectivity)) THEN
          CALL DataPoint_CheckExists(pointsConnectivity%dataPoints,dataPointUserNumber,dataPointExists, &
            & dataPointGlobalNumber,err,error,*999)
          IF(dataPointExists) THEN
            IF ((coupledMeshIndexNumber>pointsConnectivity%interface%NUMBER_OF_COUPLED_MESHES).OR.(coupledMeshIndexNumber<0)) THEN
              CALL FlagError("Interface coupled mesh index number out of range.",err,error,*999)
            ELSE
              IF(SIZE(pointsConnectivity%pointsConnectivity(dataPointGlobalNumber,coupledMeshIndexNumber)%xi,1)== &
                  & SIZE(xi,1)) THEN
                pointsConnectivity%pointsConnectivity(dataPointGlobalNumber,coupledMeshIndexNumber)% &
                  & xi(:)=xi(:) 
              ELSE
                CALL FlagError("Input xi dimension does not match full coupled mesh xi dimension.",err,error,*999)
              ENDIF
            ENDIF
          ELSE
            CALL FlagError("Data point with user number ("//TRIM(NUMBER_TO_VSTRING &
                & (dataPointUserNumber,"*",err,error))//") does not exist.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Interface points connectivity array not allocated.",err,error,*999)
        ENDIF 
      ENDIF
    ELSE
      CALL FlagError("Interface points connectivity is not associated.",err,error,*999)
    ENDIF   

    EXITS("InterfacePointsConnectivity_PointXiSet")
    RETURN
999 ERRORSEXITS("InterfacePointsConnectivity_PointXiSet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfacePointsConnectivity_PointXiSet
  
  !
  !================================================================================================================================
  !
  
  !>Update points connectivity with projection results
  SUBROUTINE InterfacePointsConnectivity_UpdateFromProjection(InterfacePointsConnectivity,dataProjection, &
      & coupledMeshIndex,err,error,*) 
  
    !Argument variables
    TYPE(InterfacePointsConnectivityType), POINTER :: InterfacePointsConnectivity !<A pointer to the interface points connectivity to finish creating
    TYPE(DataProjectionType), POINTER :: dataProjection !<The data projection that points connectivity update with
    INTEGER(INTG), INTENT(IN) :: coupledMeshIndex !<The mesh index of the the points connectivity to be updated
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointIdx
    TYPE(DataProjectionResultType), POINTER :: dataProjectionResult
    
    ENTERS("InterfacePointsConnectivity_UpdateFromProjection",err,error,*999)
    
    IF(ASSOCIATED(InterfacePointsConnectivity)) THEN
      IF(ASSOCIATED(dataProjection)) THEN
        IF(dataProjection%dataProjectionFinished) THEN
          WRITE(*,*) "InterfacePointsConnectivity_UpdateFromProjection"
          DO dataPointIdx=1,SIZE(dataProjection%dataProjectionResults,1) !Update reduced xi location, projection element number and element face/line number with projection results
            dataProjectionResult=>dataProjection%dataProjectionResults(dataPointIdx)
            InterfacePointsConnectivity%pointsConnectivity(dataPointIdx,coupledMeshIndex)%reducedXi(:)= &
              & dataProjectionResult%xi
            InterfacePointsConnectivity%pointsConnectivity(dataPointIdx,coupledMeshIndex)%coupledMeshElementNumber= &
              & dataProjectionResult%elementNumber
            InterfacePointsConnectivity%pointsConnectivity(dataPointIdx,coupledMeshIndex)%elementLineFaceNumber= &
              & dataProjectionResult%elementLineFaceNumber
          ENDDO !dataPointIdx
          CALL InterfacePointsConnectivity_FullXiCalculate(InterfacePointsConnectivity,coupledMeshIndex, &
            & err,error,*999) 
          !Update points connectivity coupledElement information
          CALL InterfacePointsConnectivity_CoupledElementsCalculate(InterfacePointsConnectivity,coupledMeshIndex,err,error,*999) 
        ELSE
          CALL FlagError("Data projection is not finished.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Interface points connectivity is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Interface points connectivity is not associated.",err,error,*999)
    ENDIF
     
    EXITS("InterfacePointsConnectivity_UpdateFromProjection")
    RETURN
999 ERRORS("InterfacePointsConnectivity_UpdateFromProjection",err,error)
    EXITS("InterfacePointsConnectivity_UpdateFromProjection")
    RETURN 1
    
  END SUBROUTINE InterfacePointsConnectivity_UpdateFromProjection
  
  !
  !================================================================================================================================
  !

  !>Calculate reduced xi locations in points connectivity from full xi
  SUBROUTINE InterfacePointsConnectivity_ReducedXiCalculate(InterfacePointsConnectivity,err,error,*) 

    !Argument variables
    TYPE(InterfacePointsConnectivityType), POINTER :: InterfacePointsConnectivity !<A pointer to the interface points connectivity to finish creating
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: meshIdx,dataPointIdx,xiIdx,interfaceMeshDimensions,coupledMeshDimensions
    TYPE(INTERFACE_TYPE), POINTER :: interface
  
    ENTERS("InterfacePointsConnectivity_ReducedXiCalculate",err,error,*999)

     IF(ASSOCIATED(InterfacePointsConnectivity)) THEN
       IF(InterfacePointsConnectivity%pointsConnectivityFinished) THEN
         CALL FlagError("Interface points connectivity has already been finished.",err,error,*999)
       ELSE
         interface=>InterfacePointsConnectivity%interface
         interfaceMeshDimensions=InterfacePointsConnectivity%interfaceMesh%NUMBER_OF_DIMENSIONS
         IF(ASSOCIATED(interface)) THEN
           DO meshIdx=1,interface%NUMBER_OF_COUPLED_MESHES
             coupledMeshDimensions=interface%COUPLED_MESHES(meshIdx)%PTR%NUMBER_OF_DIMENSIONS
             IF(interfaceMeshDimensions==coupledMeshDimensions) THEN !e.g. If 1D-2D, 2D-3D coupling, interface dimension is 1D and 2D respectively for 1st body
               DO dataPointIdx=1,InterfacePointsConnectivity%dataPoints%numberOfDataPoints
                 InterfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%reducedXi(:)= &
                   & InterfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%xi(:)
                 InterfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%elementLineFaceNumber=1 !The only local face/line for the body with lower dimension 
               ENDDO !dataPointIdx
             ELSE
               SELECT CASE(coupledMeshDimensions)
               CASE(2)
                 DO dataPointIdx=1,InterfacePointsConnectivity%dataPoints%numberOfDataPoints
                   DO xiIdx=1,coupledMeshDimensions
                     IF(ABS(InterfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%xi(xiIdx)) &
                       & < ZERO_TOLERANCE) THEN
                       InterfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%elementLineFaceNumber=4-(xiIdx-1)*2 !Calculate line number
                       EXIT
                     ELSEIF(ABS(InterfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%xi(xiIdx)-1.0_DP) &
                       & < ZERO_TOLERANCE) THEN
                       InterfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%elementLineFaceNumber=3-(xiIdx-1)*2 !Calculate line number
                       EXIT
                     ENDIF
                   ENDDO !xiIdx
                   InterfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%reducedXi= &
                     & InterfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)% &
                     & xi(OTHER_XI_DIRECTIONS2(xiIdx))  !Populate reducedXi 
                 ENDDO !dataPointIdx
               CASE(3)
                 DO dataPointIdx=1,InterfacePointsConnectivity%dataPoints%numberOfDataPoints
                   DO xiIdx=1,coupledMeshDimensions
                     IF(ABS(InterfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%xi(xiIdx)) &
                       & < ZERO_TOLERANCE) THEN
                       InterfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%elementLineFaceNumber=(xiIdx-1)*2+2 !Calculate face number
                       EXIT
                     ELSE IF(ABS(InterfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%xi(xiIdx)-1.0_DP) &
                       & > ZERO_TOLERANCE) THEN
                       InterfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%elementLineFaceNumber=(xiIdx-1)*2+1 !Calculate face number
                       EXIT
                     ENDIF
                   ENDDO !xiIdx
                   SELECT CASE(xiIdx) !Populate reducedXi 
                   CASE(1)
                     InterfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%reducedXi(1)= &
                      & InterfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%xi(2)
                     InterfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%reducedXi(2)= &
                      & InterfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%xi(3)
                   CASE(2)
                     InterfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%reducedXi(1)= &
                      & InterfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%xi(1)
                     InterfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%reducedXi(2)= &
                      & InterfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%xi(3)
                   CASE(3)
                     InterfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%reducedXi(1)= &
                      & InterfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%xi(1)
                     InterfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%reducedXi(2)= &
                      & InterfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%xi(2)
                   END SELECT
                 ENDDO !dataPointIdx
               CASE DEFAULT
                 ! Do nothing
               END SELECT
             ENDIF
           ENDDO !meshIdx
         ELSE
           CALL FlagError("Interface is not associated.",err,error,*999)
         ENDIF
       ENDIF
     ELSE
       CALL FlagError("Interface points connectivity is not associated.",err,error,*999)
     ENDIF
    
    EXITS("InterfacePointsConnectivity_ReducedXiCalculate")
    RETURN
999 ERRORS("InterfacePointsConnectivity_ReducedXiCalculate",err,error)
    EXITS("InterfacePointsConnectivity_ReducedXiCalculate")
    RETURN 1
    
  END SUBROUTINE InterfacePointsConnectivity_ReducedXiCalculate

  !
  !================================================================================================================================
  !

  !>Finalises interfaces and deallocates all memory.
  SUBROUTINE INTERFACES_FINALISE(INTERFACES,err,error,*) 

    !Argument variables
    TYPE(INTERFACES_TYPE), POINTER :: INTERFACES !<A pointer to the interfaces to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE
     
    ENTERS("INTERFACES_FINALISE",err,error,*999)

    IF(ASSOCIATED(INTERFACES)) THEN
      DO WHILE(INTERFACES%NUMBER_OF_INTERFACES>0)
        INTERFACE=>INTERFACES%INTERFACES(1)%PTR
        CALL INTERFACE_DESTROY(INTERFACE,err,error,*999)
      ENDDO
      IF(ASSOCIATED(INTERFACES%INTERFACES)) DEALLOCATE(INTERFACES%INTERFACES)
      DEALLOCATE(INTERFACES)
    ENDIF
    
    EXITS("INTERFACES_FINALISE")
    RETURN
999 ERRORSEXITS("INTERFACES_FINALISE",err,error)
    RETURN 1
  END SUBROUTINE INTERFACES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises interfaces for a region.
  SUBROUTINE INTERFACES_INITIALISE(REGION,err,error,*) 

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to initialise the interfaces for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
     
    ENTERS("INTERFACES_INITIALISE",err,error,*998)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%INTERFACES)) THEN
        LOCAL_ERROR="Region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",err,error))// &
          & " already has interfaces associated."
        CALL FlagError(LOCAL_ERROR,err,error,*999)
      ELSE
        ALLOCATE(REGION%INTERFACES,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate region interfaces.",err,error,*999)
        REGION%INTERFACES%PARENT_REGION=>REGION
        REGION%INTERFACES%NUMBER_OF_INTERFACES=0
        NULLIFY(REGION%INTERFACES%INTERFACES)
      ENDIF
    ELSE
      CALL FlagError("Region is not associated.",err,error,*999)
    ENDIF
    
    EXITS("INTERFACES_INITIALISE")
    RETURN
999 CALL INTERFACES_FINALISE(REGION%INTERFACES,err,error,*998)
998 ERRORSEXITS("INTERFACES_INITIALISE",err,error)
    RETURN 1
  END SUBROUTINE INTERFACES_INITIALISE

  !
  !================================================================================================================================
  !

  !>Initialises the interface element connectivity.
  SUBROUTINE InterfaceMeshConnectivity_ElementInitialise(INTERFACE,err,error,*)

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to initialise the mesh connectivity for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: InterfaceElementIdx,CoupledMeshIdx
    TYPE(INTERFACE_ELEMENT_CONNECTIVITY_TYPE), POINTER :: ELEMENT_CONNECTIVITY
     
    ENTERS("InterfaceMeshConnectivity_ElementInitialise",err,error,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(ASSOCIATED(INTERFACE%MESH_CONNECTIVITY)) THEN
        IF(ALLOCATED(INTERFACE%MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY)) THEN
          CALL FlagError("Interface mesh element connectivity is already allocated.",err,error,*999)
        ELSE
          IF(INTERFACE%NUMBER_OF_COUPLED_MESHES>0) THEN
            IF(INTERFACE%MESH_CONNECTIVITY%INTERFACE_MESH%NUMBER_OF_ELEMENTS>0) THEN
              INTERFACE%MESH_CONNECTIVITY%NUMBER_OF_INTERFACE_ELEMENTS=INTERFACE%MESHES%MESHES(1)%PTR%NUMBER_OF_ELEMENTS
              INTERFACE%MESH_CONNECTIVITY%NUMBER_OF_COUPLED_MESHES=INTERFACE%NUMBER_OF_COUPLED_MESHES
              ALLOCATE(INTERFACE%MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY(INTERFACE%MESH_CONNECTIVITY%NUMBER_OF_INTERFACE_ELEMENTS, &
                & INTERFACE%MESH_CONNECTIVITY%NUMBER_OF_COUPLED_MESHES),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate interface element connectivity.",err,error,*999)
              DO InterfaceElementIdx=1,INTERFACE%MESH_CONNECTIVITY%NUMBER_OF_INTERFACE_ELEMENTS
                DO CoupledMeshIdx=1,INTERFACE%MESH_CONNECTIVITY%NUMBER_OF_COUPLED_MESHES
                  ELEMENT_CONNECTIVITY=>INTERFACE%MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY(InterfaceElementIdx,CoupledMeshIdx)
                  ELEMENT_CONNECTIVITY%COUPLED_MESH_ELEMENT_NUMBER=0
                  ELEMENT_CONNECTIVITY%CONNECTED_FACE=0
                  ELEMENT_CONNECTIVITY%CONNECTED_LINE=0
                  ! Note that ELEMENT_CONNECTIVITY%XI(NumberOfCoupledMeshXiDirections,1,NumberOfInterfaceElementNodes) is allocated after the basis for the mesh connectivity has been specified in INTERFACE_MESH_CONNECTIVITY_BASIS_SET where the number of NumberOfInterfaceElementNodes can be determined.
                  !\todo see corresponding todo in regarding updating the structure of ELEMENT_CONNECTIVITY%XI
                ENDDO !CoupledMeshIdx
              ENDDO !InterfaceElementIdx
            ELSE
              CALL FlagError("Interface coupled meshes are not associated.",err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Interface coupled meshes are not associated.",err,error,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FlagError("Interface mesh connectivity is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Interface is not associated.",err,error,*999)
    ENDIF
    
    EXITS("InterfaceMeshConnectivity_ElementInitialise")
    RETURN
999 ERRORSEXITS("InterfaceMeshConnectivity_ElementInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMeshConnectivity_ElementInitialise

  !
  !================================================================================================================================
  !

  !>Finalises an interface element connectivity and deallocates all memory.
  SUBROUTINE InterfaceMeshConnectivity_ElementFinalise(INTERFACE_MESH_CONNECTIVITY,err,error,*) 

    !Argument variables
    TYPE(INTERFACE_MESH_CONNECTIVITY_TYPE) :: INTERFACE_MESH_CONNECTIVITY !<The interface element connectivity to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: InterfaceElementIdx,CoupledMeshIdx
     
    ENTERS("InterfaceMeshConnectivity_ElementFinalise",err,error,*999)

    DO InterfaceElementIdx=1,INTERFACE_MESH_CONNECTIVITY%NUMBER_OF_INTERFACE_ELEMENTS  
      DO CoupledMeshIdx=1,INTERFACE_MESH_CONNECTIVITY%NUMBER_OF_COUPLED_MESHES
        IF(ALLOCATED(INTERFACE_MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY)) THEN
          INTERFACE_MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY(InterfaceElementIdx,CoupledMeshIdx)%COUPLED_MESH_ELEMENT_NUMBER=0
          IF(ALLOCATED(INTERFACE_MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY(InterfaceElementIdx,CoupledMeshIdx)%XI)) THEN
            DEALLOCATE(INTERFACE_MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY(InterfaceElementIdx,CoupledMeshIdx)%XI)
          ENDIF
        ELSE
          CALL FlagError("Interface mesh connectivity element connectivity is being deallocated before allocation.", &
            & err,error,*999)
        ENDIF
      ENDDO !InterfaceElementIdx
    ENDDO !CoupledMeshIdx

    DEALLOCATE(INTERFACE_MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY)
    
    EXITS("InterfaceMeshConnectivity_ElementFinalise")
    RETURN
999 ERRORSEXITS("InterfaceMeshConnectivity_ElementFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMeshConnectivity_ElementFinalise

  !
  !================================================================================================================================
  !

END MODULE INTERFACE_ROUTINES
