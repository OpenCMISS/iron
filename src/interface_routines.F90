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
!> Contributor(s): Chris Bradley, David Nordsletten, Thiranja Prasad Babarenda Gamage
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
MODULE InterfaceRoutines

  USE BaseRoutines
  USE BasisRoutines
  USE CoordinateSystemAccessRoutines
  USE DataPointRoutines
  USE DataPointAccessRoutines
  USE DataProjectionRoutines
  USE DataProjectionAccessRoutines
  USE DecompositionAccessRoutines
  USE FieldRoutines
  USE GENERATED_MESH_ROUTINES
  USE INPUT_OUTPUT
  USE InterfaceAccessRoutines
  USE INTERFACE_CONDITIONS_ROUTINES
  USE ISO_VARYING_STRING
  USE Kinds
  USE Lists
  USE MeshRoutines
  USE MeshAccessRoutines
  USE NodeRoutines
  USE Strings
  USE Types

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  INTERFACE Interface_LabelGet
    MODULE PROCEDURE Interface_LabelGetC
    MODULE PROCEDURE Interface_LabelGetVS
  END INTERFACE Interface_LabelGet
  
  INTERFACE Interface_LabelSet
    MODULE PROCEDURE Interface_LabelSetC
    MODULE PROCEDURE Interface_LabelSetVS
  END INTERFACE Interface_LabelSet

  PUBLIC Interface_MeshAdd

  PUBLIC Interface_CreateStart,Interface_CreateFinish

  PUBLIC Interface_CoordinateSystemSet

  PUBLIC Interface_Destroy

  PUBLIC Interface_LabelGet,Interface_LabelSet

  PUBLIC Interfaces_Finalise,Interfaces_Initialise

  PUBLIC InterfaceMeshConnectivity_BasisSet
  
  PUBLIC InterfaceMeshConnectivity_CreateStart, InterfaceMeshConnectivity_CreateFinish

  PUBLIC InterfaceMeshConnectivity_Destroy

  PUBLIC InterfaceMeshConnectivity_ElementXiSet,InterfaceMeshConnectivity_ElementNumberSet
  
  PUBLIC InterfaceMeshConnectivity_NodeNumbersSet

  PUBLIC InterfacePointsConnectivity_CreateStart,InterfacePointsConnectivity_CreateFinish
  
  PUBLIC InterfacePointsConnectivity_DataReprojection

  PUBLIC InterfacePointsConnectivity_Destroy
  
  PUBLIC InterfacePointsConnectivity_ElementNumberGet,InterfacePointsConnectivity_ElementNumberSet
  
  PUBLIC InterfacePointsConnectivity_PointXiGet,InterfacePointsConnectivity_PointXiSet
  
  PUBLIC InterfacePointsConnectivity_UpdateFromProjection
  
CONTAINS

  !
  !================================================================================================================================
  !

  SUBROUTINE Interface_MeshAdd(INTERFACE,mesh,meshIndex,err,error,*)   

    !Argument variables
    TYPE(InterfaceType), POINTER :: interface !<A pointer to the interface to add a mesh to
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh to add to the interface
    INTEGER(INTG), INTENT(OUT) :: meshIndex !<On return, the index of the added mesh in the list of meshes in the interface
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: meshIdx
    LOGICAL :: meshAlreadyCoupled
    TYPE(MeshType), POINTER :: coupledMesh
    TYPE(MeshPtrType), ALLOCATABLE :: newCoupledMeshes(:)
    TYPE(RegionType), POINTER :: coupledMeshRegion,meshRegion
    TYPE(VARYING_STRING) :: localError

    ENTERS("Interface_MeshAdd",err,error,*999)

    CALL Interface_AssertNotFinished(INTERFACE,err,error,*999)
    CALL Mesh_AssertIsFinished(mesh,err,error,*999)
    NULLIFY(meshRegion)
    CALL Mesh_RegionGet(mesh,meshRegion,err,error,*999)
    
    ALLOCATE(newCoupledMeshes(INTERFACE%numberOfCoupledMeshes+1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new coupled meshes.",err,error,*999)
    !Check that the mesh is not already in the list of meshes for the interface.
    IF(INTERFACE%numberOfCoupledMeshes>0) THEN
      IF(.NOT.ALLOCATED(INTERFACE%coupledMeshes)) CALL FlagError("Interface coupled meshes is not associated.",err,error,*999)
      meshAlreadyCoupled=.FALSE.
      DO meshIdx=1,INTERFACE%numberOfCoupledMeshes
        NULLIFY(coupledMesh)
        CALL Interface_CoupledMeshGet(INTERFACE,meshIdx,coupledMesh,err,error,*999)
        NULLIFY(coupledMeshRegion)
        CALL Mesh_RegionGet(coupledMesh,coupledMeshRegion,err,error,*999)
        IF(meshRegion%userNumber==coupledMeshRegion%userNumber) THEN
          IF(mesh%userNumber==coupledMesh%userNumber) THEN
            meshAlreadyCoupled=.TRUE.
            EXIT
          ENDIF
        ENDIF
        newCoupledMeshes(meshIdx)%ptr=>interface%coupledMeshes(meshIdx)%ptr
      ENDDO !meshIdx
      IF(meshAlreadyCoupled) THEN
        localError="The supplied mesh has already been added to the list of coupled meshes at mesh index "// &
          & TRIM(NumberToVString(meshIdx,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDIF
    !Add the mesh to the list of coupled meshes
    newCoupledMeshes(INTERFACE%numberOfCoupledMeshes+1)%ptr=>mesh
    CALL MOVE_ALLOC(newCoupledMeshes,INTERFACE%coupledMeshes)
    !Increment the number of coupled meshes and return the index
    INTERFACE%numberOfCoupledMeshes=INTERFACE%numberOfCoupledMeshes+1
    meshIndex=INTERFACE%numberOfCoupledMeshes
     
    EXITS("Interface_MeshAdd")
    RETURN
999 IF(ALLOCATED(newCoupledMeshes)) DEALLOCATE(newCoupledMeshes)
    ERRORSEXITS("Interface_MeshAdd",err,error)
    RETURN 1
    
  END SUBROUTINE Interface_MeshAdd

  !
  !================================================================================================================================
  !

  !>Finishes the creation of an interface. \see OPENCMISS::Iron::cmfe_InterfaceCreateFinish
  SUBROUTINE Interface_CreateFinish(INTERFACE,err,error,*) 

    !Argument variables
    TYPE(InterfaceType), POINTER :: INTERFACE !<A pointer to the interface to finish the creation of
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
     
    ENTERS("Interface_CreateFinish",err,error,*999)

    CALL Interface_AssertNotFinished(INTERFACE,err,error,*999)
    IF(INTERFACE%numberOfCoupledMeshes<2) THEN
      localError="Invalid mesh coupling. Only "//TRIM(NumberToVString(INTERFACE%numberOfCoupledMeshes,"*",err,error))// &
        & " have been coupled. The number of coupled meshes must be >= 2."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    INTERFACE%interfaceFinished=.TRUE.
    
    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Interface :",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  User number = ",INTERFACE%userNumber,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Global number = ",INTERFACE%globalNumber,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Label = ",INTERFACE%LABEL,err,error,*999)
      IF(ASSOCIATED(INTERFACE%interfaces)) THEN
        IF(ASSOCIATED(INTERFACE%interfaces%parentRegion)) THEN
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Parent region user number = ",INTERFACE%interfaces% &
            & parentRegion%userNumber,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Parent region label = ",INTERFACE%interfaces% &
            & parentRegion%LABEL,err,error,*999)        
        ELSE
          CALL FlagError("Interfaces parent region is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Interface interfaces is not associated.",err,error,*999)
      ENDIF
    ENDIF
    
    EXITS("Interface_CreateFinish")
    RETURN
999 ERRORSEXITS("Interface_CreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE Interface_CreateFinish

  !
  !================================================================================================================================
  !

  !>Starts the creation of an interface on a parent region. \see OPENCMISS::Iron::cmfe_InterfaceCreateStart
  SUBROUTINE Interface_CreateStart(userNumber,parentRegion,INTERFACE,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the interface to create
    TYPE(RegionType), POINTER :: parentRegion !<A pointer to the parent region to create the interface on.
    TYPE(InterfaceType), POINTER :: INTERFACE !<On exit, a pointer to the created interface. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: interfaceIdx
    TYPE(InterfaceType), POINTER :: newInterface
    TYPE(InterfacePtrType), ALLOCATABLE :: newInterfaces(:)
    TYPE(VARYING_STRING) :: localError,localString

    NULLIFY(newInterface)
    
    ENTERS("Interface_CreateStart",err,error,*998)

    IF(.NOT.ASSOCIATED(parentRegion)) CALL FlagError("Parent region is not associated.",err,error,*999)
    IF(ASSOCIATED(INTERFACE)) CALL FlagError("Interface is already associated.",err,error,*999)
    
    NULLIFY(INTERFACE)
    CALL Interface_UserNumberFind(userNumber,parentRegion,INTERFACE,err,error,*998)
    IF(ASSOCIATED(INTERFACE)) THEN
      localError="Interface number "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " has already been created on region number "//TRIM(NumberToVString(parentRegion%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*998)
    ENDIF        
    NULLIFY(INTERFACE)
    !Allocate and set default interface properties.
    CALL Interface_Initialise(newInterface,err,error,*999)
    newInterface%userNumber=userNumber
    newInterface%globalNumber=parentRegion%interfaces%numberOfInterfaces+1
    localString="Interface_"//NumberToVString(userNumber,"*",err,error)
    IF(err/=0) GOTO 999
    newInterface%label=CHAR(localString)
    newInterface%interfaces=>parentRegion%interfaces
    newInterface%parentRegion=>parentRegion
    !Add new interface into list of interfaces in the parent region
    ALLOCATE(newInterfaces(parentRegion%interfaces%numberOfInterfaces+1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new interfaces.",err,error,*999)
    DO interfaceIdx=1,parentRegion%interfaces%numberOfInterfaces
      newInterfaces(interfaceIdx)%ptr=>parentRegion%interfaces%interfaces(interfaceIdx)%ptr
    ENDDO !interfaceIdx
    newInterfaces(parentRegion%interfaces%numberOfInterfaces+1)%ptr=>newInterface
    CALL MOVE_ALLOC(newInterfaces,parentRegion%interfaces%interfaces)
    parentRegion%interfaces%numberOfInterfaces=parentRegion%interfaces%numberOfInterfaces+1
    INTERFACE=>newInterface
    
    EXITS("Interface_CreateStart")
    RETURN
999 IF(ALLOCATED(newInterfaces)) DEALLOCATE(newInterfaces)
    CALL Interface_Finalise(INTERFACE,err,error,*998)
998 ERRORSEXITS("Interface_CreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE Interface_CreateStart

  !
  !================================================================================================================================
  !

  !>Sets the coordinate system of an interface.  \see OPENCMISS::Iron::cmfe_Interface_CoordinateSystemSet
  SUBROUTINE Interface_CoordinateSystemSet(INTERFACE,coordinateSystem,err,error,*)

    !Argument variables
    TYPE(InterfaceType), POINTER :: INTERFACE !<A pointer to the interface to set the coordinate system for
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem !<The coordinate system to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string7
    !Local Variables

    ENTERS("Interface_CoordinateSystemSet",err,error,*999)

    CALL Interface_AssertNotFinished(INTERFACE,err,error,*999)
    CALL CoordinateSystem_AssertIsFinished(coordinateSystem,err,error,*999)
    
    INTERFACE%coordinateSystem=>coordinateSystem
     
    EXITS("Interface_CoordinateSystemSet")
    RETURN
999 ERRORSEXITS("Interface_CoordinateSystemSet",err,error)
    RETURN 1
    
  END SUBROUTINE Interface_CoordinateSystemSet
  
!
  !================================================================================================================================
  !

  !>Destroys an interface. \see OPENCMISS::Iron::cmfe_InterfaceDestroy
  SUBROUTINE Interface_Destroy(INTERFACE,err,error,*) 

    !Argument variables
    TYPE(InterfaceType), POINTER :: INTERFACE !<A pointer to the interface to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: interfaceIdx,interfacePosition
    TYPE(InterfacePtrType), ALLOCATABLE :: newInterfaces(:)
    TYPE(InterfacesType), POINTER :: interfaces
     
    ENTERS("Interface_Destroy",err,error,*999)

    IF(.NOT.ASSOCIATED(INTERFACE)) CALL FlagError("Interface is not associated.",err,error,*999)
    
    NULLIFY(interfaces)
    CALL Interface_InterfacesGet(INTERFACE,interfaces,err,error,*999)
    
    interfacePosition=INTERFACE%globalNumber
    !Destroy all the interface condition components
    CALL Interface_Finalise(INTERFACE,err,error,*999)
        
    !Remove the interface from the list of interfaces
    IF(interfaces%numberOfInterfaces>1) THEN
      ALLOCATE(newInterfaces(interfaces%numberOfInterfaces-1),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new interface conditions.",err,error,*999)
      DO interfaceIdx=1,interfaces%numberOfInterfaces
        IF(interfaceIdx<interfacePosition) THEN
          newInterfaces(interfaceIdx)%ptr=>interfaces%interfaces(interfaceIdx)%ptr
        ELSE IF(interfaceIdx>interfacePosition) THEN
          interfaces%interfaces(interfaceIdx)%ptr%globalNumber=interfaces%interfaces(interfaceIdx)%ptr%globalNumber-1
          newInterfaces(interfaceIdx-1)%ptr=>interfaces%interfaces(interfaceIdx)%ptr
        ENDIF
      ENDDO !interfaceIdx
      CALL MOVE_ALLOC(newInterfaces,interfaces%interfaces)
      interfaces%numberOfInterfaces=interfaces%numberOfInterfaces-1
    ELSE
      DEALLOCATE(interfaces%interfaces)
      interfaces%numberOfInterfaces=0
    ENDIF
     
    EXITS("Interface_Destroy")
    RETURN
999 ERRORSEXITS("Interface_Destroy",err,error)
    RETURN 1
    
  END SUBROUTINE Interface_Destroy

  !
  !================================================================================================================================
  !

  !>Finalises an interface and deallocates all memory.
  SUBROUTINE Interface_Finalise(INTERFACE,err,error,*) 

    !Argument variables
    TYPE(InterfaceType), POINTER :: INTERFACE !<A pointer to the interface to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
     
    ENTERS("Interface_Finalise",err,error,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(ALLOCATED(INTERFACE%coupledMeshes)) DEALLOCATE(INTERFACE%coupledMeshes)
      CALL Interface_DecompositionConnectivityFinalise(INTERFACE%decompositionConnectivity,err,error,*999)
      CALL Interface_MeshConnectivityFinalise(INTERFACE%meshConnectivity,err,error,*999)
      CALL Interface_PointsConnectivityFinalise(INTERFACE%pointsConnectivity,err,error,*999)
      IF(ASSOCIATED(INTERFACE%NODES)) CALL Nodes_Destroy(INTERFACE%NODES,err,error,*999)
      CALL GENERATED_MESHES_FINALISE(INTERFACE%generatedMeshes,err,error,*999)
      CALL Meshes_Finalise(INTERFACE%MESHES,err,error,*999)
      CALL Fields_Finalise(INTERFACE%FIELDS,err,error,*999)
      CALL INTERFACE_CONDITIONS_FINALISE(INTERFACE%interfaceConditions,err,error,*999)
      CALL DataPointSets_Finalise(interface%dataPointSets,err,error,*999)
      DEALLOCATE(INTERFACE)
    ENDIF
    
    EXITS("Interface_Finalise")
    RETURN
999 ERRORSEXITS("Interface_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE Interface_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises an interface.
  SUBROUTINE Interface_Initialise(interface,err,error,*) 

    !Argument variables
    TYPE(InterfaceType), POINTER :: INTERFACE !<A pointer to the interface to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
     
    ENTERS("Interface_Initialise",err,error,*998)

    IF(ASSOCIATED(INTERFACE)) CALL FlagError("Interface is already associated.",err,error,*998)
    
    ALLOCATE(INTERFACE,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate interface.",err,error,*999)
    INTERFACE%userNumber=0
    INTERFACE%globalNumber=0
    INTERFACE%interfaceFinished=.FALSE.
    INTERFACE%label=""
    NULLIFY(INTERFACE%interfaces)
    NULLIFY(INTERFACE%parentRegion)
    INTERFACE%numberOfCoupledMeshes=0
    NULLIFY(INTERFACE%decompositionConnectivity)
    NULLIFY(INTERFACE%meshConnectivity)
    NULLIFY(INTERFACE%pointsConnectivity)
    NULLIFY(INTERFACE%nodes)
    NULLIFY(INTERFACE%meshes)
    NULLIFY(INTERFACE%generatedMeshes)
    NULLIFY(INTERFACE%fields)
    NULLIFY(INTERFACE%interfaceConditions)
    NULLIFY(INTERFACE%coordinateSystem)
    NULLIFY(INTERFACE%dataPointSets)
    CALL DataPointSets_Initialise(INTERFACE,err,error,*999)
    CALL Meshes_Initialise(INTERFACE,err,error,*999)
    CALL GENERATED_MESHES_INITIALISE(INTERFACE,err,error,*999)
    CALL Fields_Initialise(INTERFACE,err,error,*999)
    CALL INTERFACE_CONDITIONS_INITIALISE(INTERFACE,err,error,*999)
    
    EXITS("Interface_Initialise")
    RETURN
999 CALL Interface_Finalise(INTERFACE,dummyErr,dummyError,*998)
998 ERRORSEXITS("Interface_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE Interface_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises the decomposition connectivity and deallocates all memory
  SUBROUTINE Interface_DecompositionConnectivityFinalise(interfaceDecompositionConnectivity,err,error,*)

    !Argument variables
    TYPE(InterfaceDecompositionConnectivityType), POINTER :: interfaceDecompositionConnectivity !<The interface decomposition connectivity to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: coupledDecompositionIdx,elementIdx
     
    ENTERS("Interface_DecompositionConnectivityFinalise",err,error,*999)

    IF(ASSOCIATED(interfaceDecompositionConnectivity)) THEN
      IF(ALLOCATED(interfaceDecompositionConnectivity%elementConnectivity)) THEN
        DO coupledDecompositionIdx=1,SIZE(interfaceDecompositionConnectivity%elementConnectivity,2)
          DO elementIdx=1,SIZE(interfaceDecompositionConnectivity%elementConnectivity,1)
            CALL InterfaceElementConnectivity_Finalise(interfaceDecompositionConnectivity% &
              & elementConnectivity(elementIdx,coupledDecompositionIdx),err,error,*999)
          ENDDO !elementIdx
        ENDDO !coupledDecompositionIdx
        DEALLOCATE(interfaceDecompositionConnectivity%elementConnectivity)
      ENDIF
      DEALLOCATE(interfaceDecompositionConnectivity)
    ENDIF
       
    EXITS("Interface_DecompositionConnectivityFinalise")
    RETURN
999 ERRORSEXITS("Interface_DecompositionConnectivityFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE Interface_DecompositionConnectivityFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the interface mesh connectivity.
  SUBROUTINE Interface_DecompositionConnectivityInitialise(interface,decomposition,err,error,*)

    !Argument variables
    TYPE(InterfaceType), POINTER :: interface !<A pointer to the interface to initialise the decomposition connectivity for
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer to the decomposition to initialise the decomposition connectivity for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("Interface_DecompositionConnectivityInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(interface)) CALL FlagError("Interface is not associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*998)
    IF(ASSOCIATED(interface%decompositionConnectivity)) &
      & CALL FlagError("Interface decomposition connectivity is already associated.",err,error,*998)
    
    ALLOCATE(interface%decompositionConnectivity,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate interface decomposition connectivity.",err,error,*999)
    interface%decompositionConnectivity%interface=>interface
    interface%decompositionConnectivity%interfaceDecomposition=>decomposition
    NULLIFY(interface%decompositionConnectivity%basis)
    INTERFACE%decompositionConnectivity%numberOfInterfaceElements=0
    INTERFACE%decompositionConnectivity%totalNumberOfInterfaceElements=0
    INTERFACE%decompositionConnectivity%numberOfGlobalInterfaceElements=0
    INTERFACE%decompositionConnectivity%numberOfCoupledDecompositions=0
     
    EXITS("Interface_DecompositionConnectivityInitialise")
    RETURN
999 CALL Interface_DecompositionConnectivityFinalise(interface%decompositionConnectivity,dummyErr,dummyError,*998)
998 ERRORS("Interface_DecompositionConnectivityInitialise",err,error)
    EXITS("Interface_DecompositionConnectivityInitialise")
    RETURN 1
    
  END SUBROUTINE Interface_DecompositionConnectivityInitialise
  
  !
  !================================================================================================================================
  !

  !>Returns the label of an interface for a character label. \see OPENCMISS::Iron::cmfe_InterfaceLabelGet
  SUBROUTINE Interface_LabelGetC(interface,label,err,error,*)

    !Argument variables
    TYPE(InterfaceType), POINTER :: interface !<A pointer to the interface to get the label for
    CHARACTER(LEN=*), INTENT(OUT) :: label !<On return the interface label.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: cLength,vsLength

    ENTERS("Interface_LabelGetC",err,error,*999)

    IF(.NOT.ASSOCIATED(interface)) CALL FlagError("Interface is not associated.",err,error,*999)
    
    cLength=LEN(label)
    vsLength=LEN_TRIM(interface%label)
    IF(cLength>vsLength) THEN
      label=CHAR(interface%label,vsLength)
    ELSE
      label=CHAR(interface%label,cLength)
    ENDIF
    
    EXITS("Interface_LabelGetC")
    RETURN
999 ERRORSEXITS("Interface_LabelGetC",err,error)
    RETURN 1
    
  END SUBROUTINE Interface_LabelGetC

   !
  !================================================================================================================================
  !

  !>Returns the label of an interface for a varying string label. \see OPENCMISS::Iron::cmfe_InterfaceLabelGet
  SUBROUTINE Interface_LabelGetVS(INTERFACE,label,err,error,*)

    !Argument variables
    TYPE(InterfaceType), POINTER :: INTERFACE !<A pointer to the interface to get the label for
    TYPE(VARYING_STRING), INTENT(OUT) :: label !<On return the interface label.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Interface_LabelGetVS",err,error,*999)

    IF(.NOT.ASSOCIATED(INTERFACE)) CALL FlagError("Interface is not associated.",err,error,*999)
    
    !\todo The following line crashes the AIX compiler unless it has a VAR_STR(CHAR()) around it
    label=VAR_STR(CHAR(INTERFACE%label))
   
    EXITS("Interface_LabelGetVS")
    RETURN
999 ERRORSEXITS("Interface_LabelGetVS",err,error)
    RETURN 1
    
  END SUBROUTINE Interface_LabelGetVS

  !
  !================================================================================================================================
  !

  !>Sets the label of an interface for a character label. \see OPENCMISS::Iron::cmfe_InterfaceLabelSet
  SUBROUTINE Interface_LabelSetC(INTERFACE,label,err,error,*)

    !Argument variables
    TYPE(InterfaceType), POINTER :: INTERFACE !<A pointer to the interface to set the label for 
    CHARACTER(LEN=*), INTENT(IN) :: label !<The label to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Interface_LabelSetC",err,error,*999)

    IF(.NOT.ASSOCIATED(INTERFACE)) CALL FlagError("Interface is not associated.",err,error,*999)
    
    INTERFACE%label=label
     
    EXITS("Interface_LabelSetC")
    RETURN
999 ERRORSEXITS("Interface_LabelSetC",err,error)
    RETURN 1
    
  END SUBROUTINE Interface_LabelSetC

  !
  !================================================================================================================================
  !

  !>Sets the label of an interface for a varying string label. \see OPENCMISS::Iron::cmfe_InterfaceLabelSet
  SUBROUTINE Interface_LabelSetVS(INTERFACE,label,err,error,*)

    !Argument variables
    TYPE(InterfaceType), POINTER :: INTERFACE !<A pointer to the interface to set the label for 
    TYPE(VARYING_STRING), INTENT(IN) :: label !<The label to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Interface_LabelSetVS",err,error,*999)

    IF(.NOT.ASSOCIATED(INTERFACE)) CALL FlagError("Interface is not associated.",err,error,*999)
    
    INTERFACE%label=label
    
    EXITS("Interface_LabelSetVS")
    RETURN
999 ERRORSEXITS("Interface_LabelSetVS",err,error)
    RETURN 1
    
  END SUBROUTINE Interface_LabelSetVS

  !
  !================================================================================================================================
  !

  !>Finalises a meshes connectivity for an interface.
  SUBROUTINE InterfaceMeshConnectivity_CreateFinish(interfaceMeshConnectivity,err,error,*) 

    !Argument variables
    TYPE(InterfaceMeshConnectivityType), POINTER :: interfaceMeshConnectivity !<A pointer to the interface meshes connectivity to finish creating
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: interfaceElementIdx,CoupledMeshIdx
    TYPE(VARYING_STRING) :: localError

    ENTERS("InterfaceMeshConnectivity_CreateFinish",err,error,*999)

    CALL InterfaceMeshConnectivity_AssertNotFinished(interfaceMeshConnectivity,err,error,*999)
    
    !Check if connectivity from each interface element to an appropriate element in the coupled meshes has been setup
    DO interfaceElementIdx=1,interfaceMeshConnectivity%numberOfInterfaceElements
      DO coupledMeshIdx=1,interfaceMeshConnectivity%numberOfCoupledMeshes
        IF (interfaceMeshConnectivity%elementConnectivity(interfaceElementIdx,coupledMeshIdx)%coupledElementNumber==0) THEN
          localError="The connectivity from interface element " &
            //TRIM(NumberToVString(interfaceElementIdx,"*",err,error))//" to an element in coupled mesh " & 
            //TRIM(NumberToVString(coupledMeshIdx,"*",err,error))//" has not been defined."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDDO !CoupledMeshIdx
    ENDDO !interfaceElementIdx
    !Calculate line or face numbers for coupled mesh elements that are connected to the interface mesh
    CALL InterfaceMeshConnectivity_ConnectedLinesFacesCalculate(interfaceMeshConnectivity,err,error,*999)
    interfaceMeshConnectivity%meshConnectivityFinished=.TRUE.

    EXITS("InterfaceMeshConnectivity_CreateFinish")
    RETURN
999 ERRORSEXITS("InterfaceMeshConnectivity_CreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMeshConnectivity_CreateFinish

  !
  !================================================================================================================================
  !

  !>Initialises a meshes connectivity for an interface.
  SUBROUTINE InterfaceMeshConnectivity_CreateStart(INTERFACE,mesh,interfaceMeshConnectivity,err,error,*) 

    !Argument variables
    TYPE(InterfaceType), POINTER :: INTERFACE !<A pointer to the interface to create the meshes connectivity for
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh to initialise the mesh connectivity for
    TYPE(InterfaceMeshConnectivityType), POINTER :: interfaceMeshConnectivity !<On return, a pointer to the created meshes connectivity
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("InterfaceMeshConnectivity_CreateStart",err,error,*999)

    CALL Interface_AssertIsFinished(INTERFACE,err,error,*999)
    IF(ASSOCIATED(INTERFACE%meshConnectivity)) &
      & CALL FlagError("The interface already has a meshes connectivity associated.",err,error,*999)
    
    !Initialise the meshes connectivity
    CALL Interface_MeshConnectivityInitialise(INTERFACE,MESH,err,error,*999)
    !Return the pointer
    interfaceMeshConnectivity=>INTERFACE%meshConnectivity
    
    EXITS("InterfaceMeshConnectivity_CreateStart")
    RETURN
999 ERRORSEXITS("InterfaceMeshConnectivity_CreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMeshConnectivity_CreateStart

  !
  !================================================================================================================================
  !

  !>Finalises a meshes connectivity and deallocates all memory
  SUBROUTINE InterfaceMeshConnectivity_Destroy(interfaceMeshConnectivity,err,error,*) 

    !Argument variables
    TYPE(InterfaceMeshConnectivityType), POINTER :: interfaceMeshConnectivity !<A pointer to the interface meshes connectivity to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("InterfaceMeshConnectivity_Destroy",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceMeshConnectivity)) CALL FlagError("Interface mesh connectivity is not associated.",err,error,*999)
    
    CALL Interface_MeshConnectivityFinalise(interfaceMeshConnectivity,err,error,*999)
       
    EXITS("InterfaceMeshConnectivity_Destroy")
    RETURN
999 ERRORSEXITS("InterfaceMeshConnectivity_Destroy",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMeshConnectivity_Destroy

  !
  !================================================================================================================================
  !

  !>Finalises the meshes connectivity and deallocates all memory
  SUBROUTINE Interface_MeshConnectivityFinalise(interfaceMeshConnectivity,err,error,*)

    !Argument variables
    TYPE(InterfaceMeshConnectivityType), POINTER :: interfaceMeshConnectivity !<A pointer to the interface mesh connectivity to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: coupledMeshIdx,elementIdx
     
    ENTERS("Interface_MeshConnectivityFinalise",err,error,*999)

    IF(ASSOCIATED(interfaceMeshConnectivity)) THEN
      IF(ALLOCATED(interfaceMeshConnectivity%elementConnectivity)) THEN
        DO coupledMeshIdx=1,SIZE(interfaceMeshConnectivity%elementConnectivity,2)
          DO elementIdx=1,SIZE(interfaceMeshConnectivity%elementConnectivity,1)
            CALL InterfaceElementConnectivity_Finalise(interfaceMeshConnectivity% &
              & elementConnectivity(elementIdx,coupledMeshIdx),err,error,*999)
          ENDDO !elementIdx
        ENDDO !coupledMeshIdx
        DEALLOCATE(interfaceMeshConnectivity%elementConnectivity)
      ENDIF
      IF(ALLOCATED(interfaceMeshConnectivity%coupledNodes)) DEALLOCATE(interfaceMeshConnectivity%coupledNodes)
      DEALLOCATE(interfaceMeshConnectivity)
    ENDIF
       
    EXITS("Interface_MeshConnectivityFinalise")
    RETURN
999 ERRORSEXITS("Interface_MeshConnectivityFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE Interface_MeshConnectivityFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the interface mesh connectivity.
  SUBROUTINE Interface_MeshConnectivityInitialise(INTERFACE,MESH,err,error,*)

    !Argument variables
    TYPE(InterfaceType), POINTER :: INTERFACE !<A pointer to the interface to initialise the mesh connectivity for
    TYPE(MeshType), POINTER :: MESH
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: coupledMeshIdx,dummyErr,interfaceElementIdx
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("Interface_MeshConnectivityInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(INTERFACE)) CALL FlagError("Interface is not associated.",err,error,*998)    
    IF(.NOT.ASSOCIATED(mesh)) CALL FlagError("Mesh is not associated.",err,error,*998)    
    IF(ASSOCIATED(INTERFACE%meshConnectivity)) &
      & CALL FlagError("Interface mesh connectivity is already associated.",err,error,*998)
    IF(mesh%numberOfElements<1) THEN
      localError="The number of elements in the specified mesh of "// &
        & TRIM(NumberToVString(mesh%numberOfElements,"*",err,error))// &
        & " is invalid. The number of elements should be >= 1."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(INTERFACE%numberOfCoupledMeshes<1) THEN
      localError="The number of coupled meshes of "//TRIM(NumberToVString(INTERFACE%numberOfCoupledMeshes,"*",err,error))// &
        & " for interface number "//TRIM(NumberToVString(INTERFACE%userNumber,"*",err,error))// &
        & " is invalid. The number of coupled meshes should be >= 1."
      CALL FlagError(localError,err,error,*999)
    ENDIF
      
    ALLOCATE(INTERFACE%meshConnectivity,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate interface mesh connectivity.",err,error,*999)
    INTERFACE%meshConnectivity%INTERFACE=>INTERFACE
    INTERFACE%meshConnectivity%meshConnectivityFinished=.FALSE.
    INTERFACE%meshConnectivity%interfaceMesh=>mesh
    NULLIFY(INTERFACE%meshConnectivity%basis)
    INTERFACE%meshConnectivity%numberOfInterfaceElements=mesh%numberOfElements
    INTERFACE%meshConnectivity%numberOfCoupledMeshes=INTERFACE%numberOfCoupledMeshes
    ALLOCATE(INTERFACE%meshConnectivity%elementConnectivity(mesh%numberOfElements,INTERFACE%numberOfCOupledMeshes),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate interface mesh connectivity element connectivity.",err,error,*999)
    DO coupledMeshIdx=1,INTERFACE%numberOfCoupledMeshes
      DO interfaceElementIdx=1,mesh%numberOfElements
        CALL InterfaceElementConnectivity_Initialise(INTERFACE%meshConnectivity%elementConnectivity(interfaceElementIdx, &
          & coupledMeshIdx),err,error,*999)
      ENDDO !interfaceElementIdx
    ENDDO !coupledMeshIdx
    
    EXITS("Interface_MeshConnectivityInitialise")
    RETURN
999 CALL Interface_MeshConnectivityFinalise(INTERFACE%meshConnectivity,dummyErr,dummyError,*998)
998 ERRORSEXITS("Interface_MeshConnectivityInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE Interface_MeshConnectivityInitialise
  
  !
  !================================================================================================================================
  !

  !>Reproject data points for points connectivity
  SUBROUTINE InterfacePointsConnectivity_DataReprojection(interface,interfaceCondition,err,error,*) 

    !Argument variables
    TYPE(InterfaceType), POINTER :: interface !<A pointer to the interface where data reprojection is performed
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface where data reprojection is performed
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: fixedBodyIdx,projectionBodyIdx,dataPointIdx
    INTEGER(INTG) :: elementNumber,numberOfGeometricComponents
    INTEGER(INTG) :: coupledMeshFaceLineNumber,component
    TYPE(DataPointsType), POINTER :: dataPoints
    TYPE(DataProjectionType), POINTER :: dataProjection
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(FieldType), POINTER :: dependentFieldFixed,dependentFieldProjection
    TYPE(FieldInterpolatedPointPtrType), POINTER :: interpolatedPoints(:)
    TYPE(FieldInterpolatedPointType), POINTER :: interpolatedPoint
    TYPE(FieldInterpolationParametersPtrType), POINTER :: interpolationParameters(:)
    TYPE(InterfacePointsConnectivityType), POINTER :: interfacePointsConnectivity
 
    ENTERS("InterfacePointsConnectivity_DataReprojection",err,error,*999)
    
    NULLIFY(interpolatedPoints)
    NULLIFY(interpolationParameters)
    fixedBodyIdx=2 !\todo: need to generalise
    projectionBodyIdx=1

    IF(.NOT.ASSOCIATED(interface)) CALL FlagError("Interface is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)

    NULLIFY(interfacePointsConnectivity)
    CALL Interface_PointsConnectivityGet(INTERFACE,interfacePointsConnectivity,err,error,*999)
    NULLIFY(dataPoints)
    CALL InterfacePointsConnectivity_DataPointsGet(interfacePointsConnectivity,dataPoints,err,error,*999)
    NULLIFY(dataProjection)
    CALL DataPoints_DataProjectionIndexGet(dataPoints,projectionBodyIdx+1,dataProjection,err,error,*999)

    !Evaluate data points positions
    dependentFieldFixed=>interfaceCondition%dependent%fieldVariables(fixedBodyIdx)%ptr%field
    NULLIFY(decomposition)
    CALL Field_DecompositionGet(dependentFieldFixed,decomposition,err,error,*999)
    NULLIFY(decompositionTopology)
    CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*999)
    NULLIFY(decompositionElements)
    CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology,decompositionElements,err,error,*999)
    IF(.NOT.ASSOCIATED(dependentFieldFixed)) CALL FlagError("Fixed dependent field is not associated.",err,error,*999)
    numberOfGeometricComponents=dependentFieldFixed%geometricField%VARIABLES(1)%numberOfComponents
    CALL Field_InterpolationParametersInitialise(dependentFieldFixed,interpolationParameters,err,error,*999, &
      & FIELD_GEOMETRIC_COMPONENTS_TYPE)
    CALL Field_InterpolatedPointsInitialise(interpolationParameters,interpolatedPoints,err,error,*999, &
      & FIELD_GEOMETRIC_COMPONENTS_TYPE)
    interpolatedPoint=>interpolatedPoints(FIELD_U_VARIABLE_TYPE)%ptr
    DO dataPointIdx=1,dataPoints%numberOfDataPoints
      elementNumber=interfacePointsConnectivity%pointsConnectivity(dataPointIdx,fixedBodyIdx)%coupledElementNumber
      coupledMeshFaceLineNumber=decompositionElements%elements(elementNumber)% &
        & elementFaces(interfacePointsConnectivity%pointsConnectivity(dataPointIdx,fixedBodyIdx)% &
        & elementLineFaceNumber)
      CALL Field_InterpolationParametersFaceGet(FIELD_VALUES_SET_TYPE,coupledMeshFaceLineNumber, &
        & interpolationParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
      CALL Field_InterpolateXi(NO_PART_DERIV,interfacePointsConnectivity%pointsConnectivity(dataPointIdx, &
        & fixedBodyIdx)%reducedXi(:),interpolatedPoint,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE) !Interpolate contact data points on each surface
      DO component=1,numberOfGeometricComponents
        dataPoints%dataPoints(dataPointIdx)%position(component) = interpolatedPoint%values(component,NO_PART_DERIV)
      ENDDO !component
    ENDDO !dataPointIdx
    
    !Data reprojection and update points connectivity information with the projection results
    dependentFieldProjection=>interfaceCondition%DEPENDENT%fieldVariables(projectionBodyIdx)%ptr%FIELD
    IF(ASSOCIATED(dependentFieldProjection)) THEN
      !Projection the data points (with know spatial positions) on the projection dependent field 
      CALL DataProjection_DataPointsProjectionEvaluate(dataProjection,FIELD_VALUES_SET_TYPE,err,error,*999)
      CALL InterfacePointsConnectivity_UpdateFromProjection(InterfacePointsConnectivity,dataProjection, &
        & projectionBodyIdx,err,error,*999) 
    ELSE
      CALL FlagError("Projection dependent field is not associated.",err,error,*999)
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

  !>Finalise the points connectivity coupled mesh elements 
  SUBROUTINE Interface_PointsConnectivityFinalise(interfacePointsConnectivity,err,error,*) 

    !Argument variables
    TYPE(InterfacePointsConnectivityType), POINTER :: interfacePointsConnectivity !<A pointer to the interface points connectivity whose coupled mesh elements is to be finalised
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointIdx,elementIdx,coupledMeshIdx
  
    ENTERS("Interface_PointsConnectivityFinalise",err,error,*999)

    IF(ASSOCIATED(interfacePointsConnectivity)) THEN
      IF(ALLOCATED(interfacePointsConnectivity%pointsConnectivity)) THEN
        DO coupledMeshIdx=1,SIZE(interfacePointsConnectivity%pointsConnectivity,2)
          DO dataPointidx=1,SIZE(interfacePointsConnectivity%pointsConnectivity,1)
            CALL InterfacePointConnectivity_Finalise(interfacePointsConnectivity% &
              & pointsConnectivity(dataPointIdx,coupledMeshIdx),err,error,*999)
          ENDDO !dataPointIdx
        ENDDO !coupledMeshIdx
        DEALLOCATE(interfacePointsConnectivity%pointsConnectivity)
      ENDIF
      IF(ALLOCATED(interfacePointsConnectivity%coupledElements)) THEN
        DO coupledMeshIdx=1,SIZE(interfacePointsConnectivity%coupledElements,2)
          DO elementIdx=1,SIZE(interfacePointsConnectivity%coupledElements,1)
            CALL InterfaceCoupledElements_Finalise(interfacePointsConnectivity% &
              & coupledElements(elementIdx,coupledMeshIdx),err,error,*999)
          ENDDO !elementIdx
        ENDDO !coupledMeshIdx
        DEALLOCATE(interfacePointsConnectivity%coupledElements)
      END IF
      IF(ALLOCATED(interfacePointsConnectivity%maxNumberOfCoupledElements)) &
        & DEALLOCATE(interfacePointsConnectivity%maxNumberOfCoupledElements)
      DEALLOCATE(interfacePointsConnectivity)
     ENDIF
    
    EXITS("Interface_PointsConnectivityFinalise")
    RETURN
999 ERRORS("Interface_PointsConnectivityFinalise",err,error)
    EXITS("Interface_PointsConnectivityFinalise")
    RETURN 1
    
  END SUBROUTINE Interface_PointsConnectivityFinalise
  
  !
  !================================================================================================================================
  !

  !>Initialises the interface mesh connectivity.
  SUBROUTINE Interface_PointsConnectivityInitialise(interface,interfaceMesh,dataPoints,err,error,*)

    !Argument variables
    TYPE(InterfaceType), POINTER :: interface !<A pointer to the interface to initialise points connectivity for
    TYPE(MeshType), POINTER :: interfaceMesh !<A pointer to the interface mesh to initialise points connectivity with
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points to initialise the points connectivity with
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: coupledMeshIdx,coupledMeshDimension,dataPointIdx,elementIdx,interfaceMeshDimension
    INTEGER(INTG) :: dummyErr 
    TYPE(MeshType), POINTER :: coupledMesh
    TYPE(VARYING_STRING)  :: dummyError,localError
     
    ENTERS("Interface_PointsConnectivityInitialise",err,error,*999)

    IF(.NOT.ASSOCIATED(INTERFACE)) CALL FlagError("Interface is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(interfaceMesh)) CALL FlagError("Interface mesh is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(dataPoints)) CALL FlagError("Data points is not associated.",err,error,*999)
    IF(ASSOCIATED(INTERFACE%pointsConnectivity)) &
      & CALL FlagError("Interface has already got points connectivity associated.",err,error,*998)
    IF(INTERFACE%numberOfCoupledMeshes<=0) THEN
      localError="The number of interface coupled meshes of "// &
        & TRIM(NumberToVString(INTERFACE%numberOfCoupledMeshes,"*",err,error))// &
        & " is invalid. The number of coupled meshes must be > 0."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(dataPoints%numberOfDataPoints<=0) THEN
      localError="The number of data points of "// &
        & TRIM(NumberToVString(dataPoints%numberOfDataPoints,"*",err,error))// &
        & " is invalid. The number of data points must be > 0."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    !Initialise the poins connectivity
    ALLOCATE(INTERFACE%pointsConnectivity,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate interface points connectivity.",err,error,*998)
    INTERFACE%pointsConnectivity%INTERFACE=>interface
    INTERFACE%pointsConnectivity%pointsConnectivityFinished=.FALSE.
    INTERFACE%pointsConnectivity%interfaceMesh=>interfaceMesh
    INTERFACE%pointsConnectivity%dataPoints=>dataPoints
    interfaceMeshDimension=interfaceMesh%numberOfDimensions
    ALLOCATE(INTERFACE%pointsConnectivity%pointsConnectivity(dataPoints%numberOfDataPoints,INTERFACE%numberOfCoupledMeshes), &
      & STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate interface point connectivity.",err,error,*999)
    DO coupledMeshIdx=1,INTERFACE%numberOfCoupledMeshes
      NULLIFY(coupledMesh)
      CALL Interface_CoupledMeshGet(INTERFACE,coupledMeshIdx,coupledMesh,err,error,*999)
      coupledMeshDimension=coupledMesh%numberOfDimensions
      DO dataPointIdx=1,dataPoints%numberOfDataPoints        
        CALL InterfacePointConnectivity_Initialise(INTERFACE%pointsConnectivity%pointsConnectivity(dataPointIdx, &
          & coupledMeshIdx),coupledMeshDimension,interfaceMeshDimension,err,error,*999)
      ENDDO!dataPointIdx
    ENDDO!meshIdx
    ALLOCATE(interface%pointsConnectivity%coupledElements(interfaceMesh%numberOfElements,interface%numberOfCoupledMeshes),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate points connectivity coupled element.",err,error,*999)
    DO coupledMeshIdx=1,interface%numberOfCoupledMeshes
      DO elementIdx=1,interfaceMesh%numberOfElements
        CALL InterfaceCoupledElements_Initialise(interface%pointsConnectivity%coupledElements(elementIdx,coupledMeshIdx), &
          & err,error,*999)
      ENDDO !elementIdx
    ENDDO !coupledMeshIdx
    ALLOCATE(INTERFACE%pointsConnectivity%maxNumberOfCoupledElements(INTERFACE%numberOfCoupledMeshes),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate interface max number of coupled mesh elements.",err,error,*999)
    INTERFACE%pointsConnectivity%maxNumberOfCoupledElements=0

    EXITS("Interface_PointsConnectivityInitialise")
    RETURN
999 CALL Interface_PointsConnectivityFinalise(interface%pointsConnectivity,dummyErr,dummyError,*998) 
998 ERRORSEXITS("Interface_PointsConnectivityInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE Interface_PointsConnectivityInitialise
  
  !
  !================================================================================================================================
  !

  !>Sets the interface mesh connectivity basis
  SUBROUTINE InterfaceMeshConnectivity_BasisSet(interfaceMeshConnectivity,basis,err,error,*)

    !Argument variables
    TYPE(InterfaceMeshConnectivityType), POINTER :: interfaceMeshConnectivity !<A pointer to interface mesh connectivity to set the element number of elements for.
    TYPE(BasisType), POINTER :: basis
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(InterfaceElementConnectivityType), POINTER :: elementConnectivity
    INTEGER(INTG) :: interfaceElementIdx,coupledMeshIdx,numberOfInterfaceElementNodes,numberOfCoupledMeshXiDirections

    ENTERS("InterfaceMeshConnectivity_BasisSet",err,error,*999)

    CALL InterfaceMeshConnectivity_AssertNotFinished(interfaceMeshConnectivity,err,error,*999)
    IF(.NOT.ASSOCIATED(basis)) CALL FlagError("Basis is not associated.",err,error,*999)
    IF(ASSOCIATED(interfaceMeshConnectivity%basis)) CALL FlagError("Mesh connectivity basis already associated.",err,error,*999)
    
    interfaceMeshConnectivity%basis=>basis
    !Now that the mesh connectivity basis is set the number of interface element nodes can be determined and now MESH_CONNECTIVITY%elementConnectivity(interfaceElementIdx,coupledMeshIdx)%XI can be allocated
    !\todo numberOfCoupledMeshXiDirections currently set to the number of interface mesh xi directions + 1. Restructure elementConnectivity type see below
    numberOfCoupledMeshXiDirections=interfaceMeshConnectivity%interfaceMesh%numberOfDimensions+1
    numberOfInterfaceElementNodes=interfaceMeshConnectivity%basis%numberOfNodes
    DO interfaceElementIdx=1,interfaceMeshConnectivity%numberOfInterfaceElements
      DO coupledMeshIdx=1,interfaceMeshConnectivity%numberOfCoupledMeshes
        elementConnectivity=>interfaceMeshConnectivity%elementConnectivity(interfaceElementIdx,coupledMeshIdx)
        !\todo Update mesh component index to look at the number of mesh components in each element. 
        !\todo Currently this defaults to the first mesh component ie %xi(NumberOfInterfaceMeshXi,1,numberOfInterfaceElementNodes)). 
        !\todo The interface mesh types will also need to be restructured.
        !eg. INTERFACE%meshConnectivity%elementConnectivity(interfaceElementIdx,coupledMeshIdx)%MESH_COMPONENT(MeshComponentIdx)%xi(numberOfCoupledMeshXiDirections,numberOfInterfaceElementNodes) and adding appropriate initialize/finialize routines
        ALLOCATE(elementConnectivity%xi(numberOfCoupledMeshXiDirections,1,numberOfInterfaceElementNodes),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate interface element connectivity.",err,error,*999)
        elementConnectivity%xi=0.0_DP
      ENDDO !coupledMeshIdx
    ENDDO !interfaceElementIdx

    EXITS("InterfaceMeshConnectivity_BasisSet")
    RETURN
999 ERRORSEXITS("InterfaceMeshConnectivity_BasisSet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMeshConnectivity_BasisSet

  !
  !================================================================================================================================
  !
  
  !>Calculate line or face numbers for coupled mesh elements that are connected to the interface mesh
  SUBROUTINE InterfaceMeshConnectivity_ConnectedLinesFacesCalculate(interfaceMeshConnectivity,err,error,*)

    !Argument variables
    TYPE(InterfaceMeshConnectivityType), POINTER :: interfaceMeshConnectivity !<A pointer to interface mesh connectivity to calculate line or face numbers for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    REAL(DP) :: areaCoordinates(4,4),xi,xiCoordinates(4,4),xiDifference
    INTEGER(INTG) :: constantXiIdx,coupledMeshIdx,coupledMeshXiIdx,interfaceElementIdx,localNodeIdx,numberOfInterfaceElementNodes, &
      & numberOfInterfaceMeshXi,numberOfMeshXi,numberOfMeshXiCoordinates,xiCoordIdx,zeroXiCoordIdx
    LOGICAL :: differentXi
    TYPE(BasisType), POINTER :: interfaceBasis
    TYPE(InterfaceType), POINTER :: interface
    TYPE(InterfaceElementConnectivityType), POINTER :: elementConnectivity
    TYPE(MeshType), POINTER :: coupledMesh, interfaceMesh
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("InterfaceMeshConnectivity_ConnectedLinesFacesCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceMeshConnectivity)) CALL FlagError("Interface mesh connectivity is not associated.",err,error,*999)

    NULLIFY(interfaceBasis)
    CALL InterfaceMeshConnectivity_BasisGet(interfaceMeshConnectivity,interfaceBasis,err,error,*999)
    NULLIFY(interfaceMesh)
    CALL InterfaceMeshConnectivity_InterfaceMeshGet(interfaceMeshConnectivity,interfaceMesh,err,error,*999)
    NULLIFY(INTERFACE)
    CALL InterfaceMeshConnectivity_InterfaceGet(interfaceMeshConnectivity,INTERFACE,err,error,*999)
    
    numberOfInterfaceElementNodes=interfaceBasis%numberOfNodes
    numberOfInterfaceMeshXi=interfaceMesh%numberOfDimensions     
    
!!\TODO: We need to check the basis in order to interpret the xi and decided to use area coordinates of normal coordinates. This really should be via the coupled mesh element being considered. For know assume either a LHTP basis or a simplex basis based on the interface basis. 
    SELECT CASE(interfaceBasis%type)
    CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
      DO coupledMeshIdx=1,interfaceMeshConnectivity%numberOfCoupledMeshes
        NULLIFY(coupledMesh)
        CALL Interface_CoupledMeshGet(INTERFACE,coupledMeshIdx,coupledMesh,err,error,*999)
        numberOfMeshXi=coupledMesh%numberOfDimensions
        DO interfaceElementIdx=1,interfaceMeshConnectivity%numberOfInterfaceElements
          elementConnectivity=>interfaceMeshConnectivity%elementConnectivity(interfaceElementIdx,coupledMeshIdx)
          SELECTCASE(numberOfInterfaceMeshXi)
          CASE(1) !Lines
            DO coupledMeshXiIdx=1,numberOfMeshXi
              ! Calculate difference between first node and last node of an element
              xiDifference=elementConnectivity%xi(coupledMeshXiIdx,1,1)- &
                & elementConnectivity%xi(coupledMeshXiIdx,1,numberOfInterfaceElementNodes)              
              IF(ABS(xiDifference)<ZERO_TOLERANCE) THEN
                IF(ABS(elementConnectivity%xi(coupledMeshXiIdx,1,numberOfInterfaceElementNodes))<ZERO_TOLERANCE) THEN
                  elementConnectivity%connectedLineFace=3-(coupledMeshXiIdx-1)*2
                ELSE
                  elementConnectivity%connectedLineFace=4-(coupledMeshXiIdx-1)*2
                ENDIF
              ENDIF
            ENDDO !coupledMeshXiIdx
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
                elementConnectivity%connectedLineFace=constantXiIdx
              ELSE
                elementConnectivity%connectedLineFace=constantXiIdx+3
              ENDIF
            ENDIF
          CASE DEFAULT 
            localError="The number of interface mesh dimension of "// &
              & TRIM(NumberToVString(numberOfInterfaceMeshXi,"*",err,error))//" is invalid"
            CALL FlagError(localError,err,error,*999)
          ENDSELECT
        ENDDO !interfaceElementIdx
      ENDDO !coupledMeshIdx
    CASE(BASIS_SIMPLEX_TYPE)
      DO coupledMeshIdx=1,interfaceMeshConnectivity%numberOfCoupledMeshes
        NULLIFY(coupledMesh)
        CALL Interface_CoupledMeshGet(INTERFACE,coupledMeshIdx,coupledMesh,err,error,*999)
        numberOfMeshXi=coupledMesh%numberOfDimensions
        numberOfMeshXiCoordinates=numberOfMeshXi+1
        DO interfaceElementIdx=1,interfaceMeshConnectivity%numberOfInterfaceElements
          elementConnectivity=>interfaceMeshConnectivity%elementConnectivity(interfaceElementIdx,coupledMeshIdx)
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
            elementConnectivity%connectedLineFace=zeroXiCoordIdx
          CASE(2) !Faces
            elementConnectivity%connectedLineFace=zeroXiCoordIdx
          CASE DEFAULT
            localError="The number of interface mesh dimension of "// &
              & TRIM(NumberToVString(numberOfInterfaceMeshXi,"*",err,error))//" is invalid"
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
    
  !>Sets the mapping from an xi position of a coupled mesh element to a node of an interface mesh element
  SUBROUTINE InterfaceMeshConnectivity_ElementXiSet(interfaceMeshConnectivity,interfaceMeshElementNumber, &
    & coupledMeshIndex,coupledElementNumber,interfaceMeshLocalNodeNumber,interfaceMeshComponentNumber,xi, &
    & err,error,*)

    !Argument variables
    TYPE(InterfaceMeshConnectivityType), POINTER :: interfaceMeshConnectivity !<A pointer to the interface mesh connectivity for the interface mesh
    INTEGER(INTG), INTENT(IN) :: interfaceMeshElementNumber !<The interface mesh element number to which the specified coupled mesh element would be connected
    INTEGER(INTG), INTENT(IN) :: coupledMeshIndex !<The index of the coupled mesh at the interface to set the element connectivity for
    INTEGER(INTG), INTENT(IN) :: coupledElementNumber !<The coupled mesh element to define the element xi connectivity from
    INTEGER(INTG), INTENT(IN) :: interfaceMeshLocalNodeNumber !<The interface mesh node to assign the coupled mesh element xi to
    INTEGER(INTG), INTENT(IN) :: interfaceMeshComponentNumber !<The interface mesh node's component to assign the coupled mesh element xi to
    REAL(DP), INTENT(IN) :: xi(:) !<xi(xi_idx). The xi value for the xi_idx'th xi direction in the coupled mesh element.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: coupledMeshNumberOfXi
    TYPE(BasisType), POINTER :: basis
    TYPE(InterfaceType), POINTER :: interface
    TYPE(InterfaceElementConnectivityType), POINTER :: elementConnectivity
    TYPE(MeshType), POINTER :: coupledMesh,interfaceMesh
    TYPE(MeshElementsType), POINTER :: coupledMeshElements
    TYPE(MeshTopologyType), POINTER :: coupledMeshTopology
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("InterfaceMeshConnectivity_ElementXiSet",err,error,*999)

    CALL InterfaceMeshConnectivity_AssertNotFinished(interfaceMeshConnectivity,err,error,*999)
    IF(.NOT.ALLOCATED(interfaceMeshConnectivity%elementConnectivity)) &
      & CALL FlagError("Interface elements connectivity array not allocated.",err,error,*999)
    IF((interfaceMeshElementNumber<=0).OR. &
      & (interfaceMeshElementNumber>interfaceMeshConnectivity%numberOfInterfaceElements)) THEN
      localError="The specified interface mesh element number of "// &
        & TRIM(NumberToVString(interfaceMeshElementNumber,"*",err,error))// &
        & " is invalid. The element number should be >= 1 and <= "// &
        & TRIM(NumberToVString(interfaceMeshConnectivity%numberOfInterfaceElements,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF    
    NULLIFY(INTERFACE)
    CALL InterfaceMeshConnectivity_InterfaceGet(interfaceMeshConnectivity,interface,err,error,*999)
    NULLIFY(coupledMesh)
    CALL Interface_CoupledMeshGet(interface,coupledMeshIndex,coupledMesh,err,error,*999)
    IF((coupledElementNumber<=0).OR.(coupledElementNumber>coupledMesh%numberOfElements)) THEN
      localError="The specified coupled element number of "// &
        & TRIM(NumberToVString(coupledElementNumber,"*",err,error))// &
        & " is invalid. The element number should be >= 1 and <= "// &
        & TRIM(NumberToVString(coupledMesh%numberOfElements,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    NULLIFY(interfaceMesh)
    CALL InterfaceMeshConnectivity_InterfaceMeshGet(interfaceMeshConnectivity,interfaceMesh,err,error,*999)
    IF((interfaceMeshComponentNumber<=0).OR.(interfaceMeshComponentNumber>interfaceMesh%numberOfComponents)) THEN
      localError="The specified interface mesh component number of "// &
        & TRIM(NumberToVString(interfaceMeshComponentNumber,"*",err,error))// &
        & " is invalid. The component number should be >= 1 and <= "// &
        & TRIM(NumberToVString(interfaceMesh%numberOfComponents,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)      
    ENDIF
    NULLIFY(basis)
    CALL InterfaceMeshConnectivity_BasisGet(interfaceMeshConnectivity,basis,err,error,*999)
    IF((interfaceMeshLocalNodeNumber<=0).OR.(interfaceMeshLocalNodeNumber>basis%numberOfNodes))THEN
      localError="The specified interface mesh local node number of "// &
        & TRIM(NumberToVString(interfaceMeshLocalNodeNumber,"*",err,error))// &
        & " is invalid. The local node number should be >= 1 and <= "// &
        & TRIM(NumberToVString(basis%numberOfNodes,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)      
    ENDIF
    elementConnectivity=>interfaceMeshConnectivity%elementConnectivity(interfaceMeshElementNumber,coupledMeshIndex)
    IF(elementConnectivity%coupledElementNumber/=coupledElementNumber) THEN
      localError="The specified coupled element number of "//TRIM(NumberToVString(coupledElementNumber,"*",err,error))// &
        & " does not match the element connectivity coupled element number of "// &
        & TRIM(NumberToVString(elementConnectivity%coupledElementNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    NULLIFY(coupledMeshTopology)
    CALL Mesh_MeshTopologyGet(coupledMesh,1,coupledMeshTopology,err,error,*999)
    NULLIFY(coupledMeshElements)
    CALL MeshTopology_MeshElementsGet(coupledMeshTopology,coupledMeshElements,err,error,*999)
    coupledMeshNumberOfXi = coupledMeshElements%elements(coupledElementNumber)%basis%numberOfXi
    IF(SIZE(xi,1)<coupledMeshNumberOfXi) THEN
      localError="The size of the specified xi array of "//TRIM(NumberToVString(SIZE(xi,1),"*",err,error))// &
        & " is too small. The size needs to be >= "//TRIM(NumberToVString(coupledMeshNumberOfXi,"*",err,error))//","
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    !\todo the elementConnectivity%xi array needs to be restructured to efficiently utilize memory when coupling bodies with 2xi directions to bodies with 3xi directions using an interface.
    elementConnectivity%xi(1:coupledMeshNumberOfXi,interfaceMeshComponentNumber,interfaceMeshLocalNodeNumber)= &
      & xi(1:coupledMeshNumberOfXi)
 
    EXITS("InterfaceMeshConnectivity_ElementXiSet")
    RETURN
999 ERRORSEXITS("InterfaceMeshConnectivity_ElementXiSet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMeshConnectivity_ElementXiSet

  !
  !================================================================================================================================
  !
  
  !>Sets the connectivity between an element in a coupled mesh to an element in the interface mesh
  SUBROUTINE InterfaceMeshConnectivity_ElementNumberSet(interfaceMeshConnectivity,interfaceMeshElementNumber, &
      & coupledMeshIndex,coupledElementNumber,err,error,*)

    !Argument variables
    TYPE(InterfaceMeshConnectivityType), POINTER :: interfaceMeshConnectivity !<A pointer to the interface mesh connectivity for the interface mesh
    INTEGER(INTG), INTENT(IN) :: interfaceMeshElementNumber !<The interface mesh element number to which the specified coupled mesh element would be connected
    INTEGER(INTG), INTENT(IN) :: coupledMeshIndex !<The index of the coupled mesh at the interface to set the element connectivity for
    INTEGER(INTG), INTENT(IN) :: coupledElementNumber !<The coupled mesh element to be connected to the interface
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("InterfaceMeshConnectivity_ElementNumberSet",err,error,*999)

    CALL InterfaceMeshConnectivity_AssertNotFinished(interfaceMeshConnectivity,err,error,*999)
    IF((interfaceMeshElementNumber<=0).OR.(interfaceMeshElementNumber>interfaceMeshConnectivity%numberOfInterfaceElements)) THEN
      localError="The specified interface mesh element number of "// &
        & TRIM(NumberToVString(interfaceMeshElementNumber,"*",err,error))// &
        & " is invalid. The element number should be >=1 and <= "// &
        & TRIM(NumberToVString(interfaceMeshConnectivity%numberOfInterfaceElements,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF((coupledMeshIndex<0).OR.(coupledMeshIndex>interfaceMeshConnectivity%numberOfCoupledMeshes)) THEN
      localError="The specified coupled mesh index of "// &
        & TRIM(NumberToVString(coupledMeshIndex,"*",err,error))// &
        & " is invalid. The coupled mesh index should be >=1 and <= "// &
        & TRIM(NumberToVString(interfaceMeshConnectivity%numberOfCoupledMeshes,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(interfaceMeshConnectivity%elementConnectivity)) &
      & CALL FlagError("The interface mesh connectivity element connectivity is not associated.",err,error,*999)

    interfaceMeshConnectivity%elementConnectivity(interfaceMeshElementNumber,coupledMeshIndex)%coupledElementNumber= &
      & coupledElementNumber
    
    EXITS("InterfaceMeshConnectivity_ElementNumberSet")
    RETURN
999 ERRORSEXITS("InterfaceMeshConnectivity_ElementNumberSet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMeshConnectivity_ElementNumberSet

  !
  !================================================================================================================================
  !
  
  !>Sets the connectivity between an element in a coupled mesh to an element in the interface mesh
  SUBROUTINE InterfaceMeshConnectivity_NodeNumbersSet(interfaceMeshConnectivity,interfaceMeshNodeNumbers,&
      & firstCoupledMeshIndex,firstCoupledMeshNodeNumbers,secondCoupledMeshIndex,secondCoupledMeshNodeNumbers, &
      & err,error,*)        
        
    !Argument variables
    TYPE(InterfaceMeshConnectivityType), POINTER :: interfaceMeshConnectivity !<A pointer to the interface mesh connectivity for the interface mesh
    INTEGER(INTG), INTENT(IN) :: interfaceMeshNodeNumbers(:) !<The interface mesh element number to which the specified coupled mesh element would be connected
    INTEGER(INTG), INTENT(IN) :: firstCoupledMeshIndex,secondCoupledMeshIndex !<The index of the coupled mesh at the interface to set the element connectivity for
    INTEGER(INTG), INTENT(IN) :: firstCoupledMeshNodeNumbers(:),secondCoupledMeshNodeNumbers(:) !<The coupled mesh element to be connected to the interface
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: nodeIndex
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("InterfaceMeshConnectivity_NodeNumbersSet",err,error,*999)

    CALL InterfaceMeshConnectivity_AssertNotFinished(interfaceMeshConnectivity,err,error,*999)
    IF(ALLOCATED(interfaceMeshConnectivity%coupledNodes)) &
      & CALL FlagError("Interface mesh connectivity coupled nodes is already allocated.",err,error,*999)
    IF(SIZE(interfaceMeshNodeNumbers,1)/=SIZE(firstCoupledMeshNodeNumbers,1)) THEN
      localError="The size of the specfied interface node numbers array of "// &
        & TRIM(NumberToVString(SIZE(interfaceMeshNodeNumbers,1),"*",err,error))// &
        & " does not match the size of the first coupled mesh node numbers array of "// &
        & TRIM(NumberToVString(SIZE(firstCoupledMeshNodeNumbers,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(interfaceMeshNodeNumbers,1)/=SIZE(secondCoupledMeshNodeNumbers,1)) THEN
      localError="The size of the specfied interface node numbers array of "// &
        & TRIM(NumberToVString(SIZE(interfaceMeshNodeNumbers,1),"*",err,error))// &
        & " does not match the size of the second coupled mesh node numbers array of "// &
        & TRIM(NumberToVString(SIZE(secondCoupledMeshNodeNumbers,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    !Default to two coupled meshes
    ALLOCATE(interfaceMeshConnectivity%coupledNodes(2,SIZE(interfaceMeshNodeNumbers(:))))
    DO nodeIndex=1,SIZE(interfaceMeshNodeNumbers(:))
      interfaceMeshConnectivity%coupledNodes(firstCoupledMeshIndex,nodeIndex)=firstCoupledMeshNodeNumbers(nodeIndex)
      interfaceMeshConnectivity%coupledNodes(secondCoupledMeshIndex,nodeIndex)=secondCoupledMeshNodeNumbers(nodeIndex)
    ENDDO !nodeIndex
    
    EXITS("InterfaceMeshConnectivity_NodeNumbersSet")
    RETURN
999 ERRORSEXITS("InterfaceMeshConnectivity_NodeNumbersSet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMeshConnectivity_NodeNumbersSet

  !
  !================================================================================================================================
  !

  !>Finalise the points connectivity coupled mesh elements 
  SUBROUTINE InterfaceCoupledElements_Finalise(coupledElements,err,error,*) 

    !Argument variables
    TYPE(InterfaceCoupledElementsType) :: coupledElements !<Coupled elements to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
  
    ENTERS("InterfaceCoupledElements_Finalise",err,error,*999)

    coupledElements%numberOfCoupledElements=0
    IF(ALLOCATED(coupledElements%elementNumbers)) DEALLOCATE(coupledElements%elementNumbers)
    
    EXITS("InterfaceCoupledElements_Finalise")
    RETURN
999 ERRORS("InterfaceCoupledElements_Finalise",err,error)
    EXITS("InterfaceCoupledElements_Finalise")
    RETURN 1
    
  END SUBROUTINE InterfaceCoupledElements_Finalise
  
  !
  !================================================================================================================================
  !

  !>Initialise the coupled mesh elements for points connectivity
  SUBROUTINE InterfaceCoupledElements_Initialise(coupledElements,err,error,*) 

    !Argument variables
    TYPE(InterfaceCoupledElementsType) :: coupledElements !<The coupled elements to initalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
  
    ENTERS("InterfaceCoupledElements_Initialise",err,error,*999)

    coupledElements%numberOfCoupledElements=0
    
    EXITS("InterfaceCoupledElements_Initialise")
    RETURN
999 ERRORS("InterfaceCoupledElements_Initialise",err,error)
    EXITS("InterfaceCoupledElements_Initialise")
    RETURN 1
    
  END SUBROUTINE InterfaceCoupledElements_Initialise
  
  !
  !================================================================================================================================
  !

  !>Finalise the interface element connectivity and deallocate all memory
  SUBROUTINE InterfaceElementConnectivity_Finalise(elementConnectivity,err,error,*)

    !Argument variables
    TYPE(InterfaceElementConnectivityType) :: elementConnectivity
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: interfaceElementIdx,coupledMeshIdx
     
    ENTERS("InterfaceElementConnectivity_Finalise",err,error,*999)

    elementConnectivity%coupledElementNumber=0
    elementConnectivity%connectedLineFace=0
    IF(ALLOCATED(elementConnectivity%xi)) DEALLOCATE(elementConnectivity%xi)
    
    EXITS("InterfaceElementConnectivity_Finalise")
    RETURN
999 ERRORSEXITS("InterfaceElementConnectivity_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceElementConnectivity_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the interface element connectivity.
  SUBROUTINE InterfaceElementConnectivity_Initialise(elementConnectivity,err,error,*)

    !Argument variables
    TYPE(InterfaceElementConnectivityType) :: elementConnectivity
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
     
    ENTERS("InterfaceElementConnectivity_Initialise",err,error,*999)

    elementConnectivity%coupledElementNumber=0
    elementConnectivity%connectedLineFace=0
    ! Note that elementConnectivity%xi(numberOfCoupledMeshXiDirections,1,numberOfInterfaceElementNodes) is allocated after the basis for the mesh connectivity has been specified in Interface_MeshConnectivityBasisSet where the number of numberOfInterfaceElementNodes can be determined.
    !\todo see corresponding todo in regarding updating the structure of elementConnectivity%XI
    
    EXITS("InterfaceElementConnectivity_Initialise")
    RETURN
999 ERRORSEXITS("InterfaceElementConnectivity_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceElementConnectivity_Initialise

  !
  !================================================================================================================================
  !

  !>Finalise interface point connectivity
  SUBROUTINE InterfacePointConnectivity_Finalise(interfacePointConnectivity,err,error,*) 

    !Argument variables
    TYPE(InterfacePointConnectivityType) :: interfacePointConnectivity !<A pointer to interface point connectivity to be finalised
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("InterfacePointConnectivity_Finalise",err,error,*999)

    interfacePointConnectivity%coupledElementNumber=0
    interfacePointConnectivity%elementLineFaceNumber=0
    IF(ALLOCATED(interfacePointConnectivity%xi)) DEALLOCATE(interfacePointConnectivity%xi)
    IF(ALLOCATED(interfacePointConnectivity%reducedXi)) DEALLOCATE(interfacePointConnectivity%reducedXi)
    
    EXITS("InterfacePointConnectivity_Finalise")
    RETURN
999 ERRORSEXITS("InterfacePointConnectivity_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE InterfacePointConnectivity_Finalise
  
  !
  !================================================================================================================================
  !

  !>Initialises the interface mesh connectivity.
  SUBROUTINE InterfacePointConnectivity_Initialise(interfacePointConnectivity,coupledMeshDimension,interfaceMeshDimension, &
      & err,error,*)

    !Argument variables
    TYPE(InterfacePointConnectivityType) :: interfacePointConnectivity !<An interface point connectivity to be initliased for
    INTEGER(INTG), INTENT(IN) :: coupledMeshDimension,interfaceMeshDimension
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr !<The error code
    TYPE(VARYING_STRING)  :: dummyError !<The error string
     
    ENTERS("InterfacePointsConnectivity_Initialise",err,error,*999)

    interfacePointConnectivity%coupledElementNumber=0
    interfacePointConnectivity%elementLineFaceNumber=0
    !Allocate memory for coupled mesh full and reduced xi location
    ALLOCATE(interfacePointConnectivity%xi(coupledMeshDimension),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate interface point connectivity full xi.",err,error,*999)
    interfacePointConnectivity%xi=0.0_DP
    ALLOCATE(interfacePointConnectivity%reducedXi(interfaceMeshDimension),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate interface point connectivity reduced xi.",err,error,*999)
    interfacePointConnectivity%reducedXi=0.0_DP
       
    EXITS("InterfacePointConnectivity_Initialise")
    RETURN
999 CALL InterfacePointConnectivity_Finalise(interfacePointConnectivity,dummyErr,dummyError,*998)
998 ERRORSEXITS("InterfacePointConnectivity_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE InterfacePointConnectivity_Initialise
  
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
    INTEGER(INTG) :: elementIdx,dataPointIdx,globalDataPointNumber,globalElementNumber,numberOfElementDataPoints, &
      & numberOfCoupledElements,coupledElementIdx
    INTEGER(INTG), ALLOCATABLE :: elementNumbers(:)
    TYPE(ListType), POINTER :: elementNumbersList
    TYPE(MeshTopologyType), POINTER :: meshTopology
    TYPE(MeshDataPointsType), POINTER :: meshDataPoints
  
    ENTERS("InterfacePointsConnectivity_CoupledElementsCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(interfacePointsConnectivity)) &
      & CALL FlagError("Interface points connectivity is not associated.",err,error,*999)
    
    IF(.NOT.ALLOCATED(interfacePointsConnectivity%coupledElements)) &
      & CALL FlagError("Interface points connectivity coupled elements is not allocated.",err,error,*999)
    
    interfacePointsConnectivity%maxNumberOfCoupledElements(coupledMeshIdx)=0; !Initialise the number of coupled mesh elements
    NULLIFY(meshTopology)
    CALL Mesh_MeshTopologyGet(interfacePointsConnectivity%interfaceMesh,1,meshTopology,err,error,*999)
    NULLIFY(meshDataPoints)
    CALL MeshTopology_MeshDataPointsGet(meshTopology,meshDataPoints,err,error,*999)
    DO elementIdx=1,SIZE(interfacePointsConnectivity%coupledElements,1) !Number of interface elements
      numberOfElementDataPoints=meshDataPoints%elementDataPoints(elementIdx)%numberOfProjectedData !Get the number of data points in interface mesh element
      !Set up list
      NULLIFY(elementNumbersList)
      CALL List_CreateStart(elementNumbersList,err,error,*999)
      CALL List_DataTypeSet(elementNumbersList,LIST_INTG_TYPE,err,error,*999)
      CALL List_InitialSizeSet(elementNumbersList,numberOfElementDataPoints,err,error,*999)
      CALL List_CreateFinish(elementNumbersList,err,error,*999)
      DO dataPointIdx=1,numberOfElementDataPoints
        globalDataPointNumber=meshDataPoints%elementDataPoints(elementIdx)%dataIndices(dataPointIdx)%globalNumber
        globalElementNumber=interfacePointsConnectivity%pointsConnectivity(globalDataPointNumber,coupledMeshIdx)% &
          & coupledElementNumber
        CALL List_ItemAdd(elementNumbersList,globalElementNumber,err,error,*999)
      ENDDO !dataPointIdx
      CALL List_RemoveDuplicates(elementNumbersList,err,error,*999)
      CALL List_DetachAndDestroy(elementNumbersList,numberOfCoupledElements,elementNumbers,err,error,*999)
      IF(ALLOCATED(interfacePointsConnectivity%coupledElements(elementIdx,coupledMeshIdx)%elementNumbers)) & 
        & DEALLOCATE(interfacePointsConnectivity%coupledElements(elementIdx,coupledMeshIdx)%elementNumbers) !for updating coupledElements after projection
      ALLOCATE(interfacePointsConnectivity%coupledElements(elementIdx,coupledMeshIdx)% &
        & elementNumbers(numberOfCoupledElements),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate coupled mesh element numbers.",err,error,*999)
      DO coupledElementIdx=1,numberOfCoupledElements
        interfacePointsConnectivity%coupledElements(elementIdx,coupledMeshIdx)%elementNumbers(coupledElementIdx)= &
          & elementNumbers(coupledElementIdx)
      ENDDO !coupledElementIdx
      interfacePointsConnectivity%coupledElements(elementIdx,coupledMeshIdx)%numberOfCoupledElements=numberOfCoupledElements
      IF(ALLOCATED(elementNumbers)) DEALLOCATE(elementNumbers)
      IF(interfacePointsConnectivity%maxNumberOfCoupledElements(coupledMeshIdx)<numberOfCoupledElements) THEN
        interfacePointsConnectivity%maxNumberOfCoupledElements(coupledMeshIdx)=numberOfCoupledElements
      ENDIF
    ENDDO !elementIdx
    
    EXITS("InterfacePointsConnectivity_CoupledElementsCalculate")
    RETURN
999 ERRORS("InterfacePointsConnectivity_CoupledElementsCalculate",err,error)
    EXITS("InterfacePointsConnectivity_CoupledElementsCalculate")
    RETURN 1
    
  END SUBROUTINE InterfacePointsConnectivity_CoupledElementsCalculate
  
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
    INTEGER(INTG) :: coupledMeshIdx
    TYPE(InterfaceType), POINTER :: interface
  
    ENTERS("InterfacePointsConnectivity_CreateFinish",err,error,*999)

    CALL InterfacePointsConnectivity_AssertNotFinished(interfacePointsConnectivity,err,error,*999)

    NULLIFY(INTERFACE)
    CALL InterfacePointsConnectivity_InterfaceGet(interfacePointsConnectivity,INTERFACE,err,error,*999)
    CALL InterfacePointsConnectivity_ReducedXiCalculate(interfacePointsConnectivity,err,error,*999) 
    DO coupledMeshIdx=1,interface%numberOfCoupledMeshes
      CALL InterfacePointsConnectivity_CoupledElementsCalculate(interfacePointsConnectivity,coupledMeshIdx,err,error,*999) 
    ENDDO !coupledMeshIdx
    interfacePointsConnectivity%pointsConnectivityFinished=.TRUE.   
    
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
    TYPE(InterfaceType), POINTER :: interface !<A pointer to the interface to create the points connectivity for
    TYPE(MeshType), POINTER :: interfaceMesh !<A pointer to the interface mesh for which the points connectivity is created
    TYPE(DataPointsType), POINTER :: dataPoints !<A pointer to the data points in the interface 
    TYPE(InterfacePointsConnectivityType), POINTER :: interfacePointsConnectivity !<On return, a pointer to the created points connectivity
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("InterfacePointsConnectivity_CreateStart",err,error,*999)

    CALL Interface_AssertIsFinished(INTERFACE,err,error,*999)
    CALL Mesh_AssertIsFinished(interfaceMesh,err,error,*999)
    CALL DataPoints_AssertIsFinished(dataPoints,err,error,*999)
    IF(ASSOCIATED(INTERFACE%pointsConnectivity)) &
      & CALL FlagError("The interface already has a points connectivity associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(INTERFACE,dataPoints%INTERFACE)) &
      & CALL FlagError("The specified data points are not defined on the specified interface.",err,error,*999)
    
    CALL Interface_PointsConnectivityInitialise(INTERFACE,interfaceMesh,dataPoints,err,error,*999)
    !Return the pointer
    interfacePointsConnectivity=>interface%pointsConnectivity
    
    EXITS("InterfacePointsConnectivity_CreateStart")
    RETURN
999 ERRORSEXITS("InterfacePointsConnectivity_CreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE InterfacePointsConnectivity_CreateStart
  
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

    IF(.NOT.ASSOCIATED(interfacePointsConnectivity)) &
      & CALL FlagError("Interface points connectivity is not associated.",err,error,*999)
    
    CALL Interface_PointsConnectivityFinalise(interfacePointsConnectivity,err,error,*999) 
     
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
    INTEGER(INTG) :: dataPointGlobalNumber
    LOGICAL :: dataPointExists

    ENTERS("InterfacePointsConnectivity_ElementNumberGet",err,error,*999)
    
    IF(.NOT.ASSOCIATED(pointsConnectivity))  CALL FlagError("Interface points connectivity is not associated.",err,error,*999)
    
    CALL DataPoint_CheckExists(pointsConnectivity%dataPoints,dataPointUserNumber,dataPointExists, &
        & dataPointGlobalNumber,err,error,*999)
    IF(dataPointExists) THEN
      coupledMeshUserElementNumber=pointsConnectivity%pointsConnectivity(dataPointGlobalNumber,coupledMeshIndexNumber)% &
        & coupledElementNumber
    ELSE
      CALL FlagError("Data point with user number of "//TRIM(NumberToVString(dataPointUserNumber,"*",err,error))// &
        & " does not exist.",err,error,*999)
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
    TYPE(DataPointsType), POINTER :: dataPoints
    TYPE(InterfaceType), POINTER :: interface
    TYPE(MeshType), POINTER :: coupledMesh
    TYPE(MeshElementsType), POINTER :: meshElements
    TYPE(MeshTopologyType), POINTER :: meshTopology
    TYPE(VARYING_STRING) :: localError

    ENTERS("InterfacePointsConnectivity_ElementNumberSet",err,error,*999)

    CALL InterfacePointsConnectivity_AssertNotFinished(pointsConnectivity,err,error,*999)
    IF(.NOT.ALLOCATED(pointsConnectivity%pointsConnectivity)) &
      & CALL FlagError("Interface points connectivity points connectivity is not allocated.",err,error,*999)

    NULLIFY(dataPoints)
    CALL InterfacePointsConnectivity_DataPointsGet(pointsConnectivity,dataPoints,err,error,*999)
    NULLIFY(INTERFACE)
    CALL InterfacePointsConnectivity_InterfaceGet(pointsConnectivity,INTERFACE,err,error,*999)
    
    CALL DataPoint_CheckExists(dataPoints,dataPointUserNumber,dataPointExists,dataPointGlobalNumber,err,error,*999)
    IF(.NOT.dataPointExists) THEN
      localError="A data point with a user number of "//TRIM(NumberToVString(dataPointUserNumber,"*",err,error))// &
        & " does not exist."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF((coupledMeshIndexNumber<=0).OR.(coupledMeshIndexNumber<=dataPoints%numberOfDataPoints)) THEN
      localError="The specified coupled mesh index number of "// &
        & TRIM(NumberToVString(coupledMeshIndexNumber,"*",err,error))// &
        & " is invalid. The coupled mesh index must be >= 1 and <= "// &
        & TRIM(NumberToVString(dataPoints%numberOfDataPoints,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    NULLIFY(coupledMesh)
    CALL Interface_CoupledMeshGet(interface,coupledMeshIndexNumber,coupledMesh,err,error,*999)
    NULLIFY(meshTopology)
    CALL Mesh_MeshTopologyGet(coupledMesh,meshComponentNumber,meshTopology,err,error,*999)
    NULLIFY(meshElements)
    CALL MeshTopology_MeshElementsGet(meshTopology,meshElements,err,error,*999)
    CALL MeshElements_ElementCheckExists(meshElements,coupledMeshUserElementNumber,elementExists,elementMeshNumber, &
      & err,error,*999) !Make sure user element exists       
    IF(.NOT.elementExists) THEN
      localError="An element with a user number of "//TRIM(NumberToVString(coupledMeshUserElementNumber,"*",err,error))// &
        & " does not exist."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    pointsConnectivity%pointsConnectivity(dataPointGlobalNumber,coupledMeshIndexNumber)%coupledElementNumber=elementMeshNumber
    
    EXITS("InterfacePointsConnectivity_ElementNumberSet")
    RETURN
999 ERRORS("InterfacePointsConnectivity_ElementNumberSet",err,error)
    EXITS("InterfacePointsConnectivity_ElementNumberSet")
    RETURN 1
    
  END SUBROUTINE InterfacePointsConnectivity_ElementNumberSet
  
  !
  !================================================================================================================================
  !

  !>Calculate full xi locations in points connectivity from reduced xi
  SUBROUTINE InterfacePointsConnectivity_FullXiCalculate(interfacePointsConnectivity,coupledMeshIdx,err,error,*) 

    !Argument variables
    TYPE(InterfacePointsConnectivityType), POINTER :: interfacePointsConnectivity !<A pointer to the interface points connectivity to finish creating
    INTEGER(INTG), INTENT(IN) :: coupledMeshIdx !<Coupled mesh index to calculate the full xi location for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointIdx,interfaceMeshDimensions,coupledMeshDimensions
    TYPE(InterfaceType), POINTER :: interface
    TYPE(InterfacePointConnectivityType), POINTER :: pointConnectivity
    TYPE(MeshType), POINTER :: interfaceMesh,coupledMesh
    TYPE(VARYING_STRING) :: localError
  
    ENTERS("InterfacePointsConnectivity_FullXiCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(interfacePointsConnectivity)) &
      & CALL FlagError("Interface points connectivity is not associated.",err,error,*999)
    NULLIFY(interfaceMesh)
    CALL InterfacePointsConnectivity_InterfaceMeshGet(interfacePointsConnectivity,interfaceMesh,err,error,*999)
    NULLIFY(INTERFACE)
    CALL InterfacePointsConnectivity_InterfaceGet(interfacePointsConnectivity,INTERFACE,err,error,*999)
    NULLIFY(coupledMesh)
    CALL Interface_CoupledMeshGet(INTERFACE,coupledMeshIdx,coupledMesh,err,error,*999)    
    interfaceMeshDimensions=interfaceMesh%numberOfDimensions
    coupledMeshDimensions=coupledMesh%numberOfDimensions
    IF(interfaceMeshDimensions==coupledMeshDimensions) THEN !e.g. If 1D-2D, 2D-3D coupling, interface dimension is 1D and 2D respectively for 1st body
      DO dataPointIdx=1,InterfacePointsConnectivity%dataPoints%numberOfDataPoints
        interfacePointsConnectivity%pointsConnectivity(dataPointIdx,coupledMeshIdx)%xi(:)= &
          & interfacePointsConnectivity%pointsConnectivity(dataPointIdx,coupledMeshIdx)%reducedXi(:)
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
              & TRIM(NumberToVString(pointConnectivity%elementLineFaceNumber,"*",err,error))//" ."
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
              & TRIM(NumberToVString(pointConnectivity%elementLineFaceNumber,"*",err,error))//" ."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ENDDO !dataPointIdx
      CASE DEFAULT
        localError="Invalid coupled mesh dimension "// &
          & TRIM(NumberToVString(coupledMeshDimensions,"*",err,error))//" ."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDIF
     
    EXITS("InterfacePointsConnectivity_FullXiCalculate")
    RETURN
999 ERRORSEXITS("InterfacePointsConnectivity_FullXiCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE InterfacePointsConnectivity_FullXiCalculate
  
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
    REAL(DP), INTENT(OUT) :: xi(:) !<xi(xiIdx). The full xi location in the coupled mesh element.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointGlobalNumber
    LOGICAL :: dataPointExists
    TYPE(InterfaceType), POINTER :: INTERFACE
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("InterfacePointsConnectivity_PointXiGet",err,error,*999)
    
    ! Preliminary error checks to verify user input information
    IF(.NOT.ASSOCIATED(pointsConnectivity)) CALL FlagError("Points connectivity is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(pointsConnectivity%pointsConnectivity)) &
      & CALL FlagError("Interface points connectivity points connectivity is not allocated.",err,error,*999)
    NULLIFY(INTERFACE)
    CALL InterfacePointsConnectivity_InterfaceGet(pointsConnectivity,INTERFACE,err,error,*999)
    IF(coupledMeshIndexNumber<=0.OR.coupledMeshIndexNumber>interface%numberOfCoupledMeshes) THEN
      localError="The specified coupled mesh index number of "// &
        & TRIM(NumberToVString(coupledMeshIndexNumber,"*",err,error))// &
        & " is invalid. The couple mesh index should be >= 1 and <= "// &
        & TRIM(NumberToVString(INTERFACE%numberOfCoupledMeshes,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF    
    CALL DataPoint_CheckExists(pointsConnectivity%dataPoints,dataPointUserNumber,dataPointExists, &
      & dataPointGlobalNumber,err,error,*999)
    IF(.NOT.dataPointExists) THEN
      localError="A data point with an user number of "// &
        & TRIM(NumberToVString(dataPointUserNumber,"*",err,error))// &
        & " does not exist."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(xi,1)<SIZE(pointsConnectivity%pointsConnectivity(dataPointUserNumber,coupledMeshIndexNumber)%xi,1)) THEN
      localError="The size of the specified xi array of "// &
        & TRIM(NumberToVString(SIZE(xi,1),"*",err,error))// &
        & " is too small. The size should be >= "// &
        & TRIM(NumberToVString(SIZE(pointsConnectivity%pointsConnectivity(dataPointUserNumber, &
        & coupledMeshIndexNumber)%xi,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    xi=pointsConnectivity%pointsConnectivity(dataPointUserNumber,coupledMeshIndexNumber)%xi(:)

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
    REAL(DP), INTENT(IN) :: xi(:) !<xi(xiIdx). The full xi location in the coupled mesh element.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointGlobalNumber
    LOGICAL :: dataPointExists
    TYPE(InterfaceType), POINTER :: INTERFACE
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("InterfacePointsConnectivity_PointXiSet",err,error,*999)
    
    ! Preliminary error checks to verify user input information
    CALL InterfacePointsConnectivity_AssertNotFinished(pointsConnectivity,err,error,*999)
    IF(.NOT.ALLOCATED(pointsConnectivity%pointsConnectivity)) &
      & CALL FlagError("Interface points connectivity points connectivity is not allocated.",err,error,*999)
    NULLIFY(INTERFACE)
    CALL InterfacePointsConnectivity_InterfaceGet(pointsConnectivity,INTERFACE,err,error,*999)
    IF(coupledMeshIndexNumber<=0.OR.coupledMeshIndexNumber>interface%numberOfCoupledMeshes) THEN
      localError="The specified coupled mesh index number of "// &
        & TRIM(NumberToVString(coupledMeshIndexNumber,"*",err,error))// &
        & " is invalid. The couple mesh index should be >= 1 and <= "// &
        & TRIM(NumberToVString(INTERFACE%numberOfCoupledMeshes,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF    
    CALL DataPoint_CheckExists(pointsConnectivity%dataPoints,dataPointUserNumber,dataPointExists, &
      & dataPointGlobalNumber,err,error,*999)
    IF(.NOT.dataPointExists) THEN
      localError="A data point with an user number of "// &
        & TRIM(NumberToVString(dataPointUserNumber,"*",err,error))// &
        & " does not exist."
      CALL FlagError(localError,err,error,*999)
    ENDIF    
    IF(SIZE(xi,1)<SIZE(pointsConnectivity%pointsConnectivity(dataPointGlobalNumber,coupledMeshIndexNumber)%xi,1)) THEN
      localError="The size of the specified xi array of "// &
        & TRIM(NumberToVString(SIZE(xi,1),"*",err,error))// &
        & " is too small. The size should be >= "// &
        & TRIM(NumberToVString(SIZE(pointsConnectivity%pointsConnectivity(dataPointUserNumber, &
        & coupledMeshIndexNumber)%xi,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    pointsConnectivity%pointsConnectivity(dataPointGlobalNumber,coupledMeshIndexNumber)%xi(:)=xi(:) 
 
    EXITS("InterfacePointsConnectivity_PointXiSet")
    RETURN
999 ERRORSEXITS("InterfacePointsConnectivity_PointXiSet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfacePointsConnectivity_PointXiSet
  
  !
  !================================================================================================================================
  !
  
  !>Update points connectivity with projection results
  SUBROUTINE InterfacePointsConnectivity_UpdateFromProjection(interfacePointsConnectivity,dataProjection, &
      & coupledMeshIndex,err,error,*) 
  
    !Argument variables
    TYPE(InterfacePointsConnectivityType), POINTER :: interfacePointsConnectivity !<A pointer to the interface points connectivity to finish creating
    TYPE(DataProjectionType), POINTER :: dataProjection !<The data projection that points connectivity update with
    INTEGER(INTG), INTENT(IN) :: coupledMeshIndex !<The mesh index of the the points connectivity to be updated
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointIdx
    TYPE(DataProjectionResultType), POINTER :: dataProjectionResult
    
    ENTERS("InterfacePointsConnectivity_UpdateFromProjection",err,error,*999)
    
    IF(.NOT.ASSOCIATED(interfacePointsConnectivity)) &
      &  CALL FlagError("Interface points connectivity is not associated.",err,error,*999)
    CALL DataProjection_AssertIsFinished(dataProjection,err,error,*999)
    
    !Update reduced xi location, projection element number and element face/line number with projection results
    DO dataPointIdx=1,SIZE(dataProjection%dataProjectionResults,1) 
      dataProjectionResult=>dataProjection%dataProjectionResults(dataPointIdx)
      InterfacePointsConnectivity%pointsConnectivity(dataPointIdx,coupledMeshIndex)%reducedXi(:)= &
        & dataProjectionResult%xi
      InterfacePointsConnectivity%pointsConnectivity(dataPointIdx,coupledMeshIndex)%coupledElementNumber= &
        & dataProjectionResult%elementNumber
      InterfacePointsConnectivity%pointsConnectivity(dataPointIdx,coupledMeshIndex)%elementLineFaceNumber= &
        & dataProjectionResult%elementLineFaceNumber
    ENDDO !dataPointIdx
    CALL InterfacePointsConnectivity_FullXiCalculate(InterfacePointsConnectivity,coupledMeshIndex,err,error,*999) 
    !Update points connectivity coupledElement information
    CALL InterfacePointsConnectivity_CoupledElementsCalculate(InterfacePointsConnectivity,coupledMeshIndex,err,error,*999) 
     
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
  SUBROUTINE InterfacePointsConnectivity_ReducedXiCalculate(interfacePointsConnectivity,err,error,*) 

    !Argument variables
    TYPE(InterfacePointsConnectivityType), POINTER :: interfacePointsConnectivity !<A pointer to the interface points connectivity to finish creating
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: meshIdx,dataPointIdx,xiIdx,interfaceMeshDimensions,coupledMeshDimensions
    TYPE(DataPointsType), POINTER :: dataPoints
    TYPE(InterfaceType), POINTER :: INTERFACE
    TYPE(MeshType), POINTER :: coupledMesh,interfaceMesh
  
    ENTERS("InterfacePointsConnectivity_ReducedXiCalculate",err,error,*999)

    CALL InterfacePointsConnectivity_AssertNotFinished(interfacePointsConnectivity,err,error,*999)
    NULLIFY(INTERFACE)
    CALL InterfacePointsConnectivity_InterfaceGet(interfacePointsConnectivity,INTERFACE,err,error,*999)
    NULLIFY(interfaceMesh)
    CALL InterfacePointsConnectivity_InterfaceMeshGet(interfacePointsConnectivity,interfaceMesh,err,error,*999)
    interfaceMeshDimensions=interfaceMesh%numberOfDimensions
    NULLIFY(dataPoints)
    CALL InterfacePointsConnectivity_DataPointsGet(interfacePointsConnectivity,dataPoints,err,error,*999)
    DO meshIdx=1,INTERFACE%numberOfCoupledMeshes
      NULLIFY(coupledMesh)
      CALL Interface_CoupledMeshGet(INTERFACE,meshIdx,coupledMesh,err,error,*999)
      coupledMeshDimensions=coupledMesh%numberOfDimensions
      IF(interfaceMeshDimensions==coupledMeshDimensions) THEN
        !e.g. If 1D-2D, 2D-3D coupling, interface dimension is 1D and 2D respectively for 1st body
        DO dataPointIdx=1,dataPoints%numberOfDataPoints
          interfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%reducedXi(:)= &
            & interfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%xi(:)
          interfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%elementLineFaceNumber=1 !The only local face/line for the body with lower dimension 
        ENDDO !dataPointIdx
      ELSE
        SELECT CASE(coupledMeshDimensions)
        CASE(2)
          DO dataPointIdx=1,dataPoints%numberOfDataPoints
            DO xiIdx=1,coupledMeshDimensions
              IF(ABS(interfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%xi(xiIdx))< ZERO_TOLERANCE) THEN
                !Calculate line number
                interfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%elementLineFaceNumber=4-(xiIdx-1)*2 
                EXIT
              ELSE IF(ABS(interfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%xi(xiIdx)-1.0_DP) < &
                & ZERO_TOLERANCE) THEN
                !Calculate line number
                interfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%elementLineFaceNumber=3-(xiIdx-1)*2 
                EXIT
              ENDIF
            ENDDO !xiIdx
            interfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%reducedXi= &
              & interfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)% &
              & xi(OTHER_XI_DIRECTIONS2(xiIdx))  !Populate reducedXi 
          ENDDO !dataPointIdx
        CASE(3)
          DO dataPointIdx=1,interfacePointsConnectivity%dataPoints%numberOfDataPoints
            DO xiIdx=1,coupledMeshDimensions
              IF(ABS(interfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%xi(xiIdx)) &
                & < ZERO_TOLERANCE) THEN
                !Calculate face number
                interfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%elementLineFaceNumber=(xiIdx-1)*2+2 
                EXIT
              ELSE IF(ABS(interfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%xi(xiIdx)-1.0_DP) &
                & > ZERO_TOLERANCE) THEN
                !Calculate face number
                interfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%elementLineFaceNumber=(xiIdx-1)*2+1 
                EXIT
              ENDIF
            ENDDO !xiIdx
!!TODO: This seems wrong as the xiIdx loop has finished and so xiIdx is meaningless???
            SELECT CASE(xiIdx) !Populate reducedXi 
            CASE(1)
              interfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%reducedXi(1)= &
                & interfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%xi(2)
              interfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%reducedXi(2)= &
                & interfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%xi(3)
            CASE(2)
              interfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%reducedXi(1)= &
                & interfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%xi(1)
              interfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%reducedXi(2)= &
                & interfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%xi(3)
            CASE(3)
              interfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%reducedXi(1)= &
                & interfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%xi(1)
              interfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%reducedXi(2)= &
                & interfacePointsConnectivity%pointsConnectivity(dataPointIdx,meshIdx)%xi(2)
            END SELECT
          ENDDO !dataPointIdx
        CASE DEFAULT
          ! Do nothing
        END SELECT
      ENDIF
    ENDDO !meshIdx
    
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
  SUBROUTINE Interfaces_Finalise(interfaces,err,error,*) 

    !Argument variables
    TYPE(InterfacesType), POINTER :: interfaces !<A pointer to the interfaces to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(InterfaceType), POINTER :: INTERFACE
     
    ENTERS("Interfaces_Finalise",err,error,*999)

    IF(ASSOCIATED(interfaces)) THEN
      DO WHILE(interfaces%numberOfInterfaces>0)
        INTERFACE=>interfaces%interfaces(1)%ptr
        CALL Interface_Destroy(INTERFACE,err,error,*999)
      ENDDO
      IF(ALLOCATED(interfaces%interfaces)) DEALLOCATE(interfaces%interfaces)
      DEALLOCATE(interfaces)
    ENDIF
    
    EXITS("Interfaces_Finalise")
    RETURN
999 ERRORSEXITS("Interfaces_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE Interfaces_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises interfaces for a region.
  SUBROUTINE Interfaces_Initialise(region,err,error,*) 

    !Argument variables
    TYPE(RegionType), POINTER :: region !<A pointer to the region to initialise the interfaces for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError,localError
     
    ENTERS("Interfaces_Initialise",err,error,*998)

    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated.",err,error,*999)
    IF(ASSOCIATED(region%interfaces)) THEN
      localError="Region number "//TRIM(NumberToVString(region%userNumber,"*",err,error))// &
        & " already has interfaces associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    ALLOCATE(region%interfaces,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate region interfaces.",err,error,*999)
    region%interfaces%parentRegion=>region
    region%interfaces%numberOfInterfaces=0
    
    EXITS("Interfaces_Initialise")
    RETURN
999 CALL Interfaces_Finalise(region%interfaces,dummyErr,dummyError,*998)
998 ERRORSEXITS("Interfaces_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE Interfaces_Initialise

  !
  !================================================================================================================================
  !
  
END MODULE InterfaceRoutines
