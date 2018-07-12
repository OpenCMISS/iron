!> \file
!> \author Chris Bradley
!> \brief This module handles all decomposition routines.
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

!>This module handles all decomposition routines.
MODULE DecompositionRoutines

  USE BaseRoutines
  USE BasisRoutines
  USE BasisAccessRoutines
  USE CmissMPI
  USE CMISS_ParMETIS
  USE ComputationRoutines
  USE ComputationAccessRoutines
  USE Constants
  USE ContextAccessRoutines
  USE DecompositionAccessRoutines
  USE DOMAIN_MAPPINGS
  USE INPUT_OUTPUT
  USE InterfaceAccessRoutines
  USE ISO_VARYING_STRING
  USE Kinds
  USE Lists
  USE MeshAccessRoutines
#ifndef NOMPIMOD
  USE MPI
#endif
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
  
  PUBLIC Decomposer_CreateStart,Decomposer_CreateFinish

  PUBLIC Decomposer_DecompositionAdd

  PUBLIC Decomposer_Destroy

  PUBLIC Decomposers_Initialise,Decomposers_Finalise

  PUBLIC Decomposition_CalculateLinesSet,Decomposition_CalculateFacesSet

  PUBLIC Decomposition_CreateStart,Decomposition_CreateFinish
  
  PUBLIC Decomposition_Destroy

  PUBLIC Decomposition_ElementBasisGet

  PUBLIC Decomposition_ElementDomainCalculate

  PUBLIC Decomposition_ElementDomainGet,Decomposition_ElementDomainSet

  PUBLIC Decomposition_MeshComponentNumberGet,Decomposition_MeshComponentNumberSet
  
  PUBLIC Decomposition_NodeDomainGet

  PUBLIC DECOMPOSITION_NUMBER_OF_DOMAINS_GET,DECOMPOSITION_NUMBER_OF_DOMAINS_SET

  PUBLIC Decomposition_TypeGet,Decomposition_TypeSet
  
  PUBLIC Decomposition_WorkGroupSet

  PUBLIC Decomposition_DataProjectionCalculate

  PUBLIC Decomposition_ElementDataPointLocalNumberGet

  PUBLIC Decomposition_ElementDataPointUserNumberGet

  PUBLIC Decomposition_NumberOfElementDataPointsGet
  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Calculate the decompositions for a decomposer.
  SUBROUTINE Decomposer_Calculate(decomposer,err,error,*)

    !Argument variables
    TYPE(DecomposerType), POINTER :: decomposer !<A pointer to the decomposer to calculate the decompositions for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: graphRootIdx
    TYPE(DecomposerGraphType), POINTER :: decomposerGraph
    TYPE(DecomposerGraphNodeType), POINTER :: graphRoot

    ENTERS("Decomposer_Calculate",err,error,*999)

    CALL Decomposer_AssertIsFinished(decomposer,err,error,*999)
    NULLIFY(decomposerGraph)
    CALL Decomposer_DecomposerGraphGet(decomposer,decomposerGraph,err,error,*999)
    
    DO graphRootIdx=1,decomposerGraph%numberOfGraphRoots
      NULLIFY(graphRoot)
      CALL DecomposerGraph_RootNodeGet(decomposerGraph,graphRootIdx,graphRoot,err,error,*999)
      CALL Decomposer_ElementDomainCalculate(graphRoot,err,error,*999)
    ENDDO !graphRootIdx
       
    EXITS("Decomposer_Calculate")
    RETURN
999 ERRORSEXITS("Decomposer_Calculate",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposer_Calculate

  !
  !================================================================================================================================
  !

  !>Finish the creation of decomposer.
  SUBROUTINE Decomposer_CreateFinish(decomposer,err,error,*)

    !Argument variables
    TYPE(DecomposerType), POINTER :: decomposer !<A pointer to the decomposer to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Decomposer_CreateFinish",err,error,*999)

    CALL Decomposer_AssertNotFinished(decomposer,err,error,*999)
    IF(decomposer%numberOfDecompositions<1) CALL FlagError("No decompositions have been added to the decomposer.",err,error,*999)
    
    !Calculate the decomposer graph
    CALL Decomposer_GraphCalculate(decomposer,err,error,*999)
    !Set the finished flag
    decomposer%decomposerFinished=.TRUE.
    !Calculate the decompositions for the decomposer
    CALL Decomposer_Calculate(decomposer,err,error,*999)    
    !Calculate the decomposition topology
    CALL Decomposer_TopologyCalculate(decomposer,err,error,*999)
           
    EXITS("Decomposer_CreateFinish")
    RETURN
999 ERRORSEXITS("Decomposer_CreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposer_CreateFinish

  !
  !================================================================================================================================
  !

  !>Start the creation of decomposer for a region.
  SUBROUTINE Decomposer_CreateStart(decomposerUserNumber,region,workGroup,decomposer,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: decomposerUserNumber !<The user number of the decomposer to create
    TYPE(RegionType), POINTER :: region !<A pointer to the region to create decomposer for
    TYPE(WorkGroupType), POINTER :: workGroup !<A pointer to the work group for the created decomposer
    TYPE(DecomposerType), POINTER :: decomposer !<On exit, a pointer to the created decomposer. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: decomposerIdx,dummyErr
    TYPE(DecomposerType), POINTER :: newDecomposer
    TYPE(DecomposerPtrType), ALLOCATABLE :: newDecomposers(:)
    TYPE(DecomposersType), POINTER :: decomposers
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("Decomposer_CreateStart",err,error,*998)

    IF(ASSOCIATED(decomposer)) CALL FlagError("Decomposer is already associated.",err,error,*998)
    CALL Region_AssertIsFinished(region,err,error,*998)
    CALL WorkGroup_AssertIsFinished(workGroup,err,error,*998)
    NULLIFY(decomposers)
    CALL Region_DecomposersGet(region,decomposers,err,error,*998)

    !Check if the user number has already been used.
    NULLIFY(newDecomposer)
    CALL Decomposer_UserNumberFind(decomposerUserNumber,region,newDecomposer,err,error,*998)
    IF(ASSOCIATED(newDecomposer)) THEN
      localError="Decomposer user number "//TRIM(NumberToVString(decomposerUserNumber,"*",err,error))// &
        & " has already been created on region number "//TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*998)
    ENDIF

    !Initialise the decomposer
    CALL Decomposer_Initialise(newDecomposer,err,error,*999)
    !Set default decomposer values
    newDecomposer%userNumber=decomposerUserNumber
    newDecomposer%decomposers=>decomposers
    newDecomposer%region=>region
    newDecomposer%workGroup=>workGroup
    !Add the decomposer to the list of decomposers for the region
    ALLOCATE(newDecomposers(decomposers%numberOfDecomposers+1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new decomposers.",err,error,*999)
    DO decomposerIdx=1,decomposers%numberOfDecomposers
      newDecomposers(decomposerIdx)%ptr=>decomposers%decomposers(decomposerIdx)%ptr
    ENDDO !decomposerIdx
    newDecomposers(decomposers%numberOfDecomposers+1)%ptr=>newDecomposer
    CALL MOVE_ALLOC(newDecomposers,decomposers%decomposers)
    decomposers%numberOfDecomposers=decomposers%numberOfDecomposers+1
    !Return the pointer
    decomposer=>newDecomposer
       
    EXITS("Decomposer_CreateStart")
    RETURN
999 IF(ALLOCATED(newDecomposers)) DEALLOCATE(newDecomposers)
    CALL Decomposer_Finalise(newDecomposer,dummyErr,dummyError,*998)
998 ERRORSEXITS("Decomposer_CreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposer_CreateStart

  !
  !================================================================================================================================
  !

  !>Add a decomposition to a decomposer.
  SUBROUTINE Decomposer_DecompositionAdd(decomposer,decomposition,decompositionIndex,err,error,*)

    !Argument variables
    TYPE(DecomposerType), POINTER :: decomposer !<A pointer to the decomposer to add decompositions to
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer to the decomposition to add to the decomposer
    INTEGER(INTG), INTENT(OUT) :: decompositionIndex !<On return, the index of the decomposition that has been added to the decomposer
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: decompositionIdx
    LOGICAL :: isSubRegion
    TYPE(DecompositionPtrType), ALLOCATABLE :: newDecompositions(:)
    TYPE(MeshType), POINTER :: mesh
    TYPE(RegionType), POINTER :: decomposerRegion,decompositionRegion
    TYPE(VARYING_STRING) :: localError

    ENTERS("Decomposer_DecompositionAdd",err,error,*999)

    CALL Decomposer_AssertNotFinished(decomposer,err,error,*999)
    CALL Decomposition_AssertIsFinished(decomposition,err,error,*999)

    NULLIFY(decomposerRegion)
    CALL Decomposer_RegionGet(decomposer,decomposerRegion,err,error,*999)
    NULLIFY(decompositionRegion)
    CALL Decomposition_RegionGet(decomposition,decompositionRegion,err,error,*999)
    CALL Region_IsSubRegion(decomposerRegion,decompositionRegion,isSubRegion,err,error,*999)
    IF(.NOT.isSubRegion) THEN
      localError="Decomposition number "//TRIM(NumberToVString(decomposition%userNumber,"*",err,error))// &
        & " is in region number "//TRIM(NumberToVString(decompositionRegion%userNumber,"*",err,error))// &
        & " which is not a sub-region of region number "//TRIM(NumberToVString(decomposerRegion%userNumber,"*",err,error))// &
        & " which has decomposer number "//TRIM(NumberToVString(decomposer%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(ASSOCIATED(decomposition%decomposer)) THEN
      localError="Decomposition number "//TRIM(NumberToVString(decomposition%userNumber,"*",err,error))// &
        & " already has an decomposer associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    !Get the decomposer mesh
    NULLIFY(mesh)
    CALL Decomposition_MeshGet(decomposition,mesh,err,error,*999)
    
    !Check to see if the decomposition has already been added or if the mesh has already been added
    DO decompositionIdx=1,decomposer%numberOfDecompositions
      IF(ASSOCIATED(decomposition,decomposer%decompositions(decompositionIdx)%ptr)) THEN
        localError="Decomposition number "//TRIM(NumberToVString(decomposition%userNumber,"*",err,error))// &
          & " has already been added to decomposer number "//TRIM(NumberToVString(decomposer%userNumber,"*",err,error))// &
          & " at position index "//TRIM(NumberToVString(decompositionIdx,"*",err,error))// &
          & ". The decomposition can not be added more than once."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      IF(ASSOCIATED(mesh,decomposer%decompositions(decompositionIdx)%ptr%mesh)) THEN
        localError="Mesh number "//TRIM(NumberToVString(mesh%userNumber,"*",err,error))// &
          & " with decomposition number "//TRIM(NumberToVString(decomposition%userNumber,"*",err,error))// &
          & " has already been added to decomposer number "//TRIM(NumberToVString(decomposer%userNumber,"*",err,error))// &
          & " at position index "//TRIM(NumberToVString(decompositionIdx,"*",err,error))// &
          & ". The mesh can not be added more than once."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDDO !decompositionIdx

    !Add the decomposition to the list of decomposer decompositions
    ALLOCATE(newDecompositions(decomposer%numberOfDecompositions+1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new decomposition.",err,error,*999)
    DO decompositionIdx=1,decomposer%numberOfDecompositions
      newDecompositions(decompositionIdx)%ptr=>decomposer%decompositions(decompositionIdx)%ptr
    ENDDO !decompositionIdx
    newDecompositions(decomposer%numberOfDecompositions+1)%ptr=>decomposition
    CALL MOVE_ALLOC(newDecompositions,decomposer%decompositions)
    decomposer%numberOfDecompositions=decomposer%numberOfDecompositions+1
    decompositionIndex=decomposer%numberOfDecompositions
       
    EXITS("Decomposer_DecompositionAdd")
    RETURN
999 ERRORSEXITS("Decomposer_DecompositionAdd",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposer_DecompositionAdd

  !
  !================================================================================================================================
  !

  !>Destroys a decomposer
  SUBROUTINE Decomposer_Destroy(decomposer,err,error,*)

    !Argument variables
    TYPE(DecomposerType), POINTER :: decomposer !<A pointer to the decomposer to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Decomposer_Destroy",err,error,*999)

    IF(.NOT.ASSOCIATED(decomposer)) CALL FlagError("Decomposer is not associated.",err,error,*999)
    
    CALL Decomposer_Finalise(decomposer,err,error,*999)
       
    EXITS("Decomposer_Destroy")
    RETURN
999 ERRORSEXITS("Decomposer_Destroy",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposer_Destroy

  !
  !================================================================================================================================
  !

  !>Calculate the element domain for a decomposer graph root node.
  SUBROUTINE Decomposer_ElementDomainCalculate(rootNode,err,error,*)

    !Argument variables
    TYPE(DecomposerGraphNodeType), POINTER :: rootNode !<A pointer to the decomposer graph root node to calculate the element domains for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: decompositionIdx,groupCommunicator,myComputationNodeNumber,numberOfComputationNodes,numberOfElements, &
      & numberOfMeshElements,numberOfElementsPerNode
    INTEGER(INTG), ALLOCATABLE :: elementOffset(:)
    LOGICAL :: isInterfaceMesh
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecomposerType), POINTER :: decomposer
    TYPE(DecomposerGraphType), POINTER :: decomposerGraph
    TYPE(MeshType), POINTER :: mesh

    ENTERS("Decomposer_ElementDomainCalculate",err,error,*999)

    NULLIFY(decomposerGraph)
    CALL DecomposerGraphNode_DecomposerGraphGet(rootNode,decomposerGraph,err,error,*999)
    NULLIFY(decomposer)
    CALL DecomposerGraph_DecomposerGet(decomposerGraph,decomposer,err,error,*999)
    CALL Decomposer_AssertIsFinished(decomposer,err,error,*999)
    
    CALL WorkGroup_GroupCommunicatorGet(decomposer%workGroup,groupCommunicator,err,error,*999)
    CALL WorkGroup_NumberOfGroupNodesGet(decomposer%workGroup,numberOfComputationNodes,err,error,*999)
    CALL WorkGroup_GroupNodeNumberGet(decomposer%workGroup,myComputationNodeNumber,err,error,*999)

    IF(numberOfComputationNodes==1) THEN
    ELSE
      !Calculate the number of elements in the combined decomposer mesh
      ALLOCATE(elementOffset(decomposer%numberOfDecompositions),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate element offset.",err,error,*999)
      numberOfElements=0
      DO decompositionIdx=1,decomposer%numberOfDecompositions
        NULLIFY(decomposition)
        CALL Decomposer_DecompositionGet(decomposer,decompositionIdx,decomposition,err,error,*999)
        NULLIFY(mesh)
        CALL Decomposition_MeshGet(decomposition,mesh,err,error,*999)
        CALL Mesh_IsInterfaceMesh(mesh,isInterfaceMesh,err,error,*999)
        elementOffset(decompositionIdx)=numberOfElements
        IF(.NOT.isInterfaceMesh) THEN
          CALL Mesh_NumberOfElementsGet(mesh,numberOfMeshElements,err,error,*999)
          numberOfElements=numberOfElements+numberOfMeshElements
        ENDIF
      ENDDO !decompositionIdx
      numberOfElementsPerNode=NINT(REAL(numberOfElements,DP)/REAL(numberOfComputationNodes,DP))
    ENDIF
      
    EXITS("Decomposer_ElementDomainCalculate")
    RETURN
999 IF(ALLOCATED(elementOffset)) DEALLOCATE(elementOffset)
    ERRORSEXITS("Decomposer_ElementDomainCalculate",err,error)    
    RETURN 1
    
  END SUBROUTINE Decomposer_ElementDomainCalculate

  !
  !================================================================================================================================
  !

  !>Calculate the decomposer graph.
  SUBROUTINE Decomposer_GraphCalculate(decomposer,err,error,*)

    !Argument variables
    TYPE(DecomposerType), POINTER :: decomposer !<A pointer to the decomposer to calculate the graph for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: coupledMeshIdx,decompositionIdx,graphLinkIdx,graphRootIdx,matchedNodeIdx,numberOfGraphRoots, &
      & numberOfMatchedGraphNodes
    LOGICAL :: coupledMeshFound,isInterfaceDecomposition,isRegionDecomposition
    TYPE(DecomposerGraphType), POINTER :: decomposerGraph
    TYPE(DecomposerGraphLinkType), POINTER :: graphLink
    TYPE(DecomposerGraphLinkPtrType), ALLOCATABLE :: newGraphLinks(:)
    TYPE(DecomposerGraphNodeType), POINTER :: graphNode,graphRoot,linkedNode,newGraphNode
    TYPE(DecomposerGraphNodePtrType), ALLOCATABLE :: matchedGraphNodes(:),newGraphRoots(:)
    TYPE(DecompositionType), POINTER :: decomposition,graphRootDecomposition,linkedNodeDecomposition
    TYPE(InterfaceType), POINTER :: INTERFACE
    TYPE(MeshType), POINTER :: coupledMesh,decompositionMesh
    TYPE(VARYING_STRING) :: localError

    ENTERS("Decomposer_GraphCalculate",err,error,*999)

    CALL Decomposer_AssertNotFinished(decomposer,err,error,*999)

    !Initialise the decomposition graph
    CALL Decomposer_GraphInitialise(decomposer,err,error,*999)
    NULLIFY(decomposerGraph)
    CALL Decomposer_DecomposerGraphGet(decomposer,decomposerGraph,err,error,*999)
    !Loop over the decompositions added to construct the graph
    !First add all the region decompositions
    DO decompositionIdx=1,decomposer%numberOfDecompositions
      NULLIFY(decomposition)
      CALL Decomposer_DecompositionGet(decomposer,decompositionIdx,decomposition,err,error,*999)
      CALL Decomposition_IsRegionDecomposition(decomposition,isRegionDecomposition,err,error,*999)
      IF(isRegionDecomposition) THEN
        !Create a new graph node for this decomposition
        NULLIFY(newGraphNode)
        CALL DecomposerGraphNode_Initialise(newGraphNode,err,error,*999)
        newGraphNode%decomposerGraph=>decomposerGraph
        newGraphNode%decomposition=>decomposition
        newGraphNode%rootNode=.TRUE.
        !Add this region decomposition to the graph roots
        ALLOCATE(newGraphRoots(decomposerGraph%numberOfGraphRoots+1),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate new graph roots.",err,error,*999)
        DO graphRootIdx=1,decomposerGraph%numberOfGraphRoots
          newGraphRoots(graphRootIdx)%ptr=>decomposerGraph%graphRoots(graphRootIdx)%ptr
        ENDDO !graphRootIdx
        newGraphRoots(decomposerGraph%numberOfGraphRoots+1)%ptr=>newGraphNode
        CALL MOVE_ALLOC(newGraphRoots,decomposerGraph%graphRoots)
        decomposerGraph%numberOfGraphRoots=decomposerGraph%numberOfGraphRoots+1
      ENDIF !isRegionDecomposition
    ENDDO !decompositionIdx (region)
    !Now loop over any interface decompositions
    DO decompositionIdx=1,decomposer%numberOfDecompositions
      NULLIFY(decomposition)
      CALL Decomposer_DecompositionGet(decomposer,decompositionIdx,decomposition,err,error,*999)
      CALL Decomposition_IsInterfaceDecomposition(decomposition,isInterfaceDecomposition,err,error,*999)
      IF(isInterfaceDecomposition) THEN
        NULLIFY(interface)
        CALL Decomposition_InterfaceGet(decomposition,interface,err,error,*999)
        !Find the meshes that are coupled in the interface in the decomposer graph
        ALLOCATE(matchedGraphNodes(interface%numberOfCoupledMeshes),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate matched graph nodes.",err,error,*999)
        numberOfMatchedGraphNodes=0
        numberOfGraphRoots=0
        DO coupledMeshIdx=1,INTERFACE%numberOfCoupledMeshes
          NULLIFY(coupledMesh)
          CALL Interface_CoupledMeshGet(INTERFACE,coupledMeshIdx,coupledMesh,err,error,*999)
          !Try and find decompositions of this coupled mesh in the decomposer graph.
          coupledMeshFound=.FALSE.
          DO graphRootIdx=1,decomposerGraph%numberOfGraphRoots
            NULLIFY(graphRoot)
            CALL DecomposerGraph_RootNodeGet(decomposerGraph,graphRootIdx,graphRoot,err,error,*999)
            !Look at linked nodes first
            DO graphLinkIdx=1,graphRoot%numberOfGraphLinks
              NULLIFY(graphLink)
              CALL DecomposerGraphNode_DecomposerGraphLinkGet(graphRoot,graphLinkIdx,graphLink,err,error,*999)              
              NULLIFY(linkedNode)
              CALL DecomposerGraphLink_LinkedNodeGet(graphLink,linkedNode,err,error,*999)
              NULLIFY(linkedNodeDecomposition)
              CALL DecomposerGraphNode_DecompositionGet(linkedNode,linkedNodeDecomposition,err,error,*999)
              NULLIFY(decompositionMesh)
              CALL Decomposition_MeshGet(linkedNodeDecomposition,decompositionMesh,err,error,*999)
              IF(ASSOCIATED(coupledMesh,decompositionMesh)) THEN
                numberOfMatchedGraphNodes=numberOfMatchedGraphNodes+1
                matchedGraphNodes(numberOfMatchedGraphNodes)%ptr=>linkedNode
                coupledMeshFound=.TRUE.
                IF(linkedNode%rootNode) numberOfGraphRoots=numberOfGraphRoots+1
                EXIT
              ENDIF
            ENDDO !graphLinkIdx
            !Now look at the root graph node
            NULLIFY(graphRootDecomposition)
            CALL DecomposerGraphNode_DecompositionGet(graphRoot,graphRootDecomposition,err,error,*999)
            NULLIFY(decompositionMesh)
            CALL Decomposition_MeshGet(graphRootDecomposition,decompositionMesh,err,error,*999)
            IF(ASSOCIATED(coupledMesh,decompositionMesh)) THEN
              numberOfMatchedGraphNodes=numberOfMatchedGraphNodes+1
              matchedGraphNodes(numberOfMatchedGraphNodes)%ptr=>graphRoot
              coupledMeshFound=.TRUE.
              IF(graphRoot%rootNode) numberOfGraphRoots=numberOfGraphRoots+1
              EXIT
            ENDIF
          ENDDO !graphRootIdx
          IF(.NOT.coupledMeshFound) THEN
            localError="Coupled mesh number "//TRIM(NumberToVString(coupledMesh%userNumber,"*",err,error))// &
              & " in interface number "//TRIM(NumberToVString(INTERFACE%userNumber,"*",err,error))// &
              & " could not be found in the decompositions added to decomposer number "// &
              & TRIM(NumberToVString(decomposer%userNumber,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ENDDO !coupledMeshIdx
        IF(numberOfMatchedGraphNodes<2) THEN
          localError="Invalid interface coupling. Only "//TRIM(NumberToVString(numberOfMatchedGraphNodes,"*",err,error))// &
            & " meshes were matched in the decomposer graph for interface number "// &
            & TRIM(NumberToVString(INTERFACE%userNumber,"*",err,error))//" and decomposer number "// &
            & TRIM(NumberToVString(decomposer%userNumber,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        graphNode=>matchedGraphNodes(1)%ptr
        DO matchedNodeIdx=2,numberOfMatchedGraphNodes
          linkedNode=>matchedGraphNodes(matchedNodeIdx)%ptr
          !Create the link
          NULLIFY(graphLink)
          CALL DecomposerGraphLink_Initialise(graphLink,err,error,*999)
          graphLink%parentGraphNode=>graphNode
          graphLink%linkedGraphNode=>linkedNode
          graphLink%decomposition=>decomposition
          ALLOCATE(newGraphLinks(graphNode%numberOfGraphLinks+1),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate new graph links.",err,error,*999)
          DO graphLinkIdx=1,graphNode%numberOfGraphLinks
            newGraphLinks(graphLinkIdx)%ptr => graphNode%graphLinks(graphLinkIdx)%ptr
          ENDDO !graphLinkIdx
          newGraphLinks(graphNode%numberOfGraphLinks+1)%ptr=>graphLink
          CALL MOVE_ALLOC(newGraphLinks,graphNode%graphLinks)
          graphNode%numberOfGraphLinks=graphNode%numberOfGraphLinks+1
          IF(linkedNode%rootNode) THEN
            !Remove this graph node from the list of graph roots.
            ALLOCATE(newGraphRoots(decomposerGraph%numberOfGraphRoots-1),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate new graph roots.",err,error,*999)
            numberOfGraphRoots=0
            DO graphRootIdx=1,decomposerGraph%numberOfGraphRoots
              NULLIFY(graphRoot)
              CALL DecomposerGraph_RootNodeGet(decomposerGraph,graphRootIdx,graphRoot,err,error,*999)
              IF(ASSOCIATED(graphRoot,linkedNode)) THEN
                linkedNode%rootNode=.FALSE.
              ELSE
                numberOfGraphRoots=numberOfGraphRoots+1
                newGraphRoots(numberOfGraphRoots)%ptr=>decomposerGraph%graphRoots(graphRootIdx)%ptr
              ENDIF
            ENDDO !graphRootIdx
            CALL MOVE_ALLOC(newGraphRoots,decomposerGraph%graphRoots)
            decomposerGraph%numberOfGraphRoots=numberOfGraphRoots
          ENDIF
        ENDDO !matchedNodeIdx
        IF(ALLOCATED(matchedGraphNodes)) DEALLOCATE(matchedGraphNodes)
      ENDIF !isInterfaceDecomposition
    ENDDO !decompositionIdx (interface)

    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Decomposer graph for decomposer number : ",decomposer%userNumber, &
        & err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of graph roots = ",decomposerGraph%numberOfGraphRoots, &
        & err,error,*999)
      DO graphRootIdx=1,decomposerGraph%numberOfGraphRoots
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Graph Root : ",graphRootIdx,err,error,*999)
        NULLIFY(graphRoot)
        CALL DecomposerGraph_RootNodeGet(decomposerGraph,graphRootIdx,graphRoot,err,error,*999)
        !Recursives output diagnostics for this graph Root
        CALL DecomposerGraphNode_Diagnostics(graphRoot,err,error,*999)        
      ENDDO !graphRootIdx
    ENDIF
       
    EXITS("Decomposer_GraphCalculate")
    RETURN
999 IF(ALLOCATED(newGraphRoots)) DEALLOCATE(newGraphRoots)
    IF(ALLOCATED(matchedGraphNodes)) DEALLOCATE(matchedGraphNodes)
    IF(ALLOCATED(newGraphLinks)) DEALLOCATE(newGraphLinks)
    ERRORSEXITS("Decomposer_GraphCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposer_GraphCalculate

  !
  !================================================================================================================================
  !

  !>Finalise the decomposer graph and deallocate all memory.
  SUBROUTINE Decomposer_GraphFinalise(decomposerGraph,err,error,*)

    !Argument variables
    TYPE(DecomposerGraphType), POINTER :: decomposerGraph !<A pointer to the decomposer graph to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: graphRootIdx

    ENTERS("Decomposer_GraphFinalise",err,error,*999)

    IF(ASSOCIATED(decomposerGraph)) THEN
      IF(ALLOCATED(decomposerGraph%graphRoots)) THEN
        DO graphRootIdx=1,SIZE(decomposerGraph%GraphRoots,1)
          CALL DecomposerGraphNode_Finalise(decomposerGraph%graphRoots(graphRootIdx)%ptr,err,error,*999)
        ENDDO !graphRootIdx
        DEALLOCATE(decomposerGraph%graphRoots)
      ENDIF
      DEALLOCATE(decomposerGraph)
    ENDIF
       
    EXITS("Decomposer_GraphFinalise")
    RETURN
999 ERRORSEXITS("Decomposer_GraphFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposer_GraphFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the decomposer graph for a decomposer.
  SUBROUTINE Decomposer_GraphInitialise(decomposer,err,error,*)

    !Argument variables
    TYPE(DecomposerType), POINTER :: decomposer !<A pointer to the decomposer to initialise the graph for. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
 
    ENTERS("Decomposer_GraphInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(decomposer)) CALL FlagError("Decomposer is not associated.",err,error,*998)
    IF(ASSOCIATED(decomposer%decomposerGraph)) CALL FlagError("Decomposer graph is already associated.",err,error,*998)

    ALLOCATE(decomposer%decomposerGraph,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate decomposer graph.",err,error,*999)
    decomposer%decomposerGraph%decomposer=>decomposer
    decomposer%decomposerGraph%numberOfGraphRoots=0
        
    EXITS("Decomposer_GraphInitialise")
    RETURN
999 CALL Decomposer_GraphFinalise(decomposer%decomposerGraph,dummyErr,dummyError,*998)
998 ERRORSEXITS("Decomposer_GraphInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposer_GraphInitialise
  
  !
  !================================================================================================================================
  !

  !>Finalise the decomposer and deallocate all memory.
  SUBROUTINE Decomposer_Finalise(decomposer,err,error,*)

    !Argument variables
    TYPE(DecomposerType), POINTER :: decomposer !<A pointer to the decomposer to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: decompositionIdx

    ENTERS("Decomposer_Finalise",err,error,*999)

    IF(ASSOCIATED(decomposer)) THEN
      IF(ALLOCATED(decomposer%decompositions)) THEN
        DO decompositionIdx=1,SIZE(decomposer%decompositions,1)
          IF(ASSOCIATED(decomposer%decompositions(decompositionIdx)%ptr)) THEN
            NULLIFY(decomposer%decompositions(decompositionIdx)%ptr%decomposer)
          ENDIF
        ENDDO !decompositionIdx
        DEALLOCATE(decomposer%decompositions)
      ENDIF
      CALL Decomposer_GraphFinalise(decomposer%decomposerGraph,err,error,*999)
      DEALLOCATE(decomposer)
    ENDIF
       
    EXITS("Decomposer_Finalise")
    RETURN
999 ERRORSEXITS("Decomposer_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposer_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises a decomposer.
  SUBROUTINE Decomposer_Initialise(decomposer,err,error,*)

    !Argument variables
    TYPE(DecomposerType), POINTER :: decomposer !<A pointer to the decomposer to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
 
    ENTERS("Decomposer_Initialise",err,error,*998)

    IF(ASSOCIATED(decomposer)) CALL FlagError("Decomposer is already associated.",err,error,*998)

    ALLOCATE(decomposer,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate decomposer.",err,error,*999)
    decomposer%userNumber=0
    NULLIFY(decomposer%decomposers)
    NULLIFY(decomposer%region)
    decomposer%decomposerFinished=.FALSE.
    decomposer%numberOfDecompositions=0
    NULLIFY(decomposer%decomposerGraph)
    NULLIFY(decomposer%workGroup)
        
    EXITS("Decomposer_Initialise")
    RETURN
999 CALL Decomposer_Finalise(decomposer,dummyErr,dummyError,*998)
998 ERRORSEXITS("Decomposer_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposer_Initialise
  
  !
  !================================================================================================================================
  !

  !>Initialises a decomposer.
  SUBROUTINE Decomposer_TopologyCalculate(decomposer,err,error,*)

    !Argument variables
    TYPE(DecomposerType), POINTER :: decomposer !<A pointer to the decomposer to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: decompositionIdx
    TYPE(DecompositionType), POINTER :: decomposition
 
    ENTERS("Decomposer_TopologyCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(decomposer)) CALL FlagError("Decomposer is not associated.",err,error,*999)

    DO decompositionIdx=1,decomposer%numberOfDecompositions
      NULLIFY(decomposition)
      CALL Decomposer_DecompositionGet(decomposer,decompositionIdx,decomposition,err,error,*999)
      !Initialise the topology information for this decomposition
      CALL Decomposition_TopologyInitialise(decomposition,err,error,*999)
      !Initialise the domains for this decomposition
      CALL Decomposition_DomainsInitialise(decomposition,err,error,*999)
      !Calculate the decomposition topology
      CALL Decomposition_TopologyCalculate(decomposition,err,error,*999)
    ENDDO !decompositionIdx
         
    EXITS("Decomposer_TopologyCalculate")
    RETURN
999 ERRORSEXITS("Decomposer_TopologyCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposer_TopologyCalculate
  
  !
  !================================================================================================================================
  !

  !>Initialises a decomposer graph link.
  SUBROUTINE DecomposerGraphLink_Initialise(decomposerGraphLink,err,error,*)

    !Argument variables
    TYPE(DecomposerGraphLinkType), POINTER :: decomposerGraphLink !<A pointer to the decomposer graph link to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
 
    ENTERS("DecomposerGraphLink_Initialise",err,error,*998)

    IF(ASSOCIATED(decomposerGraphLink)) CALL FlagError("Decomposer graph link is already associated.",err,error,*998)

    ALLOCATE(decomposerGraphLink,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate decomposer graph link.",err,error,*999)
    NULLIFY(decomposerGraphLink%parentGraphNode)
    NULLIFY(decomposerGraphLink%linkedGraphNode)
    NULLIFY(decomposerGraphLink%decomposition)
        
    EXITS("DecomposerGraphLink_Initialise")
    RETURN
999 CALL DecomposerGraphLink_Finalise(decomposerGraphLink,dummyErr,dummyError,*998)
998 ERRORSEXITS("DecomposerGraphLink_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE DecomposerGraphLink_Initialise
  
  !
  !================================================================================================================================
  !

  !>Finalise the decomposer graph link and deallocate all memory.
  SUBROUTINE DecomposerGraphLink_Finalise(decomposerGraphLink,err,error,*)

    !Argument variables
    TYPE(DecomposerGraphLinkType), POINTER :: decomposerGraphLink !<A pointer to the decomposer graph link to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DecomposerGraphLink_Finalise",err,error,*999)

    IF(ASSOCIATED(decomposerGraphLink)) THEN
      DEALLOCATE(decomposerGraphLink)
    ENDIF
       
    EXITS("DecomposerGraphLink_Finalise")
    RETURN
999 ERRORSEXITS("DecomposerGraphLink_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE DecomposerGraphLink_Finalise

  !
  !================================================================================================================================
  !

  !>Outputs diagnostics for a decomposer graph node.
  RECURSIVE SUBROUTINE DecomposerGraphNode_Diagnostics(graphNode,err,error,*)

    !Argument variables
    TYPE(DecomposerGraphNodeType), POINTER :: graphNode !<A pointer to the decomposer graph node to output the diagnostics for. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: graphLinkIdx
    TYPE(DecomposerGraphLinkType), POINTER :: graphLink
    TYPE(DecomposerGraphNodeTYpe), POINTER :: linkedNode
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(InterfaceType), POINTER :: interface
    TYPE(MeshType), POINTER :: mesh
    TYPE(RegionType), POINTER :: region
 
    ENTERS("DecomposerGraphNode_Diagnostics",err,error,*999)

    CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"",err,error,*999)
    CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Decomposer graph node information:",err,error,*999)
    NULLIFY(decomposition)
    CALL DecomposerGraphNode_DecompositionGet(graphNode,decomposition,err,error,*999)
    CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Decomposition user number = ",decomposition%userNumber,err,error,*999)
    NULLIFY(mesh)
    CALL Decomposition_MeshGet(decomposition,mesh,err,error,*999)
    CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Mesh user number = ",mesh%userNumber,err,error,*999)
    NULLIFY(region)
    CALL Mesh_RegionGet(mesh,region,err,error,*999)
    CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Region user number = ",region%userNumber,err,error,*999)
    CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of graph links = ",graphNode%numberOfGraphLinks,err,error,*999)
    DO graphLinkIdx=1,graphNode%numberOfGraphLinks
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Graph link : ",graphLinkIdx,err,error,*999)
      NULLIFY(graphLink)
      CALL DecomposerGraphNode_DecomposerGraphLinkGet(graphNode,graphLinkIdx,graphLink,err,error,*999)
      NULLIFY(decomposition)
      CALL DecomposerGraphLink_DecompositionGet(graphLink,decomposition,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Decomposition user number = ",decomposition%userNumber,err,error,*999)
      NULLIFY(mesh)
      CALL Decomposition_MeshGet(decomposition,mesh,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Mesh user number = ",mesh%userNumber,err,error,*999)
      NULLIFY(interface)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Interface user number = ",interface%userNumber,err,error,*999)
      NULLIFY(region)
      CALL Interface_ParentRegionGet(interface,region,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"            Parent region user number = ",region%userNumber,err,error,*999)
      NULLIFY(linkedNode)
      CALL DecomposerGraphLink_LinkedNodeGet(graphLink,linkedNode,err,error,*999)
      CALL DecomposerGraphNode_Diagnostics(linkedNode,err,error,*999)
    ENDDO !graphLinkIdx
        
    EXITS("DecomposerGraphNode_Diagnostics")
    RETURN
999 ERRORSEXITS("DecomposerGraphNode_Diagnostics",err,error)
    RETURN 1
    
  END SUBROUTINE DecomposerGraphNode_Diagnostics
  
  !
  !================================================================================================================================
  !

  !>Finalise the decomposer graph node and deallocate all memory.
  RECURSIVE SUBROUTINE DecomposerGraphNode_Finalise(decomposerGraphNode,err,error,*)

    !Argument variables
    TYPE(DecomposerGraphNodeType), POINTER :: decomposerGraphNode !<A pointer to the decomposer graph node to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: graphLinkIdx
    TYPE(DecomposerGraphNodeType), POINTER :: linkedNode

    ENTERS("DecomposerGraphNode_Finalise",err,error,*999)

    IF(ASSOCIATED(decomposerGraphNode)) THEN
      IF(ALLOCATED(decomposerGraphNode%graphLinks)) THEN
        DO graphLinkIdx=1,SIZE(decomposerGraphNode%graphLinks,1)
          linkedNode => decomposerGraphNode%graphLinks(graphLinkIdx)%ptr%linkedGraphNode
          CALL DecomposerGraphNode_Finalise(linkedNode,err,error,*999)
          CALL DecomposerGraphLink_Finalise(decomposerGraphNode%graphLinks(graphLinkIdx)%ptr,err,error,*999)
        ENDDO !graphLinkIdx
        DEALLOCATE(decomposerGraphNode%graphLinks)
      ENDIF
      DEALLOCATE(decomposerGraphNode)
    ENDIF
       
    EXITS("DecomposerGraphNode_Finalise")
    RETURN
999 ERRORSEXITS("DecomposerGraphNode_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE DecomposerGraphNode_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises a decomposer graph node.
  SUBROUTINE DecomposerGraphNode_Initialise(decomposerGraphNode,err,error,*)

    !Argument variables
    TYPE(DecomposerGraphNodeType), POINTER :: decomposerGraphNode !<A pointer to the decomposer graph node to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
 
    ENTERS("DecomposerGraphNode_Initialise",err,error,*998)

    IF(ASSOCIATED(decomposerGraphNode)) CALL FlagError("Decomposer graph node is already associated.",err,error,*998)

    ALLOCATE(decomposerGraphNode,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate decomposer graph node.",err,error,*999)
    NULLIFY(decomposerGraphNode%decomposerGraph)
    NULLIFY(decomposerGraphNode%decomposition)
    decomposerGraphNode%rootNode=.FALSE.
    decomposerGraphNode%numberOfGraphLinks=0
        
    EXITS("DecomposerGraphNode_Initialise")
    RETURN
999 CALL DecomposerGraphNode_Finalise(decomposerGraphNode,dummyErr,dummyError,*998)
998 ERRORSEXITS("DecomposerGraphNode_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE DecomposerGraphNode_Initialise
  
  !
  !================================================================================================================================
  !

  !>Finalise the decomposers and deallocate all memory.
  SUBROUTINE Decomposers_Finalise(decomposers,err,error,*)

    !Argument variables
    TYPE(DecomposersType), POINTER :: decomposers !<A pointer to the decomposers to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: decomposerIdx

    ENTERS("Decomposers_Finalise",err,error,*999)

    IF(ASSOCIATED(decomposers)) THEN
      DO decomposerIdx=1,SIZE(decomposers%decomposers,1)
        CALL Decomposer_Finalise(decomposers%decomposers(decomposeridx)%ptr,err,error,*999)
      ENDDO !decomposerIdx
      DEALLOCATE(decomposers)
    ENDIF
       
    EXITS("Decomposers_Finalise")
    RETURN
999 ERRORSEXITS("Decomposers_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposers_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the decomposers for a region.
  SUBROUTINE Decomposers_Initialise(region,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: region !<A pointer to the region to initialise the decomposers for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
 
    ENTERS("Decomposers_Initialise",err,error,*998)

    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated.",err,error,*998)
    IF(ASSOCIATED(region%decomposers)) CALL FlagError("Decomposers is already associated for this region.",err,error,*998)

    ALLOCATE(region%decomposers,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate decomposers.",err,error,*999)
    region%decomposers%region=>region
    region%decomposers%numberOfDecomposers=0
       
    EXITS("Decomposers_Initialise")
    RETURN
999 CALL Decomposers_Finalise(region%decomposers,dummyErr,dummyError,*998)
998 ERRORSEXITS("Decomposers_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposers_Initialise
  
  !
  !================================================================================================================================
  !

  !>Finialises a decomposition adjacent element information and deallocates all memory
  SUBROUTINE Decomposition_AdjacentElementFinalise(decompositionAdjacentElement,err,error,*)
    
    !Argument variables
    TYPE(DecompositionAdjacentElementType) :: decompositionAdjacentElement !<The decomposition adjacent element to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Decomposition_AdjacentElementFinalise",err,error,*999)

    decompositionAdjacentElement%numberOfAdjacentElements=0
    IF(ALLOCATED(decompositionAdjacentElement%adjacentElements)) DEALLOCATE(decompositionAdjacentElement%adjacentElements)
       
    EXITS("Decomposition_AdjacentElementFinalise")
    RETURN
999 ERRORSEXITS("Decomposition_AdjacentElementFinalise",err,error)    
    RETURN 1
   
  END SUBROUTINE Decomposition_AdjacentElementFinalise

  !
  !================================================================================================================================
  !
  !>Initalises a decomposition adjacent element information.
  SUBROUTINE Decomposition_AdjacentElementInitialise(decompositionAdjacentElement,err,error,*)
    
    !Argument variables
    TYPE(DecompositionAdjacentElementType) :: decompositionAdjacentElement !<The decomposition adjacent element to initialise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Decomposition_AdjacentElementInitialise",err,error,*999)

    decompositionAdjacentElement%numberOfAdjacentElements=0
       
    EXITS("Decomposition_AdjacentElementInitialise")
    RETURN
999 ERRORSEXITS("Decomposition_AdjacentElementInitialise",err,error)    
    RETURN 1
   
  END SUBROUTINE Decomposition_AdjacentElementInitialise

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a domain decomposition on a given mesh. \see OpenCMISS::Iron::cmfe_Decomposition_CreateFinish
  SUBROUTINE Decomposition_CreateFinish(decomposition,err,error,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer to the decomposition to finish creating
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: decompositionIdx
    TYPE(DecompositionsType), POINTER :: decompositions
    TYPE(MeshType), POINTER :: mesh

    ENTERS("Decomposition_CreateFinish",err,error,*999)

    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)
    
    decomposition%decompositionFinished=.TRUE.
    
    IF(diagnostics1) THEN
      NULLIFY(mesh)
      CALL Decomposition_MeshGet(decomposition,mesh,err,error,*999)
      NULLIFY(decompositions)
      CALL Decomposition_DecompositionsGet(decomposition,decompositions,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Mesh = ",mesh%userNumber,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of decompositions = ",decompositions%numberOfDecompositions, &
        & err,error,*999)
      DO decompositionIdx=1,decompositions%numberOfDecompositions
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Decomposition number = ",decompositionIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Global number        = ", &
          & decompositions%decompositions(decompositionIdx)%ptr%globalNumber,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    User number          = ", &
          & decompositions%decompositions(decompositionIdx)%ptr%userNumber,err,error,*999)
      ENDDO !decompositionIdx
    ENDIF
    
    EXITS("Decomposition_CreateFinish")
    RETURN
999 ERRORSEXITS("Decomposition_CreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_CreateFinish

  !
  !================================================================================================================================
  !

  !>Starts the creation of a domain decomposition for a given mesh. \see OpenCMISS::Iron::cmfe_Decomposition_CreateStart
  SUBROUTINE Decomposition_CreateStart(userNumber,mesh,decomposition,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the decomposition 
    TYPE(MeshType), POINTER :: mesh !<A pointer to the mesh to decompose
    TYPE(DecompositionType), POINTER :: decomposition !<On return, a pointer to the created decomposition. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: decompositionIdx,dummyErr
    TYPE(ComputationEnvironmentType), POINTER :: computationEnvironment
    TYPE(ContextType), POINTER :: context
    TYPE(DecompositionType), POINTER :: newDecomposition    
    TYPE(DecompositionPtrType), ALLOCATABLE :: newDecompositions(:)
    TYPE(DecompositionsType), POINTER :: decompositions
    TYPE(RegionType), POINTER :: region    
    TYPE(VARYING_STRING) :: dummyError,localError
    TYPE(WorkGroupType), POINTER :: worldWorkGroup

    NULLIFY(newDecomposition)

    ENTERS("Decomposition_CreateStart",err,error,*999)
    
    IF(ASSOCIATED(decomposition)) CALL FlagError("Decomposition is already associated.",err,error,*999)
    CALL Mesh_AssertIsFinished(mesh,err,error,*999)
    NULLIFY(decomposition)
    CALL Decomposition_UserNumberFind(userNumber,mesh,decomposition,err,error,*999)
    IF(ASSOCIATED(decomposition)) THEN
      localError="Decomposition number "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " has already been created on mesh number "//TRIM(NumberToVString(mesh%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    NULLIFY(decompositions)
    CALL Mesh_DecompositionsGet(mesh,decompositions,err,error,*999)

    NULLIFY(newDecomposition)
    CALL Decomposition_Initialise(newDecomposition,err,error,*999)
    !Set default decomposition properties              
    newDecomposition%globalNumber=decompositions%numberOfDecompositions+1
    newDecomposition%userNumber=userNumber
    newDecomposition%decompositions=>decompositions
    newDecomposition%mesh=>mesh
    newDecomposition%region=>mesh%region
    newDecomposition%interface=>mesh%interface
    newDecomposition%numberOfDimensions=mesh%numberOfDimensions
    newDecomposition%numberOfComponents=mesh%numberOfComponents
    !By default, the process of decompostion was done on the first mesh components. But the decomposition is the
    !same for all mesh components, since the decomposition is element-based.
    newDecomposition%meshComponentNumber=1
    !Default decomposition is all the mesh with one domain.
    newDecomposition%domainDecompositionType=DECOMPOSITION_ALL_TYPE
    NULLIFY(region)
    CALL Mesh_RegionGet(mesh,region,err,error,*999)
    NULLIFY(context)
    CALL Region_ContextGet(region,context,err,error,*999)
    NULLIFY(computationEnvironment)
    CALL Context_ComputationEnvironmentGet(context,computationEnvironment,err,error,*999)
    NULLIFY(worldWorkGroup)
    CALL ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup,err,error,*999)
    newDecomposition%workGroup=>worldWorkGroup
    newDecomposition%numberOfDomains=1
    newDecomposition%numberOfElements=mesh%numberOfElements
    ALLOCATE(newDecomposition%elementDomain(mesh%numberOfElements),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new decomposition element domain.",err,error,*999)
    newDecomposition%elementDomain=0          
    !Add new decomposition into list of decompositions on the mesh
    ALLOCATE(newDecompositions(mesh%decompositions%numberOfDecompositions+1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new decompositions.",err,error,*999)
    DO decompositionIdx=1,decompositions%numberOfDecompositions
      newDecompositions(decompositionIdx)%ptr=>decompositions%decompositions(decompositionIdx)%ptr
    ENDDO !decompositionIdx
    newDecompositions(decompositions%numberOfDecompositions+1)%ptr=>newDecomposition
    CALL MOVE_ALLOC(newDecompositions,decompositions%decompositions)
    decompositions%numberOfDecompositions=decompositions%numberOfDecompositions+1        
    decomposition=>newDecomposition
   
    EXITS("Decomposition_CreateStart")
    RETURN
999 CALL Decomposition_Finalise(newDecomposition,dummyErr,dummyError,*998)
998 IF(ALLOCATED(newDecompositions)) DEALLOCATE(newDecompositions)
    ERRORSEXITS("Decomposition_CreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_CreateStart

  !
  !================================================================================================================================
  !

  !>Destroys a domain decomposition identified by an object and deallocates all memory. \see OpenCMISS::Iron::cmfe_Decomposition_Destroy
  SUBROUTINE Decomposition_Destroy(decomposition,err,error,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: decomposition !<The decomposition to destroy.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: decompositionIdx,decompositionPosition
    LOGICAL :: found
    TYPE(DecompositionPtrType), ALLOCATABLE :: newDecompositions(:)
    TYPE(DecompositionsType), POINTER :: decompositions
    TYPE(VARYING_STRING) :: localError

    ENTERS("Decomposition_Destroy",err,error,*999)

    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)
    
    NULLIFY(decompositions)
    CALL Decomposition_DecompositionsGet(decomposition,decompositions,err,error,*999)
    !Find the decomposition identified by the user number
    found=.FALSE.
    decompositionPosition=0
    DO WHILE(decompositionPosition<decompositions%numberOfDecompositions.AND..NOT.found)
      decompositionPosition=decompositionPosition+1
      IF(decompositions%decompositions(decompositionPosition)%ptr%userNumber==decomposition%userNumber) found=.TRUE.
    ENDDO
    IF(.NOT.found) THEN
      localError="Decomposition number "//TRIM(NumberToVString(decomposition%userNumber,"*",err,error))// &
        & " could not be found in the list of decompositions."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    !Finalise the decomposition
    CALL Decomposition_Finalise(decomposition,err,error,*999)
          
    !Remove the decomposition from the list of decompositions
    IF(decompositions%numberOfDecompositions>1) THEN
      ALLOCATE(newDecompositions(decompositions%numberOfDecompositions-1),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new decompositions.",err,error,*999)
      DO decompositionIdx=1,decompositions%numberOfDecompositions
        IF(decompositionIdx<decompositionPosition) THEN
          newDecompositions(decompositionIdx)%ptr=>decompositions%decompositions(decompositionIdx)%ptr
        ELSE IF(decompositionIdx>decompositionPosition) THEN
          decompositions%decompositions(decompositionIdx)%ptr%globalNumber= &
            & decompositions%decompositions(decompositionIdx)%ptr%globalNumber-1
          newDecompositions(decompositionIdx-1)%ptr=>decompositions%decompositions(decompositionIdx)%ptr
        ENDIF
      ENDDO !decompositionIdx
      CALL MOVE_ALLOC(newDecompositions,decompositions%decompositions)
      decompositions%numberOfDecompositions=decompositions%numberOfDecompositions-1
    ELSE
      DEALLOCATE(decompositions%decompositions)
      decompositions%numberOfDecompositions=0
    ENDIF
    
    EXITS("Decomposition_Destroy")
    RETURN
999 IF(ALLOCATED(newDecompositions)) DEALLOCATE(newDecompositions)
    ERRORSEXITS("Decomposition_Destroy",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_Destroy

  !
  !================================================================================================================================
  !

  !>Get the basis for an element in the decomposition
  SUBROUTINE Decomposition_ElementBasisGet(decomposition,meshComponentNumber,elementUserNumber,basis,err,error,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer to the decomposition to get the element basis for
    INTEGER(INTG), INTENT(IN) :: meshComponentNumber !<The mesh component number to get the element basis for
    INTEGER(INTG), INTENT(IN) :: elementUserNumber !<The element user number to get the basis for
    TYPE(BasisType), POINTER, INTENT(OUT) :: basis !<On return, a pointer to the basis for the element. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    LOGICAL :: ghostElement
    INTEGER(INTG) :: elementLocalNumber
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainTopologyType), POINTER :: domainTopology

    ENTERS("Decomposition_ElementBasisGet",err,error,*998)

    IF(ASSOCIATED(basis)) CALL FlagError("Basis is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*998)

    NULLIFY(decompositionTopology)
    CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*999)
    NULLIFY(decompositionElements)
    CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology,decompositionElements,err,error,*999)
    CALL DecompositionElements_LocalElementNumberGet(decompositionElements,elementUserNumber,elementLocalNumber, &
      & ghostElement,err,error,*999)
    NULLIFY(domain)
    CALL Decomposition_DomainGet(decomposition,meshComponentNumber,domain,err,error,*999)
    NULLIFY(domainTopology)
    CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
    NULLIFY(domainElements)
    CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
    CALL DomainElements_BasisGet(domainElements,elementLocalNumber,basis,err,error,*999)
   
    EXITS("Decomposition_ElementBasisGet")
    RETURN
999 NULLIFY(basis)
998 ERRORSEXITS("Decomposition_ElementBasisGet",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_ElementBasisGet

  !
  !================================================================================================================================
  !

  !>Calculates the element domains for a decomposition of a mesh. \see OpenCMISS::Iron::cmfe_Decomposition_ElementDomainCalculate
  SUBROUTINE Decomposition_ElementDomainCalculate(decomposition,err,error,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer to the decomposition to calculate the element domains for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfElementIndices,elementIndex,elementCount,elementIdx,localNodeIdx,myComputationNodeNumber, &
      & numberOfComputationNodes,computationNodeIdx,elementStart,elementStop,myElementStart, &
      & myElementStop,numberOfElements,myNumberOfElements,mpiIError,maxNumberOfElementsPerNode,componentIdx,minNumberOfXi
    INTEGER(INTG), ALLOCATABLE :: elementCounts(:),elementPtr(:),elementIndices(:),elementDistance(:),displacements(:), &
      & receiveCounts(:)
    INTEGER(INTG) :: elementWeight(1),weightFlag,numberFlag,numberOfConstraints, &
      & numberOfCommonNodes,parmetisOptions(0:2),randomSeedsSize
    INTEGER(INTG) :: groupCommunicator
    INTEGER(INTG), ALLOCATABLE :: randomSeeds(:)
    REAL(DP) :: ubvec(1)
    REAL(DP), ALLOCATABLE :: tpwgts(:)
    REAL(DP) :: numberOfElementsPerNode
    TYPE(ContextType), POINTER :: context    
    TYPE(BasisType), POINTER :: basis
    TYPE(MeshType), POINTER :: mesh
    TYPE(VARYING_STRING) :: localError

    ENTERS("Decomposition_ElementDomainCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)
    
    NULLIFY(mesh)
    CALL Decomposition_MeshGet(decomposition,mesh,err,error,*999)
 
    componentIdx=decomposition%meshComponentNumber
    
    CALL WorkGroup_GroupCommunicatorGet(decomposition%workGroup,groupCommunicator,err,error,*999)
    CALL WorkGroup_NumberOfGroupNodesGet(decomposition%workGroup,numberOfComputationNodes,err,error,*999)
    CALL WorkGroup_GroupNodeNumberGet(decomposition%workGroup,myComputationNodeNumber,err,error,*999)
    
    SELECT CASE(decomposition%domainDecompositionType)
    CASE(DECOMPOSITION_ALL_TYPE)
      !Do nothing. Decomposition checked below.
    CASE(DECOMPOSITION_CALCULATED_TYPE)
      !Calculate the general decomposition

      IF(decomposition%numberOfDomains==1) THEN
        decomposition%elementDomain=0
      ELSE              
        numberOfElementsPerNode=REAL(mesh%numberOfElements,DP)/REAL(numberOfComputationNodes,DP)
        elementStart=1
        elementStop=0
        maxNumberOfElementsPerNode=-1
        ALLOCATE(receiveCounts(0:numberOfComputationNodes-1),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate recieve counts.",err,error,*999)
        ALLOCATE(displacements(0:numberOfComputationNodes-1),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate displacements.",err,error,*999)
        ALLOCATE(elementDistance(0:numberOfComputationNodes),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate element distance.",err,error,*999)
        elementDistance(0)=0
        DO computationNodeIdx=0,numberOfComputationNodes-1
          elementStart=elementStop+1
          IF(computationNodeIdx==numberOfComputationNodes-1) THEN
            elementStop=mesh%numberOfElements
          ELSE
            elementStop=elementStart+NINT(numberOfElementsPerNode,INTG)-1
          ENDIF
          IF((numberOfComputationNodes-1-computationNodeIdx)>(mesh%numberOfElements-elementStop)) &
            & elementStop=mesh%numberOfElements-(numberOfComputationNodes-1-computationNodeIdx)
          IF(elementStart>mesh%numberOfElements) elementStart=mesh%numberOfElements
          IF(elementStop>mesh%numberOfElements) elementStop=mesh%numberOfElements
          displacements(computationNodeIdx)=elementStart-1
          elementDistance(computationNodeIdx+1)=elementStop !C numbering
          numberOfElements=elementStop-elementStart+1
          receiveCounts(computationNodeIdx)=numberOfElements
          IF(numberOfElements>maxNumberOfElementsPerNode) maxNumberOfElementsPerNode=numberOfElements
          IF(computationNodeIdx==myComputationNodeNumber) THEN
            myElementStart=elementStart
            myElementStop=elementStop
            myNumberOfElements=elementStop-elementStart+1
            numberOfElementIndices=0
            DO elementIdx=myElementStart,myElementStop
              basis=>mesh%topology(componentIdx)%ptr%elements%elements(elementIdx)%basis
              numberOfElementIndices=numberOfElementIndices+basis%numberOfNodes
            ENDDO !elementIdx
          ENDIF
        ENDDO !computationNodeIdx
        
        ALLOCATE(elementPtr(0:myNumberOfElements),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate element pointer list.",err,error,*999)
        ALLOCATE(elementIndices(0:numberOfElementIndices-1),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate element indicies list.",err,error,*999)
        ALLOCATE(tpwgts(1:decomposition%numberOfDomains),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate tpwgts.",err,error,*999)
        elementIndex=0
        elementCount=0
        elementPtr(0)=0
        minNumberOfXi=99999
        DO elementIdx=myElementStart,myElementStop
          elementCount=elementCount+1
          basis=>mesh%topology(componentIdx)%ptr%elements%elements(elementIdx)%basis
          IF(basis%numberOfXi<minNumberOfXi) minNumberOfXi=basis%numberOfXi
          DO localNodeIdx=1,basis%numberOfNodes
            elementIndices(elementIndex)=mesh%topology(componentIdx)%ptr%elements%elements(elementIdx)% &
              & meshElementNodes(localNodeIdx)-1 !C numbering
            elementIndex=elementIndex+1
          ENDDO !localNodeIdx
          elementPtr(elementCount)=elementIndex !C numbering
        ENDDO !elementIdx
              
        !Set up ParMETIS variables
        NULLIFY(context)
        CALL WorkGroup_ContextGet(decomposition%workGroup,context,err,error,*999)
        CALL Context_RandomSeedsSizeGet(context,randomSeedsSize,err,error,*999)
        ALLOCATE(randomSeeds(randomSeedsSize),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate random seeds.",err,error,*999)
        CALL Context_RandomSeedsGet(context,randomSeeds,err,error,*999)              
        weightFlag=0 !No weights
        elementWeight(1)=1 !Isn't used due to weight flag
        numberFlag=0 !C Numbering as there is a bug with Fortran numbering              
        numberOfConstraints=1
        IF(minNumberOfXi==1) THEN
          numberOfCommonNodes=1
        ELSE
          numberOfCommonNodes=2
        ENDIF
        tpwgts=1.0_DP/REAL(decomposition%numberOfDomains,DP)
        ubvec=1.05_DP
        parmetisOptions(0)=1 !If zero, defaults are used, otherwise next two values are used
        parmetisOptions(1)=7 !Level of information to output
        parmetisOptions(2)=randomSeeds(1) !Seed for random number generator
        IF(ALLOCATED(randomSeeds)) DEALLOCATE(randomSeeds)
        !Call ParMETIS to calculate the partitioning of the mesh graph.
        CALL PARMETIS_PARTMESHKWAY(elementDistance,elementPtr,elementIndices,elementWeight,weightFlag,numberFlag, &
          & numberOfConstraints,numberOfCommonNodes,decomposition%numberOfDomains,tpwgts,ubvec,parmetisOptions, &
          & decomposition%numberOfEdgesCut,decomposition%elementDomain(displacements(myComputationNodeNumber)+1:), &
          & groupCommunicator,err,error,*999)
              
        !Transfer all the element domain information to the other computation nodes so that each rank has all the info
        IF(numberOfComputationNodes>1) THEN
          !This should work on a single processor but doesn't for mpich2 under windows. Maybe a bug? Avoid for now.
          CALL MPI_ALLGATHERV(MPI_IN_PLACE,maxNumberOfElementsPerNode,MPI_INTEGER,decomposition%elementDomain, &
            & receiveCounts,displacements,MPI_INTEGER,groupCommunicator,mpiIError)
          CALL MPI_ErrorCheck("MPI_ALLGATHERV",mpiIError,err,error,*999)
        ENDIF
              
        DEALLOCATE(displacements)
        DEALLOCATE(receiveCounts)
        DEALLOCATE(elementDistance)
        DEALLOCATE(elementPtr)
        DEALLOCATE(elementIndices)
        DEALLOCATE(tpwgts)

      ENDIF
            
    CASE(DECOMPOSITION_USER_DEFINED_TYPE)
      !Do nothing. Decomposition checked below.          
    CASE DEFAULT
      CALL FlagError("Invalid domain decomposition type.",err,error,*999)            
    END SELECT
    
    !Check decomposition and check that each domain has an element in it.
    ALLOCATE(elementCounts(0:numberOfComputationNodes-1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate element count.",err,error,*999)
    elementCounts=0
    DO elementIndex=1,mesh%numberOfElements
      computationNodeIdx=decomposition%elementDomain(elementIndex)
      IF(computationNodeIdx>=0.AND.computationNodeIdx<numberOfComputationNodes) THEN
        elementCounts(computationNodeIdx)=elementCounts(computationNodeIdx)+1
      ELSE
        localError="The computation node number of "//TRIM(NumberToVString(computationNodeIdx,"*",err,error))// &
          & " for element number "//TRIM(NumberToVString(elementIndex,"*",err,error))// &
          & " is invalid. The computation node number must be between 0 and "// &
          & TRIM(NumberToVString(numberOfComputationNodes-1,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDDO !elementIndex
    DO computationNodeIdx=0,numberOfComputationNodes-1
      IF(elementCounts(computationNodeIdx)==0) THEN
        localError="Invalid decomposition. There are no elements in computation node "// &
          & TRIM(NumberToVString(computationNodeIdx,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDDO !computationNodeIdx
    DEALLOCATE(elementCounts)
             
    IF(diagnostics1) THEN
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Decomposition for mesh number ",mesh%userNumber,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of domains = ",decomposition%numberOfDomains,err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Element domains:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Decomposition type = ", decomposition%domainDecompositionType, &
        & err,error,*999)
      IF(decomposition%domainDecompositionType==DECOMPOSITION_CALCULATED_TYPE) THEN
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of edges cut = ",decomposition%numberOfEdgesCut, &
          & err,error,*999)
      ENDIF
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of elements = ",decomposition%numberOfElements, &
        & err,error,*999)
      DO elementIdx=1,decomposition%numberOfElements
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Element = ",elementIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Domain = ",decomposition%elementDomain(elementIdx), &
          & err,error,*999)
      ENDDO !elementIdx
    ENDIF
    
    EXITS("Decomposition_ElementDomainCalculate")
    RETURN
999 IF(ALLOCATED(receiveCounts)) DEALLOCATE(receiveCounts)
    IF(ALLOCATED(displacements)) DEALLOCATE(displacements)
    IF(ALLOCATED(elementDistance)) DEALLOCATE(elementDistance)
    IF(ALLOCATED(elementPtr)) DEALLOCATE(elementPtr)
    IF(ALLOCATED(elementIndices)) DEALLOCATE(elementIndices)
    IF(ALLOCATED(tpwgts)) DEALLOCATE(tpwgts)
    IF(ALLOCATED(randomSeeds)) DEALLOCATE(randomSeeds)
    ERRORSEXITS("Decomposition_ElementDomainCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_ElementDomainCalculate
  
  !
  !================================================================================================================================
  !

  !>Gets the domain for a given element in a decomposition of a mesh. \todo should be able to specify lists of elements. \see OpenCMISS::Iron::cmfe_Decomposition_ElementDomainGet
  SUBROUTINE Decomposition_ElementDomainGet(decomposition,userElementNumber,domainNumber,err,error,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer to the decomposition to set the element domain for
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number to set the domain for.
    INTEGER(INTG), INTENT(OUT) :: domainNumber !<On return, the domain of the global element.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables`
    INTEGER(INTG) :: globalElementNumber
    TYPE(MeshType), POINTER :: mesh
    TYPE(MeshElementsType), POINTER :: meshElements
    TYPE(MeshTopologyType), POINTER :: meshTopology

    ENTERS("Decomposition_ElementDomainGet",err,error,*999)

    globalElementNumber=0
    CALL Decomposition_AssertIsFinished(decomposition,err,error,*999)       
    NULLIFY(mesh)
    CALL Decomposition_MeshGet(decomposition,mesh,err,error,*999)
    NULLIFY(meshElements)
    CALL MeshTopology_MeshElementsGet(meshTopology,meshElements,err,error,*999)
    CALL MeshElements_GlobalElementNumberGet(meshElements,userElementNumber,globalElementNumber,err,error,*999)
    domainNumber=decomposition%elementDomain(globalElementNumber)
    
    EXITS("Decomposition_ElementDomainGet")
    RETURN
999 ERRORSEXITS("Decomposition_ElementDomainGet",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_ElementDomainGet
  
  !
  !================================================================================================================================
  !

  !>Sets the domain for a given element in a decomposition of a mesh. \todo should be able to specify lists of elements. \see OpenCMISS::Iron::cmfe_Decomposition_ElementDomainSet 
  SUBROUTINE Decomposition_ElementDomainSet(decomposition,userElementNumber,domainNumber,err,error,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer to the decomposition to set the element domain for
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number to set the domain for.
    INTEGER(INTG), INTENT(IN) :: domainNumber !<The domain of the global element.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: globalElementNumber,numberOfComputationNodes
    TYPE(MeshType), POINTER :: mesh
    TYPE(MeshElementsType), POINTER :: meshElements
    TYPE(MeshTopologyType), POINTER :: meshTopology
    TYPE(VARYING_STRING) :: localError

    ENTERS("Decomposition_ElementDomainSet",err,error,*999)

    CALL Decomposition_AssertNotFinished(decomposition,err,error,*999)
    NULLIFY(mesh)
    CALL Decomposition_MeshGet(decomposition,mesh,err,error,*999)
    NULLIFY(meshTopology)
    CALL Mesh_MeshTopologyGet(mesh,0,meshTopology,err,error,*999)
    NULLIFY(meshElements)
    CALL MeshTopology_MeshElementsGet(meshTopology,meshElements,err,error,*999)
    CALL MeshElements_GlobalElementNumberGet(meshElements,userElementNumber,globalElementNumber,err,error,*999)
    CALL WorkGroup_NumberOfGroupNodesGet(decomposition%workGroup,numberOfComputationNodes,err,error,*999)
    IF(domainNumber<0.OR.domainNumber>=numberOfComputationNodes) THEN
      localError="Domain number "//TRIM(NumberToVString(domainNumber,"*",err,error))// &
        & " is invalid. The limits are 0 to "//TRIM(NumberToVString(numberOfComputationNodes,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    decomposition%elementDomain(globalElementNumber)=domainNumber
    
    EXITS("Decomposition_ElementDomainSet")
    RETURN
999 ERRORSEXITS("Decomposition_ElementDomainSet",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_ElementDomainSet
  
  !
  !================================================================================================================================
  !

  !>Finalises a decomposition and deallocates all memory.
  SUBROUTINE Decomposition_Finalise(decomposition,err,error,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer to the decomposition to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Decomposition_Finalise",err,error,*999)
    
    IF(ASSOCIATED(decomposition)) THEN
      !Destroy all the decomposition components          
      IF(ALLOCATED(decomposition%elementDomain)) DEALLOCATE(decomposition%elementDomain)
      CALL Decomposition_TopologyFinalise(decomposition%topology,err,error,*999)
      CALL Decomposition_DomainsFinalise(decomposition,err,error,*999)
      DEALLOCATE(decomposition)
    ENDIF
    
    EXITS("Decomposition_Finalise")
    RETURN
999 ERRORSEXITS("Decomposition_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_Finalise
  
  !
  !================================================================================================================================
  !

  !>Inintialises a decomposition and allocates all memory.
  SUBROUTINE Decomposition_Initialise(decomposition,err,error,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer to the decomposition to initialise. Must not be allocated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("Decomposition_Initialise",err,error,*998)

    IF(ASSOCIATED(decomposition)) CALL FlagError("Decomposition is already associated.",err,error,*998)
    
    ALLOCATE(decomposition,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate decomposition.",err,error,*999)
    decomposition%globalNumber=0
    decomposition%userNumber=0
    decomposition%decompositionFinished=.FALSE.
    NULLIFY(decomposition%decompositions)
    NULLIFY(decomposition%decomposer)
    NULLIFY(decomposition%mesh)
    NULLIFY(decomposition%region)
    NULLIFY(decomposition%interface)
    decomposition%numberOfDimensions=0
    decomposition%numberOfComponents=0
    decomposition%domainDecompositionType=DECOMPOSITION_ALL_TYPE
    NULLIFY(decomposition%workGroup)
    decomposition%numberOfDomains=0
    decomposition%numberOfEdgesCut=0
    decomposition%numberOfElements=0
    NULLIFY(decomposition%topology)
    decomposition%calculateLines=.TRUE. 
    decomposition%calculateFaces=.TRUE.
    
    EXITS("Decomposition_Initialise")
    RETURN
999 CALL Decomposition_Finalise(decomposition,dummyErr,dummyError,*998)
998 ERRORSEXITS("Decomposition_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_Initialise
  
  !
  !================================================================================================================================
  !
 
  !>Gets the mesh component number which will be used for the decomposition of a mesh. \see OpenCMISS::Iron::cmfe_Decomposition_MeshComponentGet
  SUBROUTINE Decomposition_MeshComponentNumberGet(decomposition,meshComponentNumber,err,error,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer to the decomposition to get the mesh component for
    INTEGER(INTG), INTENT(OUT) :: meshComponentNumber !<On return, the mesh component number to get
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Decomposition_MeshComponentNumberGet",err,error,*999)

    CALL Decomposition_AssertIsFinished(decomposition,err,error,*999)

    meshComponentNumber=decomposition%meshComponentNumber
  
    EXITS("Decomposition_MeshComponentNumberGet")
    RETURN
999 ERRORSEXITS("Decomposition_MeshComponentNumberGet",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_MeshComponentNumberGet
  
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the mesh component number which will be used for the decomposition of a mesh. \see OpenCMISS::Iron::cmfe_DecompositionMeshComponentSet
  SUBROUTINE Decomposition_MeshComponentNumberSet(decomposition,meshComponentNumber,err,error,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer to the decomposition to set the mesh component for
    INTEGER(INTG), INTENT(IN) :: meshComponentNumber !<The mesh component number to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Decomposition_MeshComponentNumberSet",err,error,*999)

    CALL Decomposition_AssertNotFinished(decomposition,err,error,*999)
    IF(meshComponentNumber<=0.OR.meshComponentNumber>decomposition%numberOfComponents) THEN
      localError="The specified mesh component number of "//TRIM(NumberToVString(meshComponentNumber,"*",err,error))// &
        & "is invalid. The component number must be between 1 and "// &
        & TRIM(NumberToVString(decomposition%numberOfComponents,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    decomposition%meshComponentNumber=meshComponentNumber
  
    EXITS("Decomposition_MeshComponentNumberSet")
    RETURN
999 ERRORSEXITS("Decomposition_MeshComponentNumberSet",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_MeshComponentNumberSet
  
  !
  !================================================================================================================================
  !

  !>Gets the domain for a given node in a decomposition of a mesh. \todo should be able to specify lists of elements. \see OpenCMISS::Iron::cmfe_Decomposition_NodeDomainGet
  SUBROUTINE Decomposition_NodeDomainGet(decomposition,userNodeNumber,meshComponentNumber,domainNumber,err,error,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer to the decomposition to set the element domain for
    INTEGER(INTG), INTENT(IN) :: userNodeNumber !<The global element number to set the domain for.
    INTEGER(INTG), INTENT(IN) :: meshComponentNumber !<The mesh component number to set the domain for.
    INTEGER(INTG), INTENT(OUT) :: domainNumber !<On return, the domain of the global element.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: meshNodeNumber
    TYPE(DomainType), POINTER :: domain
    TYPE(MeshType), POINTER :: mesh
    TYPE(MeshNodesType), POINTER :: meshNodes
    TYPE(MeshTopologyType), POINTER :: meshTopology
    TYPE(VARYING_STRING) :: localError

    ENTERS("Decomposition_NodeDomainGet",err,error,*999)

    meshNodeNumber=0
    CALL Decomposition_AssertIsFinished(decomposition,err,error,*999)
    NULLIFY(mesh)
    CALL Decomposition_MeshGet(decomposition,mesh,err,error,*999)
    IF(meshComponentNumber<=0.OR.meshComponentNumber>mesh%numberOfComponents) THEN
      localError="Mesh component number "//TRIM(NumberToVString(meshComponentNumber,"*",err,error))// &
        & " is invalid. The mesh component number should be >= 1 and <= "// &
        & TRIM(NumberToVString(mesh%numberOfComponents,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    NULLIFY(meshTopology)
    CALL Mesh_MeshTopologyGet(mesh,meshComponentNumber,meshTopology,err,error,*999)
    NULLIFY(meshNodes)
    CALL MeshTopology_MeshNodesGet(meshTopology,meshNodes,err,error,*999)
    CALL MeshNodes_MeshNodeNumberGet(meshNodes,userNodeNumber,meshNodeNumber,err,error,*999)
    NULLIFY(domain)
    CALL Decomposition_DomainGet(decomposition,meshComponentNumber,domain,err,error,*999)
    domainNumber=domain%nodeDomain(meshNodeNumber)
    
    EXITS("Decomposition_NodeDomainGet")
    RETURN
999 ERRORSEXITS("Decomposition_NodeDomainGet",err,error)
    RETURN 1

  END SUBROUTINE Decomposition_NodeDomainGet

  !
  !================================================================================================================================
  !

  !>Gets the number of domains for a decomposition. \see OpenCMISS::Iron::cmfe_DecompositionNumberOfDomainsGet
  SUBROUTINE DECOMPOSITION_NUMBER_OF_DOMAINS_GET(decomposition,numberOfDomains,err,error,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer to the decomposition to get the number of domains for.
    INTEGER(INTG), INTENT(OUT) :: numberOfDomains !<On return, the number of domains to get.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("decomposition_NUMBER_OF_DOMAINS_GET",err,error,*999)

    IF(ASSOCIATED(decomposition)) THEN
      IF(decomposition%decompositionFinished) THEN
        CALL FlagError("Decomposition has been finished.",err,error,*999)
      ELSE
        numberOfDomains=decomposition%numberOfDomains
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",err,error,*999)
    ENDIF

    EXITS("DECOMPOSITION_NUMBER_OF_DOMAINS_GET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_NUMBER_OF_DOMAINS_GET",err,error)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_NUMBER_OF_DOMAINS_GET


  !
  !================================================================================================================================
  !

  !>Sets/changes the number of domains for a decomposition. \see OpenCMISS::Iron::cmfe_DecompositionNumberOfDomainsSet
  SUBROUTINE DECOMPOSITION_NUMBER_OF_DOMAINS_SET(decomposition,numberOfDomains,err,error,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer to the decomposition to set the number of domains for.
    INTEGER(INTG), INTENT(IN) :: numberOfDomains !<The number of domains to set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfGroupComputationNodes
    TYPE(VARYING_STRING) :: localError

    ENTERS("DECOMPOSITION_NUMBER_OF_DOMAINS_SET",err,error,*999)

    IF(ASSOCIATED(decomposition)) THEN
      IF(decomposition%decompositionFinished) THEN
        CALL FlagError("Decomposition has been finished.",err,error,*999)
      ELSE
        SELECT CASE(decomposition%domainDecompositionType)
        CASE(DECOMPOSITION_ALL_TYPE)
          IF(numberOfDomains==1) THEN
            decomposition%numberOfDomains=1
          ELSE
            CALL FlagError("Can only have one domain for all decomposition type.",err,error,*999)
          ENDIF
        CASE(DECOMPOSITION_CALCULATED_TYPE,DECOMPOSITION_USER_DEFINED_TYPE)
          IF(numberOfDomains>=1) THEN
            !wolfye???<=?
            IF(numberOfDomains<=decomposition%numberOfElements) THEN
              !Get the number of computation nodes
              CALL WorkGroup_NumberOfGroupNodesGet(decomposition%workGroup,numberOfGroupComputationNodes, &
                & err,error,*999)
!!TODO: relax this later
              !IF(numberOfDomains==numberOfGroupComputationNodes) THEN
              decomposition%numberOfDomains=numberOfDomains             
              !ELSE
              !  localError="The number of domains ("//TRIM(NumberToVString(numberOfDomainsS,"*",err,error))// &
              !    & ") is not equal to the number of computation nodes ("// &
              !    & TRIM(NumberToVString(numberOfGroupComputationNodes,"*",err,error))//")"
              !  CALL FlagError(localError,err,error,*999)
              !ENDIF
            ELSE
              localError="The number of domains ("//TRIM(NumberToVString(numberOfDomains,"*",err,error))// &
                & ") must be <= the number of global elements ("// &
                & TRIM(NumberToVString(decomposition%numberOfElements,"*",err,error))//") in the mesh."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Number of domains must be >= 1.",err,error,*999)
          ENDIF
        CASE DEFAULT
          localError="Decomposition type "//TRIM(NumberToVString(decomposition%domainDecompositionType,"*",err,error))// &
            & " is not valid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FlagError("Decomposition is not associated.",err,error,*999)
    ENDIF

    EXITS("DECOMPOSITION_NUMBER_OF_DOMAINS_SET")
    RETURN
999 ERRORSEXITS("DECOMPOSITION_NUMBER_OF_DOMAINS_SET",err,error)
    RETURN 1
  END SUBROUTINE DECOMPOSITION_NUMBER_OF_DOMAINS_SET

  !
  !================================================================================================================================
  !

  !>Sets the workgroup to use for a decomposition on a given mesh. \see OpenCMISS::Iron::cmfe_Decomposition_WorkGroupSet
  SUBROUTINE Decomposition_WorkGroupSet(decomposition,workGroup,err,error,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer to the decomposition to set the work group for
    TYPE(WorkGroupType), POINTER :: workGroup !<A pointer to the work group to use for the decomposition
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Decomposition_WorkGroupSet",err,error,*999)

    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)
    IF(decomposition%decompositionFinished) CALL FlagError("Decomposition has already been finished.",err,error,*999)
    IF(.NOT.ASSOCIATED(workGroup)) CALL FlagError("Work group is not associated.",err,error,*999)
    IF(.NOT.workGroup%workGroupFinished) CALL FlagError("Work group has not been finished.",err,error,*999)

    decomposition%workGroup=>workGroup

    EXITS("Decomposition_WorkGroupSet")
    RETURN
999 ERRORSEXITS("Decomposition_WorkGroupSet",err,error)
    RETURN 1

  END SUBROUTINE Decomposition_WorkGroupSet

  !
  !================================================================================================================================
  !

  !>Calculates the topology for a decomposition.
  SUBROUTINE Decomposition_TopologyCalculate(decomposition,err,error,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer to the decomposition to finish creating.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: meshComponentNumber
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology

    ENTERS("Decomposition_TopologyCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)

    NULLIFY(decompositionTopology)
    CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*999)
    !Calculate the elements topology
    CALL DecompositionTopology_ElementsCalculate(decompositionTopology,err,error,*999)
    !Calculate the line topology
    IF(decomposition%calculateLines) CALL DecompositionTopology_LinesCalculate(decompositionTopology,err,error,*999)
    !Calculate the face topology
    IF(decomposition%calculateFaces) CALL DecompositionTopology_FacesCalculate(decompositionTopology,err,error,*999)
    !Calculate the data points topology
    meshComponentNumber=decomposition%meshComponentNumber
    IF(ALLOCATED(decomposition%mesh%topology(meshComponentNumber)%ptr%dataPoints%dataPoints)) &
      & CALL DecompositionTopology_DataPointsCalculate(decompositionTopology,err,error,*999)
 
    EXITS("Decomposition_TopologyCalculate")
    RETURN
999 ERRORSEXITS("Decomposition_TopologyCalculate",err,error)
    RETURN 1
  END SUBROUTINE Decomposition_TopologyCalculate

  !
  !================================================================================================================================
  !

  !>Finalises the element data points for a decomposition element.
  SUBROUTINE DecompositionElementDataPoints_Finalise(elementDataPoints,err,error,*)

    !Argument variables
    TYPE(DecompositionElementDataPointsType) :: elementDataPoints !<The element data points to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DecompositionElementDataPoints_Finalise",err,error,*999)

    elementDataPoints%numberOfProjectedData=0
    elementDataPoints%globalElementNumber=0
    IF(ALLOCATED(elementDataPoints%dataIndices)) DEALLOCATE(elementDataPoints%dataIndices)
    
    EXITS("DecompositionElementDataPoints_Finalise")
    RETURN
999 ERRORSEXITS("DecompositionElementDataPoints_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE DecompositionElementDataPoints_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the element data points for a decomposition element.
  SUBROUTINE DecompositionElementDataPoints_Initialise(elementDataPoints,err,error,*)

    !Argument variables
    TYPE(DecompositionElementDataPointsType) :: elementDataPoints !<The element data points to intialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DecompositionElementDataPoints_Initialise",err,error,*999)

    elementDataPoints%numberOfProjectedData=0
    elementDataPoints%globalElementNumber=0
    
    EXITS("DecompositionElementDataPoints_Initialise")
    RETURN
999 ERRORSEXITS("DecompositionElementDataPoints_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE DecompositionElementDataPoints_Initialise

  !
  !================================================================================================================================
  !

  !>Calculates the decomposition element topology.
  SUBROUTINE DecompositionTopology_DataPointsCalculate(decompositionTopology,err,error,*)

    !Argument variables
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology !<A pointer to the decomposition topology to calculate the elements for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: localElementIdx,globalElement,globalElementIdx,dataPointIdx,localData,meshComponentNumber, &
      & groupCommunicator
    INTEGER(INTG) :: insertStatus,mpiIError,numberOfComputationNodes,myGroupComputationNodeNumber,numberOfGhostData, &
      & numberOfLocalData
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(DecompositionDataPointsType), POINTER :: decompositionData
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainMappingType), POINTER :: elementsMapping
    TYPE(DomainMappingsType), POINTER :: domainMappings
    TYPE(MeshType), POINTER :: mesh
    TYPE(MeshTopologyType), POINTER :: meshTopology
    TYPE(MeshDataPointsType), POINTER :: meshData

    ENTERS("DecompositionTopology_DataPointsCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(decompositionTopology)) CALL FlagError("Decomposition topology is not associated.",err,error,*999)
    NULLIFY(decompositionData)
    CALL DecompositionTopology_DecompositionDataPointsGet(decompositionTopology,decompositionData,err,error,*999)
    NULLIFY(decomposition)
    CALL DecompositionTopology_DecompositionGet(decompositionTopology,decomposition,err,error,*999)
    CALL WorkGroup_GroupCommunicatorGet(decomposition%workGroup,groupCommunicator,err,error,*999)
    NULLIFY(decompositionElements)
    CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology,decompositionElements,err,error,*999)
    NULLIFY(domain)
    CALL Decomposition_DomainGet(decomposition,meshComponentNumber,domain,err,error,*999)
    NULLIFY(domainMappings)
    CALL Domain_DomainMappingsGet(domain,domainMappings,err,error,*999)
    NULLIFY(elementsMapping)
    CALL DomainMappings_ElementsMappingGet(domainMappings,elementsMapping,err,error,*999)
    meshComponentNumber=decomposition%meshComponentNumber
    NULLIFY(mesh)
    CALL Decomposition_MeshGet(decomposition,mesh,err,error,*999)
    NULLIFY(meshTopology)
    CALL Mesh_MeshTopologyGet(mesh,meshComponentNumber,meshTopology,err,error,*999)
    NULLIFY(meshData)
    CALL MeshTopology_MeshDataPointsGet(meshTopology,meshData,err,error,*999)
    CALL WorkGroup_NumberOfGroupNodesGet(decomposition%workGroup,numberOfComputationNodes,err,error,*999)
    CALL WorkGroup_GroupNodeNumberGet(decomposition%workGroup,myGroupComputationNodeNumber,err,error,*999)
    ALLOCATE(decompositionData%numberOfDomainLocal(0:numberOfComputationNodes-1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate decomposition data number of domain local.",err,error,*999)
    ALLOCATE(decompositionData%numberOfDomainGhost(0:numberOfComputationNodes-1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate decomposition data number of domain ghost.",err,error,*999)
    ALLOCATE(decompositionData%numberOfElementDataPoints(decompositionElements%numberOfGlobalElements),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate decomposition dat number of element data.",err,error,*999)
    ALLOCATE(decompositionData%elementDataPoints(decompositionElements%totalNumberOfElements),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate decomposition element data points.",err,error,*999)
    CALL Tree_CreateStart(decompositionData%dataPointsTree,err,error,*999)
    CALL Tree_InsertTypeSet(decompositionData%dataPointsTree,TREE_NO_DUPLICATES_ALLOWED,err,error,*999)
    CALL Tree_CreateFinish(decompositionData%dataPointsTree,err,error,*999)
    decompositionData%numberOfGlobalDataPoints=meshData%totalNumberOfProjectedData 
    DO globalElementIdx=1,decompositionElements%numberOfGlobalElements
      decompositionData%numberOfElementDataPoints(globalElementIdx)= &
        & meshData%elementDataPoints(globalElementIdx)%numberOfProjectedData
    ENDDO !globalElementIdx
    localData=0;
    DO localElementIdx=1,decompositionElements%totalNumberOfElements
      CALL DecompositionElementDataPoints_Initialise(decompositionData%elementDataPoints(localElementIdx),err,error,*999)
      globalElement=decompositionElements%elements(localElementIdx)%globalNumber
      decompositionData%elementDataPoints(localElementIdx)%numberOfProjectedData= &
        & meshData%elementDataPoints(globalElement)%numberOfProjectedData
      decompositionData%elementDataPoints(localElementIdx)%globalElementNumber=globalElement
      IF(localElementIdx<elementsMapping%ghostStart) THEN
        decompositionData%numberOfDataPoints=decompositionData%numberOfDataPoints+ &
          & decompositionData%elementDataPoints(localElementIdx)%numberOfProjectedData
      ENDIF
      decompositionData%totalNumberOfDataPoints=decompositionData%totalNumberOfDataPoints+ &
        & decompositionData%elementDataPoints(localElementIdx)%numberOfProjectedData
      ALLOCATE(decompositionData%elementDataPoints(localElementIdx)%dataIndices(decompositionData% &
        & elementDataPoints(localElementIdx)%numberOfProjectedData),STAT=err)
      DO dataPointIdx=1,decompositionData%elementDataPoints(localElementIdx)%numberOfProjectedData
        decompositionData%elementDataPoints(localElementIdx)%dataIndices(dataPointIdx)%userNumber= &
          & meshData%elementDataPoints(globalElement)%dataIndices(dataPointIdx)%userNumber
        decompositionData%elementDataPoints(localElementIdx)%dataIndices(dataPointIdx)%globalNumber= &
          & meshData%elementDataPoints(globalElement)%dataIndices(dataPointIdx)%globalNumber
        localData=localData+1
        decompositionData%elementDataPoints(localElementIdx)%dataIndices(dataPointIdx)%localNumber=localData
        CALL Tree_ItemInsert(decompositionData%dataPointsTree,decompositionData% &
          & elementDataPoints(localElementIdx)%dataIndices(dataPointIdx)%userNumber,localData, &
          & insertStatus,err,error,*999)
      ENDDO !dataPointIdx
    ENDDO !localElementIdx
    !Calculate number of ghost data points on the current computation domain
    numberOfLocalData=decompositionData%numberOfDataPoints
    numberOfGhostData=decompositionData%totalNumberOfDataPoints-decompositionData%numberOfDataPoints
    !Gather number of local data points on all computation nodes
    CALL MPI_ALLGATHER(numberOfLocalData,1,MPI_INTEGER,decompositionData% &
      & numberOfDomainLocal,1,MPI_INTEGER,groupCommunicator,mpiIError)
    CALL MPI_ErrorCheck("MPI_ALLGATHER",mpiIError,err,error,*999)
    !Gather number of ghost data points on all computation nodes
    CALL MPI_ALLGATHER(numberOfGhostData,1,MPI_INTEGER,decompositionData% &
      & numberOfDomainGhost,1,MPI_INTEGER,groupCommunicator,mpiIError)
    CALL MPI_ErrorCheck("MPI_ALLGATHER",mpiIError,err,error,*999)

    EXITS("DecompositionTopology_DataPointsCalculate")
    RETURN
999 ERRORSEXITS("DecompositionTopology_DataPointsCalculate",err,error)
    RETURN 1

  END SUBROUTINE DecompositionTopology_DataPointsCalculate

  !
  !================================================================================================================================
  !

  !>Calculates the decomposition element topology for a data projection (for data projections on fields).
  SUBROUTINE Decomposition_DataProjectionCalculate(decomposition,err,error,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer to the decomposition to calculate the elements for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology

    ENTERS("Decomposition_DataProjectionCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)
    NULLIFY(decompositionTopology)
    CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*999)
    
    CALL DecompositionTopology_DataPointsInitialise(decompositionTopology,err,error,*999)
    CALL DecompositionTopology_DataPointsCalculate(decompositionTopology,err,error,*999)

    EXITS("Decomposition_DataProjectionCalculate")
    RETURN
999 ERRORS("Decomposition_DataProjectionCalculate",err,error)
    EXITS("Decomposition_DataProjectionCalculate")
    RETURN 1

  END SUBROUTINE Decomposition_DataProjectionCalculate

  !
  !================================================================================================================================
  !

  !>Gets the local data point number for data points projected on an element
  SUBROUTINE Decomposition_ElementDataPointLocalNumberGet(decomposition,elementUserNumber,dataPointIndex, &
    & dataPointLocalNumber,err,error,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer to the decomposition to calculate the elements for
    INTEGER(INTG), INTENT(IN) :: elementUserNumber !<The element number to get the data point for
    INTEGER(INTG), INTENT(IN) :: dataPointIndex !<The data point index to get the number for
    INTEGER(INTG), INTENT(OUT) :: dataPointLocalNumber !<The data point local number to return
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: numberOfDataPoints,elementLocalNumber
    LOGICAL :: ghostElement
    TYPE(DecompositionDataPointsType), POINTER :: decompositionData
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(VARYING_STRING) :: localError

    ENTERS("Decomposition_ElementDataPointLocalNumberGet",err,error,*999)

    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)
    NULLIFY(decompositionTopology)
    CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*999)
    NULLIFY(decompositionElements)
    CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology,decompositionElements,err,error,*999)
    CALL DecompositionElements_LocalElementNumberGet(decompositionElements,elementUserNumber,elementLocalNumber, &
      & ghostElement,err,error,*999)
    NULLIFY(decompositionData)
    CALL DecompositionTopology_DecompositionDataPointsGet(decompositionTopology,decompositionData,err,error,*999)
    numberOfDataPoints = decompositionData%elementDataPoints(elementLocalNumber)%numberOfProjectedData
    IF(dataPointIndex<=0.OR.dataPointIndex>numberOfDataPoints) THEN
     localError="Element data point index "//TRIM(NumberToVString(dataPointIndex,"*",err,error))// &
        & " out of range for element "//TRIM(NumberToVString(elementUserNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    dataPointLocalNumber = decompositionData%elementDataPoints(elementLocalNumber)%dataIndices(dataPointIndex)%localNumber

    EXITS("Decomposition_ElementDataPointLocalNumberGet")
    RETURN
999 ERRORS("Decomposition_ElementDataPointLocalNumberGet",err,error)
    EXITS("Decomposition_ElementDataPointLocalNumberGet")
    RETURN 1
    
  END SUBROUTINE Decomposition_ElementDataPointLocalNumberGet

  !
  !================================================================================================================================
  !

  !>Gets the user data point number for data points projected on an element
  SUBROUTINE Decomposition_ElementDataPointUserNumberGet(decomposition,userElementNumber,dataPointIndex, &
    & dataPointUserNumber,err,error,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer to the decomposition topology to calculate the elements for
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The element number to get the data point for
    INTEGER(INTG), INTENT(IN) :: dataPointIndex !<The data point index to get the number for
    INTEGER(INTG), INTENT(OUT) :: dataPointUserNumber !<The data point user number to return
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: numberOfDataPoints,localElementNumber
    LOGICAL :: ghostElement
    TYPE(DecompositionDataPointsType), POINTER :: decompositionData
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(VARYING_STRING) :: localError

    ENTERS("Decomposition_ElementDataPointUserNumberGet",err,error,*999)

    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)

    NULLIFY(decompositionTopology)
    CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*999)
    NULLIFY(decompositionElements)
    CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology,decompositionElements,err,error,*999)
    NULLIFY(decompositionData)
    CALL DecompositionTopology_DecompositionDataPointsGet(decompositionTopology,decompositionData,err,error,*999)
    CALL DecompositionElements_LocalElementNumberGet(decompositionElements,userElementNumber,localElementNumber, &
       & ghostElement,err,error,*999)
    numberOfDataPoints = decompositionData%elementDataPoints(localElementNumber)%numberOfProjectedData
    IF(dataPointIndex<=0.OR.dataPointIndex>numberOfDataPoints) THEN
      localError="Element data point index "//TRIM(NumberToVString(dataPointIndex,"*",err,error))// &
        & " out of range for element "//TRIM(NumberToVString(userElementNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    dataPointUserNumber=decompositionData%elementDataPoints(localElementNumber)%dataIndices(dataPointIndex)%userNumber

    EXITS("Decomposition_ElementDataPointUserNumberGet")
    RETURN
999 ERRORS("Decomposition_ElementDataPointUserNumberGet",err,error)
    EXITS("Decomposition_ElementDataPointUserNumberGet")
    RETURN 1

  END SUBROUTINE Decomposition_ElementDataPointUserNumberGet

  !
  !================================================================================================================================
  !

  !>Gets the number of data points projected on an element
  SUBROUTINE Decomposition_NumberOfElementDataPointsGet(decomposition,userElementNumber,numberOfDataPoints,err,error,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer to the decomposition to calculate the elements for
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The element number to get the data point for
    INTEGER(INTG), INTENT(OUT) :: numberOfDataPoints !<The data point local number to return
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: localElementNumber
    LOGICAL :: ghostElement
    TYPE(DecompositionDataPointsType), POINTER :: decompositionData
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology

    ENTERS("Decomposition_NumberOfElementDataPointsGet",err,error,*999)

    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)
    
    NULLIFY(decompositionTopology)
    CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*999)
    NULLIFY(decompositionElements)
    CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology,decompositionElements,err,error,*999)
    NULLIFY(decompositionData)
    CALL DecompositionTopology_DecompositionDataPointsGet(decompositionTopology,decompositionData,err,error,*999)
    CALL DecompositionElements_LocalElementNumberGet(decompositionElements,userElementNumber,localElementNumber, &
      & ghostElement,err,error,*999)
    
    numberOfDataPoints=decompositionData%elementDataPoints(localElementNumber)%numberOfProjectedData

    EXITS("Decomposition_NumberOfElementDataPointsGet")
    RETURN
999 ERRORS("Decomposition_NumberOfElementDataPointsGet",err,error)
    EXITS("Decomposition_NumberOfElementDataPointsGet")
    RETURN 1
    
  END SUBROUTINE Decomposition_NumberOfElementDataPointsGet

  !
  !================================================================================================================================
  !

  !>Finalises the given decomposition element.
  SUBROUTINE DecompositionElement_Finalise(decompositionElement,err,error,*)

    !Argument variables
    TYPE(DecompositionElementType) :: decompositionElement !<The decomposition element to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: xiCoordIdx

    ENTERS("DecompositionElement_Finalise",err,error,*999)

    IF(ALLOCATED(decompositionElement%adjacentElements)) THEN
      DO xiCoordIdx=LBOUND(decompositionElement%adjacentElements,1),UBOUND(decompositionElement%adjacentElements,1)
        CALL Decomposition_AdjacentElementFinalise(decompositionElement%adjacentElements(xiCoordIdx),err,error,*999)
      ENDDO !xiCoordIdx
      DEALLOCATE(decompositionElement%adjacentElements)
    ENDIF
    IF(ALLOCATED(decompositionElement%elementLines)) DEALLOCATE(decompositionElement%elementLines)
    IF(ALLOCATED(decompositionElement%elementFaces)) DEALLOCATE(decompositionElement%elementFaces)

    EXITS("DecompositionElement_Finalise")
    RETURN
999 ERRORSEXITS("DecompositionElement_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE DecompositionElement_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the given decomposition element.
  SUBROUTINE DecompositionElement_Initialise(decompositionElement,err,error,*)

    !Argument variables
    TYPE(DecompositionElementType) :: decompositionElement !<The decomposition element to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DecompositionElement_Initialise",err,error,*999)

    decompositionElement%userNumber=0
    decompositionElement%localNumber=0
    decompositionElement%globalNumber=0
    decompositionElement%boundaryElement=.FALSE.

    EXITS("DecompositionElement_Initialise")
    RETURN
999 ERRORSEXITS("DecompositionElement_Initialise",err,error)
    RETURN 1   
  END SUBROUTINE DecompositionElement_Initialise

  !
  !================================================================================================================================
  !

  !>Calculates the element numbers adjacent to an element in a decomposition topology.
  SUBROUTINE DecompositionTopology_ElementsAdjacentElementsCalculate(decompositionTopology,err,error,*)

    !Argument variables
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology !<A pointer to the decomposition topology to calculate the adjacent elements for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: j,elementIdx,elementIdx1,surroundingElementIdx1,xiIdx,xiCoordIdx, &
      & localNodeIdx,localNodeIdx1,localNodeIdx2,localNodeIdx3,nodeNumber,nodeNumber1,dummyErr, &
      & faceXi(2),faceXiC(3),nodePositionIndex(4)
    INTEGER(INTG) :: xiDirection,directionIndex,xiDirCheck,xiDirSearch,numberOfNodeMatches
    INTEGER(INTG) :: candidateIdx,faceNodeIdx,nodeIdx,surroundingElementIdx,candidateElement,idx
    INTEGER(INTG) :: numberSurrounding,numberOfNodesXiC(4),numberSurroundingElements
    INTEGER(INTG), ALLOCATABLE :: nodeMatches(:),adjacentElements(:),surroundingElements(:)
    LOGICAL :: xiCollapsed,faceCollapsed(-3:3),subset
    TYPE(LIST_TYPE), POINTER :: nodeMatchList, surroundingElementsList
    TYPE(LIST_PTR_TYPE) :: adjacentElementsList(-4:4)
    TYPE(BasisType), POINTER :: basis
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(DomainType), POINTER :: DOMAIN
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(VARYING_STRING) :: dummyError,localError

    NULLIFY(nodeMatchList)
    DO xiCoordIdx=-4,4
      NULLIFY(adjacentElementsList(xiCoordIdx)%ptr)
    ENDDO !xiCoordIdx

    ENTERS("DecompositionTopology_ElementsAdjacentElementsCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(decompositionTopology)) CALL FlagError("Decomposition topology is not associated.",err,error,*999)
    
    NULLIFY(decomposition)
    CALL DecompositionTopology_DecompositionGet(decompositionTopology,decomposition,err,error,*999)
    NULLIFY(decompositionElements)
    CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology,decompositionElements,err,error,*999)
    NULLIFY(domain)
    CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*999)
    NULLIFY(domainTopology)
    CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
    NULLIFY(domainNodes)
    CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
    NULLIFY(domainElements)
    CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
    
    !Loop over the elements in the decomposition
    DO elementIdx=1,decompositionElements%totalNumberOfElements
      NULLIFY(basis)
      CALL DomainElements_BasisGet(domainElements,elementIdx,basis,err,error,*999)
     !Create a list for every xi direction (plus and minus)
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
      CALL List_ItemAdd(adjacentElementsList(0)%ptr,decompositionElements%elements(elementIdx)%localNumber,err,error,*999)
      SELECT CASE(basis%TYPE)
      CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
!!TODO: Calculate this and set it as part of the basis type
        !Determine the collapsed "faces" if any
        nodePositionIndex=1
        !Loop over the face normals of the element
        DO xiIdx=1,basis%numberOfXi
          !Determine the face xi directions that lie in this xi direction
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
                    IF(domainElements%elements(elementIdx)%elementNodes(localNodeIdx1)/= &
                      & domainElements%elements(elementIdx)%elementNodes(localNodeIdx2)) xiCollapsed=.FALSE.
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
            !If the face is collapsed then don't look in this xi direction. The exception is if the opposite face is
            !also collapsed. This may indicate that we have a funny element in non-rc coordinates that goes around the
            !central axis back to itself
            IF(faceCollapsed(xiDirection).AND..NOT.faceCollapsed(-xiDirection)) THEN
              !Do nothing - the match lists are already empty
            ELSE
              !Find the nodes to match and add them to the node match list
              SELECT CASE(basis%numberOfXi)
              CASE(1)
                localNodeIdx=basis%nodePositionIndexInv(nodePositionIndex(1),1,1,1)
                IF(localNodeIdx/=0) THEN
                  nodeNumber=domainElements%elements(elementIdx)%elementNodes(localNodeIdx)
                  CALL List_ItemAdd(nodeMatchList,nodeNumber,err,error,*999)
                ENDIF
              CASE(2)
                DO localNodeIdx1=1,numberOfNodesXiC(faceXi(1)),numberOfNodesXiC(faceXi(1))-1
                  nodePositionIndex(faceXi(1))=localNodeIdx1
                  localNodeIdx=basis%nodePositionIndexInv(nodePositionIndex(1),nodePositionIndex(2),1,1)
                  IF(localNodeIdx/=0) THEN
                    nodeNumber=domainElements%elements(elementIdx)%elementNodes(localNodeIdx)
                    CALL List_ItemAdd(nodeMatchList,nodeNumber,err,error,*999)
                  ENDIF
                ENDDO !localNodeIdx1
              CASE(3)
                DO localNodeIdx1=1,numberOfNodesXiC(faceXi(1)),numberOfNodesXiC(faceXi(1))-1
                  nodePositionIndex(faceXi(1))=localNodeIdx1
                  DO localNodeIdx2=1,numberOfNodesXiC(faceXi(2)),numberOfNodesXiC(faceXi(2))-1
                    nodePositionIndex(faceXi(2))=localNodeIdx2
                    localNodeIdx=basis%nodePositionIndexInv(nodePositionIndex(1),nodePositionIndex(2),nodePositionIndex(3),1)
                    IF(localNodeIdx/=0) THEN
                      nodeNumber=domainElements%elements(elementIdx)%elementNodes(localNodeIdx)
                      CALL List_ItemAdd(nodeMatchList,nodeNumber,err,error,*999)
                    ENDIF
                  ENDDO !localNodeIdx2
                ENDDO !localNodeIdx1
              CASE DEFAULT
                localError="The number of xi directions in the basis of "// &
                  & TRIM(NumberToVString(basis%numberOfXi,"*",err,error))//" is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDIF
            CALL List_RemoveDuplicates(nodeMatchList,err,error,*999)
            CALL List_DetachAndDestroy(nodeMatchList,numberOfNodeMatches,nodeMatches,err,error,*999)
            numberSurrounding=0
            IF(numberOfNodeMatches>0) THEN
              !nodeMatches now contain the list of corner nodes in the current face with normal_xi of xiIdx.
              !Look at the surrounding elements of each of these nodes, if there is a repeated element that
              !is not the current element elementIdx, it's an adjacent element.
              candidateIdx=0
              NULLIFY(surroundingElementsList)
              CALL List_CreateStart(surroundingElementsList,err,error,*999)
              CALL List_DataTypeSet(surroundingElementsList,LIST_INTG_TYPE,err,error,*999)
              CALL List_InitialSizeSet(surroundingElementsList,2,err,error,*999)
              CALL List_CreateFinish(surroundingElementsList,err,error,*999)
              DO faceNodeIdx=1,numberOfNodeMatches
                !Dump all the surrounding elements into an array, see if any are repeated
                nodeIdx=nodeMatches(faceNodeIdx)
                DO surroundingElementIdx=1,domainNodes%nodes(nodeIdx)%numberOfSurroundingElements
                  candidateElement=domainNodes%nodes(nodeIdx)%surroundingElements(surroundingElementIdx)
                  IF(candidateElement/=elementIdx) THEN
                    candidateIdx=candidateIdx+1
                    CALL List_ItemAdd(surroundingElementsList,candidateElement,err,error,*999)
                  ENDIF
                ENDDO
              ENDDO !faceNodeIdx
              CALL List_DetachAndDestroy(surroundingElementsList,numberSurroundingElements,surroundingElements, &
                & err,error,*999)
              DO idx=1,candidateIdx
                elementIdx1=surroundingElements(idx)
                IF(COUNT(surroundingElements(1:numberSurroundingElements)==elementIdx1)>=basis%numberOfXi) THEN
                  !Found it, just exit
                  CALL List_ItemAdd(adjacentElementsList(xiDirection)%ptr,elementIdx1,err,error,*999)
                  numberSurrounding=numberSurrounding+1
                  EXIT
                ENDIF
              ENDDO
            ENDIF
            IF(ALLOCATED(nodeMatches)) DEALLOCATE(nodeMatches)
            IF(ALLOCATED(surroundingElements)) DEALLOCATE(surroundingElements)
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
                  & nodePositionIndex(4))
                IF(localNodeIdx/=0) THEN
                  nodeNumber=domainElements%elements(elementIdx)%elementNodes(localNodeIdx)
                  CALL List_ItemAdd(nodeMatchList,nodeNumber,err,error,*999)
                ENDIF
              ENDDO !localNodeIdx3
            ENDDO !localNodeIdx2
          ENDDO !localNodeIdx1
          CALL List_RemoveDuplicates(nodeMatchList,err,error,*999)
          CALL List_DetachAndDestroy(nodeMatchList,numberOfNodeMatches,nodeMatches,err,error,*999)
          IF(numberOfNodeMatches>0) THEN
            !Find list of elements surrounding those nodes
            DO nodeIdx=1,numberOfNodeMatches
              nodeNumber1=nodeMatches(nodeIdx)
              DO surroundingElementIdx1=1,domainNodes%nodes(nodeNumber1)%numberOfSurroundingElements
                elementIdx1=domainNodes%nodes(nodeNumber1)%surroundingElements(surroundingElementIdx1)
                IF(elementIdx1/=elementIdx) THEN !Don't want the current element
                  ! grab the nodes list for current and this surrouding elements
                  ! current face : nodeMatches
                  ! candidate elem : decompositionTopology%ELEMENTS%elements(elementIdx1)%meshElementNodes 
                  ! if all of current face belongs to the candidate element, we will have found the neighbour
                  CALL List_SubsetOf(nodeMatches(1:numberOfNodeMatches),domainElements%elements(elementIdx1)% &
                    & elementNodes,subset,err,error,*999)
                  IF(subset) THEN
                    CALL List_ItemAdd(adjacentElementsList(xiCoordIdx)%ptr,elementIdx1,err,error,*999)
                  ENDIF
                ENDIF
              ENDDO !surroundingElementIdx1
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
        localError="The basis type of "//TRIM(NumberToVString(basis%TYPE,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)                     
      END SELECT
      !Set the surrounding elements for this element
      ALLOCATE(decompositionElements%elements(elementIdx)%adjacentElements(-basis%numberOfXiCoordinates: &
        basis%numberOfXiCoordinates),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate adjacent elements.",err,error,*999)
      DO xiCoordIdx=-basis%numberOfXiCoordinates,basis%numberOfXiCoordinates
        CALL Decomposition_AdjacentElementInitialise(decompositionElements%elements(elementIdx)%adjacentElements(xiCoordIdx), &
          & err,error,*999)
        CALL List_DetachAndDestroy(adjacentElementsList(xiCoordIdx)%ptr,decompositionElements%elements(elementIdx)% &
          & adjacentElements(xiCoordIdx)%numberOfAdjacentElements,adjacentElements,err,error,*999)
        ALLOCATE(decompositionElements%elements(elementIdx)%adjacentElements(xiCoordIdx)%adjacentElements( &
          decompositionElements%elements(elementIdx)%adjacentElements(xiCoordIdx)%numberOfAdjacentElements),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate element adjacent elements.",err,error,*999)
        decompositionElements%elements(elementIdx)%adjacentElements(xiCoordIdx)%adjacentElements(1: &
          & decompositionElements%elements(elementIdx)%adjacentElements(xiCoordIdx)%numberOfAdjacentElements)= &
          adjacentElements(1:decompositionElements%elements(elementIdx)%adjacentElements(xiCoordIdx)%numberOfAdjacentElements)
        IF(ALLOCATED(adjacentElements)) DEALLOCATE(adjacentElements)
      ENDDO !xiCoordIdx
    ENDDO !elementIdx           

    IF(DIAGNOSTICS1) THEN
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Total number of elements = ",decompositionElements%totalNumberOfElements, &
        & err,error,*999)
      DO elementIdx=1,decompositionElements%totalNumberOfElements
        NULLIFY(basis)
        CALL DomainElements_BasisGet(domainElements,elementIdx,basis,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Local element number : ",elementIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of xi coordinates = ",basis%numberOfXiCoordinates,err,error,*999)
        DO xiCoordIdx=-basis%numberOfXiCoordinates,basis%numberOfXiCoordinates
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Xi coordinate : ",xiCoordIdx,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Number of adjacent elements = ", &
            & decompositionElements%elements(elementIdx)%adjacentElements(xiCoordIdx)%numberOfAdjacentElements,err,error,*999)
          IF(decompositionElements%elements(elementIdx)%adjacentElements(xiCoordIdx)%numberOfAdjacentElements>0) THEN
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,decompositionElements%elements(elementIdx)% &
              & adjacentElements(xiCoordIdx)%numberOfAdjacentElements,8,8,decompositionElements%elements(elementIdx)% &
              & adjacentElements(xiCoordIdx)%adjacentElements, &
              & '("        Adjacent elements :",8(X,I6))','(30x,8(X,I6))',err,error,*999)
          ENDIF
        ENDDO !xiCoordIdx
      ENDDO !elementIdx
    ENDIF

    EXITS("DecompositionTopology_ElementsAdjacentElementsCalculate")
    RETURN
999 IF(ALLOCATED(nodeMatches)) DEALLOCATE(nodeMatches)
    IF(ALLOCATED(surroundingElements)) DEALLOCATE(surroundingElements)
    IF(ALLOCATED(adjacentElements)) DEALLOCATE(adjacentElements)
    IF(ASSOCIATED(nodeMatchList)) CALL List_Destroy(nodeMatchList,dummyErr,dummyError,*998)
998 IF(ASSOCIATED(surroundingElementsList)) CALL List_Destroy(surroundingElementsList,dummyErr,dummyError,*997)
997 DO xiCoordIdx=-4,4
      IF(ASSOCIATED(adjacentElementsList(xiCoordIdx)%ptr)) &
        & CALL List_Destroy(adjacentElementsList(xiCoordIdx)%ptr,dummyErr,dummyError,*996)
    ENDDO !ni
996 ERRORS("DecompositionTopology_ElementsAdjacentElementsCalculate",err,error)
    EXITS("DecompositionTopology_ElementsAdjacentElementsCalculate")
    RETURN 1

  END SUBROUTINE DecompositionTopology_ElementsAdjacentElementsCalculate

  !
  !================================================================================================================================
  !

  !>Calculates the decomposition element topology.
  SUBROUTINE DecompositionTopology_ElementsCalculate(decompositionTopology,err,error,*)

    !Argument variables
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology !<A pointer to the decomposition topology to calculate the elements for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: globalElement,insertStatus,localElementIdx
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainMappingType), POINTER :: domainElementsMapping
    TYPE(DomainMappingsType), POINTER :: domainMappings
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(MeshType), POINTER :: mesh
    TYPE(MeshElementsType), POINTER :: meshElements
    TYPE(MeshTopologyType), POINTER :: meshTopology

    ENTERS("DecompositionTopology_ElementsCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(decompositionTopology)) CALL FlagError("Decomposition topology is not associated.",err,error,*999)
    NULLIFY(decompositionElements)
    CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology,decompositionElements,err,error,*999)
    NULLIFY(decomposition)
    CALL DecompositionTopology_DecompositionGet(decompositionTopology,decomposition,err,error,*999)
    NULLIFY(domain)
    CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*999)
    NULLIFY(domainTopology)
    CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
    NULLIFY(domainMappings)
    CALL Domain_DomainMappingsGet(domain,domainMappings,err,error,*999)
    NULLIFY(domainElements)
    CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
    NULLIFY(domainElementsMapping)
    CALL DomainMappings_ElementsMappingGet(domainMappings,domainElementsMapping,err,error,*999)
    NULLIFY(mesh)
    CALL Decomposition_MeshGet(decomposition,mesh,err,error,*999)
    NULLIFY(meshTopology)
    CALL Mesh_MeshTopologyGet(mesh,decomposition%meshComponentNumber,meshTopology,err,error,*999)
    NULLIFY(meshElements)
    CALL MeshTopology_MeshElementsGet(meshTopology,meshElements,err,error,*999)
    
    !Allocate the element topology arrays
    ALLOCATE(decompositionElements%elements(domainElements%totalNumberOfElements),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate decomposition elements elements.",err,error,*999)
    decompositionElements%numberOfElements=domainElements%numberOfElements
    decompositionElements%totalNumberOfElements=domainElements%totalNumberOfElements
    decompositionElements%numberOfGlobalElements=domainElements%numberOfGlobalElements
    CALL Tree_CreateStart(decompositionElements%elementsTree,err,error,*999)
    CALL Tree_InsertTypeSet(decompositionElements%elementsTree,TREE_NO_DUPLICATES_ALLOWED,err,error,*999)
    CALL Tree_CreateFinish(decompositionElements%elementsTree,err,error,*999)
    DO localElementIdx=1,decompositionElements%totalNumberOfElements
      CALL DecompositionElement_Initialise(decompositionElements%elements(localElementIdx),err,error,*999)
      globalElement=domainElementsMapping%localToGlobalMap(localElementIdx)
      decompositionElements%elements(localElementIdx)%userNumber=meshElements%elements(globalElement)%userNumber
      decompositionElements%elements(localElementIdx)%localNumber=localElementIdx
      CALL Tree_ItemInsert(decompositionElements%elementsTree,decompositionElements% &
        & elements(localElementIdx)%userNumber,localElementIdx,insertStatus,err,error,*999)
      decompositionElements%elements(localElementIdx)%globalNumber=globalElement
      decompositionElements%elements(localElementIdx)%boundaryElement=meshElements%elements(globalElement)%boundaryElement
    ENDDO !localElementIdx
    !Calculate the elements surrounding the elements in the decomposition topology
    CALL DecompositionTopology_ElementsAdjacentElementsCalculate(decompositionTopology,err,error,*999)

    EXITS("DecompositionTopology_ElementsCalculate")
    RETURN
999 ERRORSEXITS("DecompositionTopology_ElementsCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE DecompositionTopology_ElementsCalculate

  !
  !================================================================================================================================
  !

  !>Finalises the elements in the given decomposition topology. 
  SUBROUTINE DecompositionTopology_ElementsFinalise(decompositionElements,err,error,*)

    !Argument variables
    TYPE(DecompositionElementsType), POINTER :: decompositionElements !<A pointer to the decomposition elements to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementIdx

    ENTERS("DecompositionTopology_ElementsFinalise",err,error,*999)

    IF(ASSOCIATED(decompositionElements)) THEN
      DO elementIdx=1,SIZE(decompositionElements%elements,1)
        CALL DecompositionElement_Finalise(decompositionElements%elements(elementIdx),err,error,*999)
      ENDDO !elementIdx
      IF(ALLOCATED(decompositionElements%elements)) DEALLOCATE(decompositionElements%elements)
      IF(ASSOCIATED(decompositionElements%elementsTree)) CALL Tree_Destroy(decompositionElements%elementsTree,err,error,*999)
      DEALLOCATE(decompositionElements)
    ENDIF

    EXITS("DecompositionTopology_ElementsFinalise")
    RETURN
999 ERRORSEXITS("DecompositionTopology_ElementsFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE DecompositionTopology_ElementsFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the element data structures for a decomposition topology.
  SUBROUTINE DecompositionTopology_ElementsInitialise(decompositionTopology,err,error,*)

    !Argument variables
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology !<A pointer to the decomposition topology to initialise the elements for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("DecompositionTopology_ElementsInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(decompositionTopology)) CALL FlagError("Decomposition topology is not associated.",err,error,*998)
    IF(ASSOCIATED(decompositionTopology%elements)) &
      & CALL FlagError("Decomposition already has topology elements associated.",err,error,*998)
    NULLIFY(decomposition)
    CALL DecompositionTopology_DecompositionGet(decompositionTopology,decomposition,err,error,*998)
    
    ALLOCATE(decompositionTopology%elements,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate topology elements.",err,error,*999)
    decompositionTopology%elements%numberOfElements=0
    decompositionTopology%elements%totalNumberOfElements=0
    decompositionTopology%elements%numberOfGlobalElements=0
    decompositionTopology%elements%decompositionTopology=>decompositionTopology
    NULLIFY(decompositionTopology%elements%elementsTree)

    EXITS("DecompositionTopology_ElementsInitialise")
    RETURN
999 CALL DecompositionTopology_ElementsFinalise(decompositionTopology%elements,dummyErr,dummyError,*998)
998 ERRORSEXITS("DecompositionTopology_ElementsInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE DecompositionTopology_ElementsInitialise

  !
  !================================================================================================================================
  !

  !>Finalises the topology in the given decomposition.
  SUBROUTINE Decomposition_TopologyFinalise(decompositionTopology,err,error,*)

    !Argument variables
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology !<A pointer to the decomposition topology to finalise 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Decomposition_TopologyFinalise",err,error,*999)

    IF(ASSOCIATED(decompositionTopology)) THEN
      CALL DecompositionTopology_ElementsFinalise(decompositionTopology%elements,err,error,*999)
      CALL DecompositionTopology_LinesFinalise(decompositionTopology%lines,err,error,*999)
      CALL DecompositionTopology_FacesFinalise(decompositionTopology%faces,err,error,*999)
      DEALLOCATE(decompositionTopology)
    ENDIF

    EXITS("Decomposition_TopologyFinalise")
    RETURN
999 ERRORSEXITS("Decomposition_TopologyFinalise",err,error)
    RETURN 1

  END SUBROUTINE Decomposition_TopologyFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the topology for a given decomposition.
  SUBROUTINE Decomposition_TopologyInitialise(decomposition,err,error,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer to the decomposition to initialise the topology for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,meshComponentNumber
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("Decomposition_TopologyInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*998)
    IF(ASSOCIATED(decomposition%topology)) CALL FlagError("Decomposition already has topology associated.",err,error,*998)
    
    !Allocate decomposition topology
    ALLOCATE(decomposition%topology,STAT=err)
    IF(err/=0) CALL FlagError("Decomposition topology could not be allocated.",err,error,*999)
    decomposition%topology%decomposition=>decomposition
    NULLIFY(decomposition%topology%elements)
    NULLIFY(decomposition%topology%lines)
    NULLIFY(decomposition%topology%faces)
    NULLIFY(decomposition%topology%dataPoints)
    !Initialise the topology components
    CALL DecompositionTopology_ElementsInitialise(decomposition%TOPOLOGY,err,error,*999)
    IF(decomposition%calculateLines) CALL DecompositionTopology_LinesInitialise(decomposition%topology,err,error,*999)
    IF(decomposition%calculateFaces) CALL DecompositionTopology_FacesInitialise(decomposition%topology,err,error,*999)
    meshComponentNumber=decomposition%meshComponentNumber
    IF(ALLOCATED(decomposition%mesh%topology(meshComponentNumber)%ptr%dataPoints%dataPoints)) THEN
      CALL DecompositionTopology_DataPointsInitialise(decomposition%topology,err,error,*999)
    ENDIF

    EXITS("Decomposition_TopologyInitialise")
    RETURN
999 CALL Decomposition_TopologyFinalise(decomposition%topology,dummyErr,dummyError,*998)
998 ERRORSEXITS("Decomposition_TopologyInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_TopologyInitialise

  !
  !================================================================================================================================
  !

  !>Finalise the data point data structures for a decomposition topology.
  SUBROUTINE DecompositionTopology_DataPointsFinalise(decompositionDataPoints,err,error,*)

    !Argument variables
    TYPE(DecompositionDataPointsType), POINTER :: decompositionDataPoints !<A pointer to the decomposition data points to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementIdx

    ENTERS("DecompositionTopology_DataPointsFinalise",err,error,*999)

    IF(ASSOCIATED(decompositionDataPoints)) THEN
      IF(ALLOCATED(decompositionDataPoints%numberOfDomainLocal)) DEALLOCATE(decompositionDataPoints%numberOfDomainLocal)
      IF(ALLOCATED(decompositionDataPoints%numberOfDomainGhost)) DEALLOCATE(decompositionDataPoints%numberOfDomainGhost)
      IF(ALLOCATED(decompositionDataPoints%numberOfElementDataPoints)) &
        & DEALLOCATE(decompositionDataPoints%numberOfElementDataPoints)
      DO elementIdx=1,SIZE(decompositionDataPoints%elementDataPoints,1)
        CALL DecompositionElementDataPoints_Finalise(decompositionDataPoints%elementDataPoints(elementIdx),err,error,*999)
      ENDDO !elementIdx
      IF(ALLOCATED(decompositionDataPoints%elementDataPoints)) DEALLOCATE(decompositionDataPoints%elementDataPoints)
      IF(ASSOCIATED(decompositionDataPoints%dataPointsTree)) &
        & CALL Tree_Destroy(decompositionDataPoints%dataPointsTree,err,error,*999)
      DEALLOCATE(decompositionDataPoints)
    ENDIF

    EXITS("DecompositionTopology_DataPointsFinalise")
    RETURN
999 ERRORSEXITS("DecompositionTopology_DataPointsFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE DecompositionTopology_DataPointsFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the data point data structures for a decomposition topology.
  SUBROUTINE DecompositionTopology_DataPointsInitialise(decompositionTopology,err,error,*)

    !Argument variables
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology !<A pointer to the decomposition topology to initialise the lines for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("DecompositionTopology_DataPointsInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(decompositionTopology)) CALL FlagError("Decomposition topology is not associated.",err,error,*998)
    IF(ASSOCIATED(decompositionTopology%dataPoints)) &
      & CALL FlagError("Decomposition already has topology data points associated.",err,error,*998)
     
    ALLOCATE(decompositionTopology%dataPoints,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate topology data points.",err,error,*999)
    decompositionTopology%dataPoints%numberOfDataPoints=0
    decompositionTopology%dataPoints%totalNumberOfDataPoints=0
    decompositionTopology%dataPoints%numberOfGlobalDataPoints=0
    NULLIFY(decompositionTopology%dataPoints%dataPointsTree)
    decompositionTopology%dataPoints%decompositionTopology=>decompositionTopology

    EXITS("DecompositionTopology_DataPointsInitialise")
    RETURN
999 CALL DecompositionTopology_DataPointsFinalise(decompositionTopology%dataPoints,dummyErr,dummyError,*998)
998 ERRORSEXITS("DecompositionTopology_DataPointsInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE DecompositionTopology_DataPointsInitialise

  !
  !================================================================================================================================
  !

  !>Finalises a line in the given decomposition topology and deallocates all memory.
  SUBROUTINE DecompositionLine_Finalise(decompositionLine,err,error,*)

    !Argument variables
    TYPE(DecompositionLineType) :: decompositionLine !<The decomposition line to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DecompositionLine_Finalise",err,error,*999)

    decompositionLine%number=0
    decompositionLine%xiDirection=0
    decompositionLine%numberOfSurroundingElements=0
    IF(ALLOCATED(decompositionLine%surroundingElements)) DEALLOCATE(decompositionLine%surroundingElements)
    IF(ALLOCATED(decompositionLine%elementLines)) DEALLOCATE(decompositionLine%elementLines)
    decompositionLine%adjacentLines=0

    EXITS("DecompositionLine_Finalise")
    RETURN
999 ERRORSEXITS("DecompositionLine_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE DecompositionLine_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the line data structure for a decomposition topology line.
  SUBROUTINE DecompositionLine_Initialise(decompositionLine,err,error,*)

    !Argument variables
    TYPE(DecompositionLineType) :: decompositionLine !<The decomposition line to initialise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DecompositionLine_Initialise",err,error,*999)

    decompositionLine%number=0
    decompositionLine%xiDirection=0
    decompositionLine%numberOfSurroundingElements=0
    decompositionLine%adjacentLines=0
    decompositionLine%boundaryLine=.FALSE.

    EXITS("DecompositionLine_Initialise")
    RETURN
999 ERRORSEXITS("DecompositionLine_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE DecompositionLine_Initialise

  !
  !================================================================================================================================
  !

  !>Calculates the lines in the given decomposition topology.
  SUBROUTINE DecompositionTopology_LinesCalculate(decompositionTopology,err,error,*)

    !Argument variables
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology !<A pointer to the decomposition topology to calculate the lines for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables      
    INTEGER(INTG) :: basisLocalLineIdx,basisLocalLineNodeIdx,basisNodeIdx,componentIdx,count,derivativeIdx,elementIdx, &
      & elementLocalNodeIdx,localLineIdx,localNodeIdx,lineEndNodeIdx,lineNumber,maxNumberOfLines,nodeIdx, &
      & nodesInLine(4),numberOfLines,newMaxNumberOfLines,surroundingElement,surroundingElementIdx, &
      & surroundingElementBasisLocalLineIdx,surroundingElementLocalLineIdx,versionIdx
    INTEGER(INTG), ALLOCATABLE :: nodesNumberOfLines(:)
    INTEGER(INTG), POINTER :: tempLines(:,:),newTempLines(:,:)
    REAL(DP) :: approxDimension
    LOGICAL :: found
    TYPE(BasisType), POINTER :: basis,basis2
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionElementType), POINTER :: decompositionElement
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(DecompositionLineType), POINTER :: decompositionLine,decompositionLine2
    TYPE(DecompositionLinesType), POINTER :: decompositionLines
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainElementType), POINTER :: domainElement
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainLineType), POINTER :: domainLine,domainLine2
    TYPE(DomainLinesType), POINTER :: domainLines
    TYPE(DomainNodeType), POINTER :: domainNode
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(MeshType), POINTER :: mesh

    NULLIFY(tempLines)
    NULLIFY(newTempLines)

    ENTERS("DecompositionTopology_LinesCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(decompositionTopology)) CALL FlagError("Decomposition topology is not associated.",err,error,*999)
    NULLIFY(decompositionElements)
    CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology,decompositionElements,err,error,*999)
    NULLIFY(decompositionLines)
    CALL DecompositionTopology_DecompositionLinesGet(decompositionTopology,decompositionLines,err,error,*999)
    NULLIFY(decomposition)
    CALL DecompositionTopology_DecompositionGet(decompositionTopology,decomposition,err,error,*999)
    !Process the mesh component number (component number the decomposition was calculated from) first to establish line
    !topology then process the other mesh components.
    CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*999)
    NULLIFY(domainTopology)
    CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
    NULLIFY(domainNodes)
    CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
    NULLIFY(domainElements)
    CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
    NULLIFY(domainLines)
    CALL DomainTopology_DomainLinesGet(domainTopology,domainLines,err,error,*999)
    
    !Guestimate the number of lines
    SELECT CASE(domain%numberOfDimensions)
    CASE(1)
      maxNumberOfLines=domainElements%totalNumberOfElements
    CASE(2)
      approxDimension=SQRT(REAL(domainElements%totalNumberOfElements,DP))
      !This should give the maximum and will over estimate the number of lines for a "square mesh" by approx 33%
      maxNumberOfLines=NINT(3.0_DP*approxDimension*(approxDimension+1),INTG)
    CASE(3)
      !This should give the maximum and will over estimate the number of lines for a "cube mesh" by approx 73%
      approxDimension=REAL(domainElements%totalNumberOfElements,DP)**(1.0_DP/3.0_DP)
      maxNumberOfLines=NINT(11.0_DP*approxDimension*approxDimension*(approxDimension+1),INTG)
    CASE DEFAULT
      CALL FlagError("Invalid number of dimensions for a domain.",err,error,*999)
    END SELECT
    ALLOCATE(tempLines(4,maxNumberOfLines),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate temporary lines array.",err,error,*999)
    ALLOCATE(nodesNumberOfLines(domainNodes%totalNumberOfNodes),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate nodes number of lines array.",err,error,*999)
    nodesNumberOfLines=0
    numberOfLines=0
    tempLines=0
    !Loop over the elements in the topology
    DO elementIdx=1,domainElements%totalNumberOfElements
      domainElement=>domainElements%elements(elementIdx)
      decompositionElement=>decompositionElements%elements(elementIdx)
      basis=>domainElement%basis
      ALLOCATE(decompositionElement%elementLines(basis%numberOfLocalLines),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate element element lines.",err,error,*999)
      !Loop over the local lines of the element
      DO basisLocalLineIdx=1,basis%numberOfLocalLines
        !Calculate the topology node numbers that make up the line
        nodesInLine=0
        DO basisLocalLineNodeIdx=1,basis%numberOfNodesInLocalLine(basisLocalLineIdx)
          nodesInLine(basisLocalLineNodeIdx)=domainElement%elementNodes( &
            & basis%nodeNumbersInLocalLine(basisLocalLineNodeIdx,basisLocalLineIdx))
        ENDDO !basisLocalLineNodeIdx
        !Try and find a previously created line that matches in the adjacent elements
        found=.FALSE.
        nodeIdx=nodesInLine(1)
        DO surroundingElementIdx=1,domainNodes%nodes(nodeIdx)%numberOfSurroundingElements
          surroundingElement=domainNodes%nodes(nodeIdx)%surroundingElements(surroundingElementIdx)
          IF(surroundingElement/=elementIdx) THEN
            IF(ALLOCATED(decompositionElements%elements(surroundingElement)%elementLines)) THEN
              basis2=>domainElements%elements(surroundingElement)%basis
              DO surroundingElementBasisLocalLineIdx=1,basis2%numberOfLocalLines
                localLineIdx=decompositionElements%elements(surroundingElement)% &
                  & elementLines(surroundingElementBasisLocalLineIdx)
                IF(ALL(nodesInLine(1:basis%numberOfNodesInLocalLine(basisLocalLineIdx))== &
                  & tempLines(1:basis%numberOfNodesInLocalLine(basisLocalLineIdx),localLineIdx))) THEN
                  found=.TRUE.
                  EXIT
                ENDIF
              ENDDO !surroundingElementBasisLocalLineIdx
              IF(found) EXIT
            ENDIF
          ENDIF
        ENDDO !surroundingElementIdx
        IF(found) THEN
          !Line has already been created
          decompositionElement%elementLines(basisLocalLineIdx)=localLineIdx
        ELSE
          !Line has not been created
          IF(numberOfLines==maxNumberOfLines) THEN
            !We are at maximum. Reallocate the LINES array to be 20% bigger and try again.
            newMaxNumberOfLines=NINT(1.20_DP*REAL(maxNumberOfLines,DP),INTG)
            ALLOCATE(newTempLines(4,newMaxNumberOfLines),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate new number of lines.",err,error,*999)
            newTempLines(:,1:numberOfLines)=tempLines(:,1:numberOfLines)
            newTempLines(:,numberOfLines+1:newMaxNumberOfLines)=0
            DEALLOCATE(tempLines)
            tempLines=>newTempLines
            NULLIFY(newTempLines)
            maxNumberOfLines=newMaxNumberOfLines
          ENDIF
          numberOfLines=numberOfLines+1
          tempLines(:,numberOfLines)=nodesInLine
          decompositionElement%elementLines(basisLocalLineIdx)=numberOfLines
          DO basisLocalLineNodeIdx=1,SIZE(nodesInLine,1)
            IF(nodesInLine(basisLocalLineNodeIdx)/=0) &
              & nodesNumberOfLines(nodesInLine(basisLocalLineNodeIdx))= &
              & nodesNumberOfLines(nodesInLine(basisLocalLineNodeIdx))+1
          ENDDO !basisLocalLineNodeIdx
        ENDIF
      ENDDO !basisLocalLineIdx
    ENDDO !elementIdx
    !Allocate the line arrays and set them from the lines and nodeLines arrays
    DO nodeIdx=1,domainNodes%totalNumberOfNodes
      ALLOCATE(domainNodes%nodes(nodeIdx)%nodeLines(nodesNumberOfLines(nodeIdx)),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate node lines array.",err,error,*999)
      domainNodes%nodes(nodeIdx)%numberOfNodeLines=0
    ENDDO !nodeIdx
    DEALLOCATE(nodesNumberOfLines)
    ALLOCATE(decompositionLines%lines(numberOfLines),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate decomposition topology lines.",err,error,*999)
    decompositionLines%numberOfLines=numberOfLines
    ALLOCATE(domainLines%lines(numberOfLines),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate domain topology lines.",err,error,*999)
    domainLines%numberOfLines=numberOfLines
    DO localLineIdx=1,domainLines%numberOfLines
      CALL DecompositionLine_Initialise(decompositionLines%lines(localLineIdx),err,error,*999)
      CALL DomainLine_Initialise(domainLines%lines(localLineIdx),err,error,*999)
      DO basisLocalLineNodeIdx=1,SIZE(tempLines,1)
        IF(tempLines(basisLocalLineNodeIdx,localLineIdx)/=0) THEN
          nodeIdx=tempLines(basisLocalLineNodeIdx,localLineIdx)
          domainNodes%nodes(nodeIdx)%numberOfNodeLines=domainNodes%nodes(nodeIdx)%numberOfNodeLines+1
          domainNodes%nodes(nodeIdx)%nodeLines(domainNodes%nodes(nodeIdx)%numberOfNodeLines)=localLineIdx
        ENDIF
      ENDDO !basisLocalLineNodeIdx  
    ENDDO !localLineIdx
    DO elementIdx=1,decompositionElements%totalNumberOfElements
      decompositionElement=>decompositionElements%elements(elementIdx)
      domainElement=>domainElements%elements(elementIdx)
      basis=>domainElement%basis
      DO basisLocalLineIdx=1,basis%numberOfLocalLines
        lineNumber=decompositionElement%elementLines(basisLocalLineIdx)
        decompositionLine=>decompositionLines%lines(lineNumber)
        domainLine=>domainLines%lines(lineNumber)
        decompositionLine%numberOfSurroundingElements=decompositionLine%numberOfSurroundingElements+1
        IF(.NOT.ASSOCIATED(domainLine%basis)) THEN
          decompositionLine%number=lineNumber
          domainLine%number=lineNumber
          domainLine%elementNumber=elementIdx !Needs checking
          decompositionLine%xiDirection=basis%localLineXiDirection(basisLocalLineIdx)
          IF(ALLOCATED(basis%localLineBasis)) THEN
            domainLine%basis=>basis%localLineBasis(basisLocalLineIdx)%ptr
          ELSE
            !Basis is only 1D
            domainLine%basis=>basis
          ENDIF
          ALLOCATE(domainLine%nodesInLine(basis%numberOfNodesInLocalLine(basisLocalLineIdx)),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate line nodes in line.",err,error,*999)
          ALLOCATE(domainLine%derivativesInLine(2,domainLine%basis%maximumNumberOfDerivatives, &
            & basis%numberOfNodesInLocalLine(basisLocalLineIdx)),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate line derivatives in line.",err,error,*999)
          domainLine%nodesInLine(1:basis%numberOfNodesInLocalLine(basisLocalLineIdx))= &
            & tempLines(1:basis%numberOfNodesInLocalLine(basisLocalLineIdx),lineNumber)
          DO basisLocalLineNodeIdx=1,basis%numberOfNodesInLocalLine(basisLocalLineIdx)
            !Set derivative number of u (NO_GLOBAL_DERIV) for the domain line
            domainLine%derivativesInLine(1,1,basisLocalLineNodeIdx)=NO_GLOBAL_DERIV
            !Set version number of u (NO_GLOBAL_DERIV) for the domain line
            versionIdx=domainElement%elementVersions(1,basis%nodeNumbersInLocalLine(basisLocalLineNodeIdx,basisLocalLineIdx))
            domainLine%derivativesInLine(2,1,basisLocalLineNodeIdx)=versionIdx
            IF(domainLine%basis%maximumNumberOfDerivatives>1) THEN
              derivativeIdx=domainElement%elementDerivatives( &
                & basis%derivativeNumbersInLocalLine(basisLocalLineNodeIdx,basisLocalLineIdx), &
                & basis%nodeNumbersInLocalLine(basisLocalLineNodeIdx,basisLocalLineIdx))
              domainLine%derivativesInLine(1,2,basisLocalLineNodeIdx)=derivativeIdx
              versionIdx=domainElement%elementVersions(basis%derivativeNumbersInLocalLine( &
                & basisLocalLineNodeIdx,basisLocalLineIdx),basis%nodeNumbersInLocalLine( &
                & basisLocalLineNodeIdx,basisLocalLineIdx))
              domainLine%derivativesInLine(2,2,basisLocalLineNodeIdx)=versionIdx
            ENDIF
          ENDDO !basisLocalLineNodeIdx
        ENDIF
      ENDDO !basisLocalLineIdx
    ENDDO !elementIdx
    DEALLOCATE(tempLines)
    !Calculate adjacent lines and the surrounding elements for each line
    DO localLineIdx=1,decompositionLines%numberOfLines
      decompositionLine=>decompositionLines%lines(localLineIdx)
      domainLine=>domainLines%lines(localLineIdx)
      basis=>domainLine%basis
      IF(decompositionLine%numberOfSurroundingElements==1) THEN
        decompositionLine%boundaryLine=.TRUE.
        domainLine%boundaryLine=.TRUE.
      ENDIF
      !Allocate the elements surrounding the line
      ALLOCATE(decompositionLine%surroundingElements(decompositionLine%numberOfSurroundingElements),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate line surrounding elements.",err,error,*999)
      ALLOCATE(decompositionLine%elementLines(decompositionLine%numberOfSurroundingElements),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate line element lines.",err,error,*999)
      decompositionLine%numberOfSurroundingElements=0
      decompositionLine%adjacentLines=0
      !Loop over the nodes at each end of the line
      DO lineEndNodeIdx=0,1
        found=.FALSE.
        nodeIdx=domainLine%nodesInLine(lineEndNodeIdx*(basis%numberOfNodes-1)+1)
        !Loop over the elements surrounding the node.
        DO surroundingElementIdx=1,domainNodes%nodes(nodeIdx)%numberOfSurroundingElements
          surroundingElement=domainNodes%nodes(nodeIdx)%surroundingElements(surroundingElementIdx)
          decompositionElement=>decompositionElements%elements(surroundingElement)
          domainElement=>domainElements%elements(surroundingElement)
          !Loop over the local lines of the element
          DO basisLocalLineIdx=1,domainElement%basis%numberOfLocalLines
            surroundingElementLocalLineIdx=decompositionElement%elementLines(basisLocalLineIdx)
            IF(surroundingElementLocalLineIdx/=localLineIdx) THEN
              decompositionLine2=>decompositionLines%lines(surroundingElementLocalLineIdx)
              domainLine2=>domainLines%lines(surroundingElementLocalLineIdx)
              IF(decompositionLine2%xiDirection==decompositionLine%xiDirection) THEN
                !Lines run in the same direction.
                basis2=>domainLine2%basis
                IF(lineEndNodeIdx==0) THEN
                  localNodeIdx=domainLine2%nodesInLine(basis2%numberOfNodes)
                ELSE
                  localNodeIdx=domainLine2%nodesInLine(1)
                ENDIF
                IF(localNodeIdx==nodeIdx) THEN
                  !The node at the 'other' end of this line matches the node at the current end of the line.
                  !Check it is not a coexistant line running the other way
                  IF(basis2%interpolationOrder(1)==basis%interpolationOrder(1)) THEN
                    count=0
                    DO basisNodeIdx=1,basis%numberOfNodes
                      IF(domainLine2%nodesInLine(basisNodeIdx)== &
                        & domainLine%nodesInLine(basis2%numberOfNodes-basisNodeIdx+1)) &
                        & count=count+1
                    ENDDO !basisNodeIdx
                    IF(count<basis%numberOfNodes) THEN
                      found=.TRUE.
                      EXIT
                    ENDIF
                  ELSE
                    found=.TRUE.
                    EXIT
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
          ENDDO !basisLocalLineIdx
          IF(found) EXIT
        ENDDO !surroundingElementIdx
        IF(found) decompositionLine%adjacentLines(lineEndNodeIdx)=surroundingElementLocalLineIdx
      ENDDO !lineEndNodeIdx
    ENDDO !localLineIdx
    !Set the surrounding elements
    DO elementIdx=1,decompositionElements%totalNumberOfElements
      decompositionElement=>decompositionElements%elements(elementIdx)
      domainElement=>domainElements%elements(elementIdx)
      basis=>domainElement%basis
      DO basisLocalLineIdx=1,basis%numberOfLocalLines
        lineNumber=decompositionElement%elementLines(basisLocalLineIdx)
        decompositionLine=>decompositionLines%lines(lineNumber)
        decompositionLine%numberOfSurroundingElements=decompositionLine%numberOfSurroundingElements+1
        decompositionLine%surroundingElements(decompositionLine%numberOfSurroundingElements)=elementIdx
        decompositionLine%elementLines(decompositionLine%numberOfSurroundingElements)=basisLocalLineIdx
      ENDDO !basisLocalLineIdx
    ENDDO !elementIdx
    !Now loop over the other mesh components in the decomposition and calculate the domain lines
    NULLIFY(mesh)
    CALL Decomposition_MeshGet(decomposition,mesh,err,error,*999)
    DO componentIdx=1,mesh%numberOfComponents
      IF(componentIdx/=decomposition%meshComponentNumber) THEN
        NULLIFY(domain)
        CALL Decomposition_DomainGet(decomposition,componentIdx,domain,err,error,*999)
        NULLIFY(domainTopology)
        CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
        NULLIFY(domainNodes)
        CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
        NULLIFY(domainElements)
        CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
        NULLIFY(domainLines)
        CALL DomainTopology_DomainLinesGet(domainTopology,domainLines,err,error,*999)
        ALLOCATE(domainLines%lines(decompositionLines%numberOfLines),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate domain lines lines.",err,error,*999)
        domainLines%numberOfLines=decompositionLines%numberOfLines
        ALLOCATE(nodesNumberOfLines(domainNodes%totalNumberOfNodes),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate nodes number of lines array.",err,error,*999)
        nodesNumberOfLines=0
        !Loop over the lines in the topology
        DO localLineIdx=1,decompositionLines%numberOfLines
          decompositionLine=>decompositionLines%lines(localLineIdx)
          domainLine=>domainLines%lines(localLineIdx)
          IF(decompositionLine%numberOfSurroundingElements>0) THEN
            elementIdx=decompositionLine%surroundingElements(1)
            basisLocalLineIdx=decompositionLine%elementLines(1)
            CALL DomainLine_Initialise(domainLines%lines(localLineIdx),err,error,*999)
            domainLine%NUMBER=localLineIdx
            domainElement=>domainElements%elements(elementIdx)
            basis=>domainElement%basis
            domainLine%elementNumber=domainElement%NUMBER
            IF(ALLOCATED(basis%localLineBasis)) THEN
              domainLine%basis=>basis%localLineBasis(basisLocalLineIdx)%ptr
            ELSE
              !Basis is only 1D
              domainLine%basis=>basis
            ENDIF
            ALLOCATE(domainLine%nodesInLine(basis%numberOfNodesInLocalLine(basisLocalLineIdx)),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate nodes in line.",err,error,*999)
            ALLOCATE(domainLine%derivativesInLine(2,domainLine%basis%maximumNumberOfDerivatives, &
              & basis%numberOfNodesInLocalLine(basisLocalLineIdx)),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate derivatives in line.",err,error,*999)
            DO basisLocalLineNodeIdx=1,basis%numberOfNodesInLocalLine(basisLocalLineIdx)
              elementLocalNodeIdx=basis%nodeNumbersInLocalLine(basisLocalLineNodeIdx,basisLocalLineIdx)
              nodeIdx=domainElement%elementNodes(elementLocalNodeIdx)
              domainLine%nodesInLine(basisLocalLineNodeIdx)=nodeIdx
              !Set derivative number of u (NO_GLOBAL_DERIV) for the domain line
              domainLine%derivativesInLine(1,1,basisLocalLineNodeIdx)=NO_GLOBAL_DERIV
              !Set version number of u (NO_GLOBAL_DERIV) for the domain line
              versionIdx=domainElement%elementVersions(1,basis%nodeNumbersInLocalLine(basisLocalLineNodeIdx,basisLocalLineIdx))
              domainLine%derivativesInLine(2,1,basisLocalLineNodeIdx)=versionIdx
              IF(domainLine%basis%maximumNumberOfDerivatives>1) THEN
                derivativeIdx=domainElement%elementDerivatives(basis%derivativeNumbersInLocalLine( &
                  & basisLocalLineNodeIdx,basisLocalLineIdx),elementLocalNodeIdx)
                domainLine%derivativesInLine(1,2,basisLocalLineNodeIdx)=derivativeIdx
                versionIdx=domainElement%elementVersions(basis%derivativeNumbersInLocalLine( &
                  & basisLocalLineNodeIdx,basisLocalLineIdx),basis%nodeNumbersInLocalLine( &
                  & basisLocalLineNodeIdx,basisLocalLineIdx))
                domainLine%derivativesInLine(2,2,basisLocalLineNodeIdx)=versionIdx
              ENDIF
              nodesNumberOfLines(nodeIdx)=nodesNumberOfLines(nodeIdx)+1
            ENDDO !basisLocalLineNodeIdx
          ELSE
            CALL FlagError("Line is not surrounded by any elements?",err,error,*999)
          ENDIF
        ENDDO !localLineIdx
        DO nodeIdx=1,domainNodes%totalNumberOfNodes
          ALLOCATE(domainNodes%nodes(nodeIdx)%nodeLines(nodesNumberOfLines(nodeIdx)),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate node lines.",err,error,*999)
          domainNodes%nodes(nodeIdx)%numberOfNodeLines=0
        ENDDO !nodeIdx
        DEALLOCATE(nodesNumberOfLines)
        DO localLineIdx=1,domainLines%numberOfLines
          domainLine=>domainLines%lines(localLineIdx)
          basis=>domainLine%basis
          DO basisLocalLineNodeIdx=1,basis%numberOfNodes
            nodeIdx=domainLine%nodesInLine(basisLocalLineNodeIdx)
            domainNode=>domainNodes%nodes(nodeIdx)
            domainNode%numberOfNodeLines=domainNode%numberOfNodeLines+1
            domainNode%nodeLines(domainNode%numberOfNodeLines)=localLineIdx
          ENDDO !basisLocalLineNodeIdx
        ENDDO !localLineIdx
      ENDIF
    ENDDO !componentIdx
 
    IF(DIAGNOSTICS1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Decomposition topology lines:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of mesh components = ",mesh%numberOfComponents,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of lines = ",decompositionLines%numberOfLines,err,error,*999)
      DO localLineIdx=1,decompositionLines%numberOfLines
        decompositionLine=>decompositionLines%lines(localLineIdx)
        domainLine=>domainLines%lines(localLineIdx)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Line number = ",decompositionLine%NUMBER,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Xi direction = ",decompositionLine%xiDirection,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of surrounding elements = ", &
          & decompositionLine%numberOfSurroundingElements,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,decompositionLine%numberOfSurroundingElements,4,4, &
          & decompositionLine%surroundingElements,'("      Surrounding elements :",4(X,I8))','(28X,4(X,I8))',err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,decompositionLine%numberOfSurroundingElements,4,4, &
          & decompositionLine%elementLines,'("      Element lines        :",4(X,I8))','(28X,4(X,I8))',err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,2,2,2,decompositionLine%adjacentLines, &
          & '("      Adjacent lines       :",2(X,I8))','(28X,2(X,I8))',err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Boundary line = ",decompositionLine%boundaryLine,err,error,*999)
        DO componentIdx=1,mesh%numberOfComponents
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Mesh component : ",componentIdx,err,error,*999)
          domain=>decomposition%domain(componentIdx)%ptr
          domainLine=>domain%topology%LINES%lines(localLineIdx)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Basis user number = ",domainLine%basis%userNumber, &
            & err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Basis family number = ",domainLine%basis%familyNumber, &
            & err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Basis interpolation type = ",domainLine%basis% &
            & interpolationType(1),err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Basis interpolation order = ",domainLine%basis% &
            & interpolationOrder(1),err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Number of nodes in lines = ",domainLine%basis%numberOfNodes, &
            & err,error,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,domainLine%basis%numberOfNodes,4,4,domainLine%nodesInLine, &
            & '("        Nodes in line        :",4(X,I8))','(30X,4(X,I8))',err,error,*999)
          DO basisLocalLineNodeIdx=1,domainLine%basis%numberOfNodes
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Node : ",basisLocalLineNodeIdx,err,error,*999)
            !/TODO::Loop over local_derivative index so this output makes more sense !<derivativesInLine(i,local_derivativeIdx,localNodeIdx)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1, &
              & domainLine%basis%numberOfDerivatives(basisLocalLineNodeIdx),4,4, &
              & domainLine%derivativesInLine(1,:,basisLocalLineNodeIdx),'("            Derivatives in line  :",4(X,I8))', &
              & '(34X,4(X,I8))',err,error,*999)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1, &
              & domainLine%basis%numberOfDerivatives(basisLocalLineNodeIdx),4,4, &
              & domainLine%derivativesInLine(2,:,basisLocalLineNodeIdx), &
              & '("            Derivatives Versions in line  :",4(X,I8))','(34X,4(X,I8))',err,error,*999)
          ENDDO !basisLocalLineNodeIdx
        ENDDO !componentIdx
      ENDDO !localLineIdx
    ENDIF

    EXITS("DecompositionTopology_LinesCalculate")
    RETURN
999 IF(ASSOCIATED(tempLines)) DEALLOCATE(tempLines)
    IF(ASSOCIATED(newTempLines)) DEALLOCATE(newTempLines)
    IF(ALLOCATED(nodesNumberOfLines)) DEALLOCATE(nodesNumberOfLines)
    ERRORSEXITS("DecompositionTopology_LinesCalculate",err,error)
    RETURN 1   
  END SUBROUTINE DecompositionTopology_LinesCalculate

  !
  !================================================================================================================================
  !

  !>Finalises the lines in the given decomposition topology.
  SUBROUTINE DecompositionTopology_LinesFinalise(decompositionLines,err,error,*)

    !Argument variables
    TYPE(DecompositionLinesType), POINTER :: decompositionLines !<A pointer to the decomposition lines to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: lineIdx

    ENTERS("DecompositionTopology_LinesFinalise",err,error,*999)

    IF(ASSOCIATED(decompositionLines)) THEN
      DO lineIdx=1,SIZE(decompositionLines%lines,1)
        CALL DecompositionLine_Finalise(decompositionLines%lines(lineIdx),err,error,*999)
      ENDDO !lineIdx
      IF(ALLOCATED(decompositionLines%lines)) DEALLOCATE(decompositionLines%lines)
      DEALLOCATE(decompositionLines)
    ENDIF

    EXITS("DecompositionTopology_LinesFinalise")
    RETURN
999 ERRORSEXITS("DecompositionTopology_LinesFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE DecompositionTopology_LinesFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the line data structures for a decomposition topology.
  SUBROUTINE DecompositionTopology_LinesInitialise(decompositionTopology,err,error,*)

    !Argument variables
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology !<A pointer to the decomposition topology to initialise the lines for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("DecompositionTopology_LinesInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(decompositionTopology)) CALL FlagError("Decomposition topology is not associated.",err,error,*998)
    IF(ASSOCIATED(decompositionTopology%lines)) &
      & CALL FlagError("Decomposition already has topology lines associated.",err,error,*998)
      
    ALLOCATE(decompositionTopology%lines,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate topology lines.",err,error,*999)
    decompositionTopology%lines%numberOfLines=0
    decompositionTopology%lines%decompositionTopology=>decompositionTopology

    EXITS("DecompositionTopology_LinesInitialise")
    RETURN
999 CALL DecompositionTopology_LinesFinalise(decompositionTopology%lines,dummyErr,dummyError,*998)
998 ERRORSEXITS("DecompositionTopology_LinesInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE DecompositionTopology_LinesInitialise

  !
  !================================================================================================================================
  !

  !>Finalises a decompositio face in the given decomposition topology and deallocates all memory.
  SUBROUTINE DecompositionFace_Finalise(decompositionFace,err,error,*)

    !Argument variables
    TYPE(DecompositionFaceType) :: decompositionFace !<The decomposition face to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DecompositionFace_Finalise",err,error,*999)

    decompositionFace%number=0
    decompositionFace%xiNormalDirection=0
    decompositionFace%numberOfSurroundingElements=0
    IF(ALLOCATED(decompositionFace%surroundingElements)) DEALLOCATE(decompositionFace%surroundingElements)
    IF(ALLOCATED(decompositionFace%elementFaces)) DEALLOCATE(decompositionFace%elementFaces)
 
    EXITS("DecompositionFace_Finalise")
    RETURN
999 ERRORSEXITS("DecompositionFace_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE DecompositionFace_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the face data structure for a decomposition topology face.
  SUBROUTINE DecompositionFace_Initialise(decompositionFace,err,error,*)

    !Argument variables
    TYPE(DecompositionFaceType) :: decompositionFace !<The decomposition face to initialise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DecompositionFace_Initialise",err,error,*999)

    decompositionFace%number=0
    decompositionFace%xiNormalDirection=0
    decompositionFace%numberOfSurroundingElements=0
    decompositionFace%boundaryFace=.FALSE.

    EXITS("DecompositionFace_Initialise")
    RETURN
999 ERRORSEXITS("DecompositionFace_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE DecompositionFace_Initialise

  !
  !================================================================================================================================
  !

  !>Calculates the faces in the given decomposition topology.
  SUBROUTINE DecompositionTopology_FacesCalculate(decompositionTopology,err,error,*)

    !Argument variables
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology !<A pointer to the decomposition topology to calculate the faces for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: basisLocalFaceIdx,basisLocalFaceDerivativeIdx,basisLocalFaceNodeIdx,componentIdx,derivativeIdx,elementIdx, &
      & elementLocalNodeIdx,faceIdx,faceNumber,nodeIdx,maxNumberOfFaces,nodesInFace(16),numberOfFaces,newMaxNumberOfFaces, &
      & surroundingElement,surroundingElementIdx,surroundingElementBasisLocalFaceIdx,versionIdx
    INTEGER(INTG), ALLOCATABLE :: nodesNumberOfFaces(:)
    INTEGER(INTG), POINTER :: tempFaces(:,:),newTempFaces(:,:)
    LOGICAL :: found
    TYPE(BasisType), POINTER :: basis,basis2
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionElementType), POINTER :: decompositionElement
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(DecompositionFaceType), POINTER :: decompositionFace!,decompositionFace2
    TYPE(DecompositionFacesType), POINTER :: decompositionFaces
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainElementType), POINTER :: domainElement
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainFaceType), POINTER :: domainFace!,domainFace2
    TYPE(DomainFacesType), POINTER :: domainFaces
    TYPE(DomainNodeType), POINTER :: domainNode
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(MeshType), POINTER :: mesh

    NULLIFY(tempFaces)
    NULLIFY(newTempFaces)

    ENTERS("DecompositionTopology_FacesCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(decompositionTopology)) CALL FlagError("Decomposition topology is not associated.",err,error,*999)
    NULLIFY(decompositionElements)
    CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology,decompositionElements,err,error,*999)
    NULLIFY(decompositionFaces)
    CALL DecompositionTopology_DecompositionFacesGet(decompositionTopology,decompositionFaces,err,error,*999)
    NULLIFY(decomposition)
    CALL DecompositionTopology_DecompositionGet(decompositionTopology,decomposition,err,error,*999)
    
    !Process the mesh component number (component number the decomposition was calculated from) first to establish face
    !topology then process the other mesh components.
    NULLIFY(domain)
    CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*999)
    NULLIFY(domainTopology)
    CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
    NULLIFY(domainNodes)
    CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
    NULLIFY(domainElements)
    CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
    NULLIFY(domainFaces)
    CALL DomainTopology_DomainFacesGet(domainTopology,domainFaces,err,error,*999)
        
    !Estimate the number of faces
    SELECT CASE(DOMAIN%numberOfDimensions)
    CASE(1)
      ! Faces not calculated in 1D 
    CASE(2)
      ! Faces not calculated in 2D
    CASE(3)
      !This should give the maximum and will over estimate the number of faces for a "cube mesh" by approx 33%
      maxNumberOfFaces=NINT(((REAL(domainElements%totalNumberOfElements,DP)*5.0_DP)+1.0_DP)*(4.0_DP/3.0_DP),INTG)
      ALLOCATE(tempFaces(16,maxNumberOfFaces),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate temporary faces array",err,error,*999)
      ALLOCATE(nodesNumberOfFaces(domainNodes%totalNumberOfNodes),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate nodes number of faces array",err,error,*999)
      nodesNumberOfFaces=0
      numberOfFaces=0
      tempFaces=0
      !Loop over the elements in the topology and fill temp_faces with node numbers for each element
      DO elementIdx=1,domainElements%totalNumberOfElements
        domainElement=>domainElements%elements(elementIdx)
        decompositionElement=>decompositionElements%elements(elementIdx)
        basis=>domainElement%basis
        ALLOCATE(decompositionElement%elementFaces(basis%numberOfLocalFaces),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate element faces of element",err,error,*999)
        !Loop over the local faces of the element
        DO basisLocalFaceIdx=1,basis%numberOfLocalFaces
          !Calculate the topology node numbers that make up the face
          nodesInFace=0
          !Check whether face has already been read out
          DO basisLocalFaceNodeIdx=1,basis%numberOfNodesInLocalFace(basisLocalFaceIdx)
            !Read out node numbers of local face from ELEMENT_NODES
            nodesInFace(basisLocalFaceNodeIdx)=domainElement%elementNodes( &
              & basis%nodeNumbersInLocalFace(basisLocalFaceNodeIdx,basisLocalFaceIdx))
          ENDDO !basisLocalFaceNodeIdx
          !Try and find a previously created face that matches in the adjacent elements
          found=.FALSE.
          nodeIdx=nodesInFace(1)
          DO surroundingElementIdx=1,domainNodes%nodes(nodeIdx)%numberOfSurroundingElements
            surroundingElement=domainNodes%nodes(nodeIdx)%surroundingElements(surroundingElementIdx)
            IF(surroundingElement/=elementIdx) THEN
              IF(ALLOCATED(decompositionElements%elements(surroundingElement)%elementFaces)) THEN
                basis2=>domainElements%elements(surroundingElement)%basis
                DO surroundingElementBasisLocalFaceIdx=1,basis2%numberOfLocalFaces
                  faceIdx=decompositionElements%elements(surroundingElement)%elementFaces( &
                    & surroundingElementBasisLocalFaceIdx)
                  IF(ALL(nodesInFace(1:basis%numberOfNodesInLocalFace(basisLocalFaceIdx))== &
                    & tempFaces(1:basis%numberOfNodesInLocalFace(basisLocalFaceIdx),faceIdx))) THEN
                    found=.TRUE.
                    EXIT
                  ENDIF
                ENDDO !surroundingElementBasisLocalFaceIdx
                IF(found) EXIT
              ENDIF
            ENDIF
          ENDDO !surroundingElemntIdx
          IF(found) THEN
            !Face has already been created
            decompositionElement%elementFaces(basisLocalFaceIdx)=faceIdx
          ELSE
            !Face has not been created
            IF(numberOfFaces==maxNumberOfFaces) THEN
              !We are at maximum. Reallocate the FACES array to be 20% bigger and try again.
              newMaxNumberOfFaces=NINT(1.20_DP*REAL(maxNumberOfFaces,DP),INTG)
              !\todo: Change 16 to a variable and above for nodesInFace
              ALLOCATE(newTempFaces(16,newMaxNumberOfFaces),STAT=err)
              IF(err/=0) CALL FlagError("Could not allocate new number of faces",err,error,*999)
              newTempFaces(:,1:numberOfFaces)=tempFaces(:,1:numberOfFaces)
              newTempFaces(:,numberOfFaces+1:newMaxNumberOfFaces)=0
              DEALLOCATE(tempFaces)
              tempFaces=>newTempFaces
              NULLIFY(newTempFaces)
              maxNumberOfFaces=newMaxNumberOfFaces
            ENDIF
            numberOfFaces=numberOfFaces+1
            tempFaces(:,numberOfFaces)=nodesInFace(:)
            decompositionElement%elementFaces(basisLocalFaceIdx)=numberOfFaces
            DO basisLocalFaceNodeIdx=1,SIZE(nodesInFace,1)
              IF(nodesInFace(basisLocalFaceNodeIdx)/=0) &
                & nodesNumberOfFaces(nodesInFace(basisLocalFaceNodeIdx))= &
                & nodesNumberOfFaces(nodesInFace(basisLocalFaceNodeIdx))+1
            ENDDO !basisLocalFaceNodeIdx
          ENDIF
        ENDDO !basisLocalFaceIdx
      ENDDO !elementIdx
      
      !Allocate the face arrays and set them from the FACES and nodeFaces arrays
      DO nodeIdx=1,domainNodes%totalNumberOfNodes
        ALLOCATE(domainNodes%nodes(nodeIdx)%nodeFaces(nodesNumberOfFaces(nodeIdx)),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate node faces array",err,error,*999)
        domainNodes%nodes(nodeIdx)%numberOfNodeFaces=0
      ENDDO !nodeIdx
      DEALLOCATE(nodesNumberOfFaces)
      ALLOCATE(decompositionFaces%faces(numberOfFaces),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate decomposition topology faces",err,error,*999)
      decompositionFaces%numberOfFaces=numberOfFaces
      ALLOCATE(domainFaces%faces(numberOfFaces),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate domain topology faces",err,error,*999)
      domainFaces%numberOfFaces=numberOfFaces
      DO faceIdx=1,domainFaces%numberOfFaces
        CALL DecompositionFace_Initialise(decompositionFaces%faces(faceIdx),err,error,*999)
        CALL DomainFace_Initialise(domainFaces%faces(faceIdx),err,error,*999)
        DO basisLocalFaceNodeIdx=1,SIZE(tempFaces,1)
          IF(tempFaces(basisLocalFaceNodeIdx,faceIdx)/=0) THEN
            nodeIdx=tempFaces(basisLocalFaceNodeIdx,faceIdx)
            domainNodes%nodes(nodeIdx)%numberOfNodeFaces=domainNodes%nodes(nodeIdx)%numberOfNodeFaces+1
            domainNodes%nodes(nodeIdx)%nodeFaces(domainNodes%nodes(nodeIdx)%numberOfNodeFaces)=faceIdx
          ENDIF
        ENDDO !basisLocalFaceNodeIdx  
      ENDDO !faceIdx
      
      !Set nodes in face and derivatives of nodes in face for domain faces
      DO elementIdx=1,decompositionElements%totalNumberOfElements
        decompositionElement=>decompositionElements%elements(elementIdx)
        domainElement=>domainElements%elements(elementIdx)
        basis=>domainElement%basis
        !Loop over local faces of element
        DO basisLocalFaceIdx=1,basis%numberOfLocalFaces
          faceNumber=decompositionElement%elementFaces(basisLocalFaceIdx)
          decompositionFace=>decompositionFaces%faces(faceNumber)
          domainFace=>domainFaces%faces(faceNumber)
          decompositionFace%numberOfSurroundingElements=decompositionFace%numberOfSurroundingElements+1
          IF(.NOT.ASSOCIATED(domainFace%basis)) THEN
            decompositionFace%NUMBER=faceNumber
            domainFace%NUMBER=faceNumber
            domainFace%elementNumber=elementIdx !! Needs checking
            !decompositionFace%elementNumber=decompositionElement%number
            !domainFace%elementNumber=domainElement%number
            decompositionFace%xiNormalDirection=basis%localFaceXiNormal(basisLocalFaceIdx)
            IF(ALLOCATED(basis%localFaceBasis)) THEN
              domainFace%basis=>basis%localFaceBasis(basisLocalFaceIdx)%ptr
            ELSE
              !Basis is only 2D
              domainFace%basis=>basis
            ENDIF
            ALLOCATE(domainFace%nodesInFace(basis%numberOfNodesInLocalFace(basisLocalFaceIdx)), &
              & STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate face nodes in face",err,error,*999)
            ALLOCATE(domainFace%derivativesInFace(2,domainFace%basis%maximumNumberOfDerivatives, &
              & basis%numberOfNodesInLocalFace(basisLocalFaceIdx)),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate face derivatives in face",err,error,*999)
            domainFace%derivativesInFace=0
            !Set nodes in face based upon face number
            domainFace%nodesInFace(1:basis%numberOfNodesInLocalFace(basisLocalFaceIdx))= &
              & tempFaces(1:basis%numberOfNodesInLocalFace(basisLocalFaceIdx),faceNumber)
            !Set derivatives of nodes in domain face from derivatives of nodes in element
            DO basisLocalFaceNodeIdx=1,basis%numberOfNodesInLocalFace(basisLocalFaceIdx)
              elementLocalNodeIdx=basis%nodeNumbersInLocalFace(basisLocalFaceNodeIdx, &
                & basisLocalFaceIdx)
              !Set derivative number of u (NO_GLOBAL_DERIV) for the domain face
              domainFace%derivativesInFace(1,1,basisLocalFaceNodeIdx)=NO_GLOBAL_DERIV
              !Set version number of u (NO_GLOBAL_DERIV) for the domain face
              versionIdx=domainElement%elementVersions(1,basis%nodeNumbersInLocalFace( &
                & basisLocalFaceNodeIdx,basisLocalFaceIdx))
              domainFace%derivativesInFace(2,1,basisLocalFaceNodeIdx)=versionIdx
              IF(domainFace%basis%maximumNumberOfDerivatives>1) THEN
                DO basisLocalFaceDerivativeIdx=2,domainFace%basis%maximumNumberOfDerivatives
                  derivativeIdx=domainElement%elementDerivatives(basis%derivativeNumbersInLocalFace( &
                    & basisLocalFaceDerivativeIdx,basisLocalFaceNodeIdx,basisLocalFaceIdx),elementLocalNodeIdx)
                  domainFace%derivativesInFace(1,basisLocalFaceDerivativeIdx,basisLocalFaceNodeIdx)=derivativeIdx
                  versionIdx=domainElement%elementVersions(basis%derivativeNumbersInLocalFace( &
                    & basisLocalFaceDerivativeIdx,basisLocalFaceNodeIdx,basisLocalFaceIdx),elementLocalNodeIdx)
                  domainFace%derivativesInFace(2,basisLocalFaceDerivativeIdx,basisLocalFaceNodeIdx)=versionIdx
                ENDDO !basisLocalFaceDerivativeIdx
              ENDIF
            ENDDO !basisLocalFaceNodeIdx
          ENDIF
        ENDDO !basisLocalFaceIdx
      ENDDO !elementIdx
      
      DEALLOCATE(tempFaces)
      !\todo Note: Adjacency will be left out of faces calculation for the time being
      !Calculate adjacent faces and the surrounding elements for each face
      DO faceIdx=1,decompositionFaces%numberOfFaces
        decompositionFace=>decompositionFaces%faces(faceIdx)
        domainFace=>domainFaces%faces(faceIdx)
        basis=>domainFace%basis
        IF(decompositionFace%numberOfSurroundingElements==1) THEN
          decompositionFace%boundaryFace=.TRUE.
          domainFace%boundaryFace=.TRUE.
        ENDIF
        !Allocate the elements surrounding the face
        ALLOCATE(decompositionFace%surroundingElements(decompositionFace%numberOfSurroundingElements),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate face surrounding elements",err,error,*999)          
        ALLOCATE(decompositionFace%elementFaces(decompositionFace%numberOfSurroundingElements),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate face element faces",err,error,*999)
        
        !decompositionFace%numberOfSurroundingElements=0
        !decompositionFace%ADJACENT_FACES=0
        !Loop over the nodes at each end of the face
        !DO nodeIdx1=0,1
        !  DO nodeIdx2=0,1
        !    found=.FALSE.
        !    nodeIdx=domainFace%nodesInFace((nodeIdx2*basis%numberOfNodesInXiDirection*(basis%numberOfFaces-1))&
        !      &+(nodeIdx1*(basis%numberOfNodesInXiDirection-1))+1)
        !    !Loop over the elements surrounding the node.
        !    DO surroundingElementIdx=1,domainNodes%nodes(nodeIdx)%numberOfSurroundingElements
        !      elementIdx=domainNodes%nodes(nodeIdx)%surroundingElements(surroundingElementIdx)
        !      decompositionElement=>decompositionElements%elements(elementIdx)
        !      domainElement=>domainElements%elements(elementIdx)
        !      !Loop over the local faces of the element
        !      DO basisLocalFaceIdx=1,domainElement%basis%numberOfLocalFaces
        !        faceIdx2=decompositionElement%elementFaces(basisLocalFaceIdx)
        !        IF(faceIdx2/=faceIdx) THEN
        !          decompositionFace2=>decompositionFaces%faces(faceIdx2)
        !          domainFace2=>domainFaces%faces(faceIdx2)
        !          !Check whether XI of face have same direction
        !          IF ((OTHER_XI_DIRECTIONS3(basis%localFaceXiNormal(basisLocalFaceIdx),2,1)==&
        !            &OTHER_XI_DIRECTIONS3(basis2%localFaceXiNormal(basisLocalFaceIdx),2,1)).OR.&
        !            &(OTHER_XI_DIRECTIONS3(basis%localFaceXiNormal(basisLocalFaceIdx),3,1)==&
        !            &OTHER_XI_DIRECTIONS3(basis2%localFaceXiNormal(basisLocalFaceIdx),3,1))) THEN
        !            !Loop over nodes in face of surrounding element
        !            basis2=>domainFace2%basis
        !            IF(basis2%interpolationOrder(1)==basis%interpolationOrder(1)) THEN
        !              nodeCount=0
        !              DO nodeIdx3=1,basis%numberOfNodesInXiDirection
        !                DO nodeIdx4=1,basis%numberOfNodesInXiDirection
        !                  nodeNumber2=domainFace2%nodesInFace((nodeIdx4*(basis2%numberOfFaces-1))&
        !                    &+(nodeIdx3*(basis2%numberOfNodesInXiDirection-1))+1)
        !                  IF(nodeNumber2==nodeIdx) nodeCount=nodeCount+1
        !                ENDDO !nodeIdx4
        !              ENDDO !nodeIdx3
        !              IF(nodeCount<basis%numberOfNodes) THEN
        !                found=.TRUE.
        !                EXIT
        !              ENDIF
        !            ENDIF
        !          ENDIF
        !        ENDIF
        !      ENDDO !basisLocalFaceIdx
        !      IF(found) EXIT
        !    ENDDO !surroundingElementIdx
        !    IF(found) decompositionFace%ADJACENT_FACES(nodeIdx2)=faceIdx2
        !  ENDDO !nodeIdx2
        !  IF(found) decompositionFace%ADJACENT_FACES(nodeIdx1)=faceIdx2
        !ENDDO !nodeIdx1
                                                   
      ENDDO !faceIdx
      
      !Set the surrounding elements
      DO elementIdx=1,decompositionElements%totalNumberOfElements
        decompositionElement=>decompositionElements%elements(elementIdx)
        domainElement=>domainElements%elements(elementIdx)
        basis=>domainElement%basis
        DO basisLocalFaceIdx=1,basis%numberOfLocalFaces
          faceNumber=decompositionElement%elementFaces(basisLocalFaceIdx)
          decompositionFace=>decompositionFaces%faces(faceNumber)
          DO faceIdx=1,decompositionFace%numberOfSurroundingElements
            decompositionFace%surroundingElements(faceIdx)=elementIdx
            decompositionFace%elementFaces(faceIdx)=basisLocalFaceIdx
          ENDDO
        ENDDO !basisLocalFaceIdx
      ENDDO !elementIdx
    CASE DEFAULT
      CALL FlagError("Invalid number of dimensions for a topology domain",err,error,*999)
    END SELECT
    
    !Now loop over the other mesh components in the decomposition and calculate the domain faces
    NULLIFY(mesh)
    CALL Decomposition_MeshGet(decomposition,mesh,err,error,*999)
    DO componentIdx=1,mesh%numberOfComponents
      IF(componentIdx/=decomposition%meshComponentNumber) THEN
        NULLIFY(domain)
        CALL Decomposition_DomainGet(decomposition,componentIdx,domain,err,error,*999)
        NULLIFY(domainTopology)
        CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
        NULLIFY(domainNodes)
        CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
        NULLIFY(domainElements)
        CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
        NULLIFY(domainFaces)
        CALL DomainTopology_DomainFacesGet(domainTopology,domainFaces,err,error,*999)
        ALLOCATE(domainFaces%faces(decompositionFaces%numberOfFaces),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate domain faces faces",err,error,*999)
        domainFaces%numberOfFaces=decompositionFaces%numberOfFaces
        ALLOCATE(nodesNumberOfFaces(domainNodes%totalNumberOfNodes),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate nodes number of faces array",err,error,*999)
        nodesNumberOfFaces=0
        !Loop over the faces in the topology
        DO faceIdx=1,decompositionFaces%numberOfFaces
          decompositionFace=>decompositionFaces%faces(faceIdx)
          domainFace=>domainFaces%faces(faceIdx)
          IF(decompositionFace%numberOfSurroundingElements>0) THEN
            elementIdx=decompositionFace%surroundingElements(1)
            basisLocalFaceIdx=decompositionFace%elementFaces(1)
            CALL DomainFace_Initialise(domainFaces%faces(faceIdx),err,error,*999)
            domainFace%NUMBER=faceIdx
            domainElement=>domainElements%elements(elementIdx)
            basis=>domainElement%basis
            IF(ALLOCATED(basis%localFaceBasis)) THEN
              domainFace%basis=>basis%localFaceBasis(basisLocalFaceIdx)%ptr
            ELSE
              !Basis is only 2D
              domainFace%basis=>basis
            ENDIF
            ALLOCATE(domainFace%nodesInFace(basis%numberOfNodesInLocalFace(basisLocalFaceIdx)),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate nodes in face",err,error,*999)
            ALLOCATE(domainFace%derivativesInFace(2,domainFace%basis%maximumNumberOfDerivatives, &
              & basis%numberOfNodesInLocalFace(basisLocalFaceIdx)),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate derivatives in face",err,error,*999)
            !Set derivatives of nodes in domain face from derivatives of nodes in element
            DO basisLocalFaceNodeIdx=1,basis%numberOfNodesInLocalFace(basisLocalFaceIdx)
              elementLocalNodeIdx=basis%nodeNumbersInLocalFace(basisLocalFaceNodeIdx,basisLocalFaceIdx)
              nodeIdx=domainElement%elementNodes(elementLocalNodeIdx)
              domainFace%nodesInFace(basisLocalFaceNodeIdx)=nodeIdx
              !Set derivative number of u (NO_GLOBAL_DERIV) for the domain face
              domainFace%derivativesInFace(1,1,basisLocalFaceNodeIdx)=NO_GLOBAL_DERIV
              !Set version number of u (NO_GLOBAL_DERIV) for the domain face
              versionIdx=domainElement%elementVersions(1,basis%nodeNumbersInLocalFace( &
                & basisLocalFaceNodeIdx,basisLocalFaceIdx))
              domainFace%derivativesInFace(2,1,basisLocalFaceNodeIdx)=versionIdx
              IF(domainFace%basis%maximumNumberOfDerivatives>1) THEN
                DO basisLocalFaceDerivativeIdx=2,domainFace%basis%maximumNumberOfDerivatives
                  derivativeIdx=domainElement%elementDerivatives(basis%derivativeNumbersInLocalFace( &
                    & basisLocalFaceDerivativeIdx,basisLocalFaceNodeIdx,basisLocalFaceIdx),elementLocalNodeIdx)
                  domainFace%derivativesInFace(1,basisLocalFaceDerivativeIdx,basisLocalFaceNodeIdx)=derivativeIdx
                  versionIdx=domainElement%elementVersions(basis%derivativeNumbersInLocalFace( &
                    & basisLocalFaceDerivativeIdx,basisLocalFaceNodeIdx,basisLocalFaceIdx),elementLocalNodeIdx)
                  domainFace%derivativesInFace(2,basisLocalFaceDerivativeIdx,basisLocalFaceNodeIdx)=versionIdx
                ENDDO !basisLocalFaceDerivativeIdx
              ENDIF
              nodesNumberOfFaces(nodeIdx)=nodesNumberOfFaces(nodeIdx)+1
            ENDDO !basisLocalFaceNodeIdx
          ELSE
            CALL FlagError("Face is not surrounded by any elements?",err,error,*999)
          ENDIF
        ENDDO !faceIdx
        DO nodeIdx=1,domainNodes%totalNumberOfNodes
          ALLOCATE(domainNodes%nodes(nodeIdx)%nodeFaces(nodesNumberOfFaces(nodeIdx)),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate node faces",err,error,*999)
          domainNodes%nodes(nodeIdx)%numberOfNodeFaces=0
        ENDDO !nodeIdx
        DEALLOCATE(nodesNumberOfFaces)
        DO faceIdx=1,domainFaces%numberOfFaces
          domainFace=>domainFaces%faces(faceIdx)
          basis=>domainFace%basis
          DO basisLocalFaceNodeIdx=1,basis%numberOfNodes
            nodeIdx=domainFace%nodesInFace(basisLocalFaceNodeIdx)
            domainNode=>domainNodes%nodes(nodeIdx)
            domainNode%numberOfNodeFaces=domainNode%numberOfNodeFaces+1
            !Set the face numbers a node is on
            domainNode%nodeFaces(domainNode%numberOfNodeFaces)=faceIdx
          ENDDO !basisLocalFaceNodeIdx
        ENDDO !faceIdx
      ENDIF
    ENDDO !componentIdx

    IF(DIAGNOSTICS1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Decomposition topology faces:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of mesh components = ",mesh%numberOfComponents,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of faces = ",decompositionFaces%numberOfFaces,err,error,*999)
      DO faceIdx=1,decompositionFaces%numberOfFaces
        decompositionFace=>decompositionFaces%faces(faceIdx)
        domainFace=>domainFaces%faces(faceIdx)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Face number = ",decompositionFace%NUMBER,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Xi direction (Normal to Face) = &
          &",decompositionFace%xiNormalDirection,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of surrounding elements = ", &
          & decompositionFace%numberOfSurroundingElements,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,decompositionFace%numberOfSurroundingElements,4,4, &
          & decompositionFace%surroundingElements,'("      Surrounding elements :",4(X,I8))','(28X,4(X,I8))',err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,decompositionFace%numberOfSurroundingElements,4,4, &
          & decompositionFace%elementFaces,'("      Element faces        :",4(X,I8))','(28X,4(X,I8))',err,error,*999)
        !CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,2,2,2,decompositionFace%ADJACENT_FACES, &
        !  & '("      Adjacent faces       :",2(X,I8))','(28X,2(X,I8))',err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Boundary face = ",decompositionFace%boundaryFace,err,error,*999)
        DO componentIdx=1,mesh%numberOfComponents
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Mesh component : ",componentIdx,err,error,*999)
          domain=>decomposition%domain(componentIdx)%ptr
          domainFace=>domain%topology%faces%faces(faceIdx)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Basis user number = ",domainFace%basis%userNumber, &
            & err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Basis family number = ",domainFace%basis%familyNumber, &
            & err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Basis interpolation type = ",domainFace%basis% &
            & interpolationType(1),err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Basis interpolation order = ",domainFace%basis% &
            & interpolationOrder(1),err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Number of nodes in faces = ",domainFace%basis%numberOfNodes, &
            & err,error,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,domainFace%basis%numberOfNodes,4,4,domainFace%nodesInFace, &
            & '("        Nodes in face        :",4(X,I8))','(30X,4(X,I8))',err,error,*999)
          DO basisLocalFaceNodeIdx=1,domainFace%basis%numberOfNodes
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"          Node : ",basisLocalFaceNodeIdx,err,error,*999)
            !/TODO::Loop over local_derivative index so this output makes more sense !<derivativesInLine(i,local_derivativeIdx,local_node_idx)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1, &
              & domainFace%basis%numberOfDerivatives(basisLocalFaceNodeIdx),4,4,domainFace% &
              & derivativesInFace(1,:,basisLocalFaceNodeIdx),'("            Derivatives in face  :",4(X,I8))', &
              & '(34X,4(X,I8))',err,error,*999)
            CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1, &
              & domainFace%basis%numberOfDerivatives(basisLocalFaceNodeIdx),4,4,domainFace% &
              & derivativesInFace(2,:,basisLocalFaceNodeIdx),'("            Derivatives Versions in face  :",4(X,I8))', &
              & '(34X,4(X,I8))',err,error,*999)
          ENDDO !basisLocalFaceNodeIdx
        ENDDO !componentIdx
      ENDDO !faceIdx
    ENDIF

    EXITS("DecompositionTopology_FacesCalculate")
    RETURN
999 IF(ASSOCIATED(tempFaces)) DEALLOCATE(tempFaces)
    IF(ASSOCIATED(newTempFaces)) DEALLOCATE(newTempFaces)
    IF(ALLOCATED(nodesNumberOfFaces)) DEALLOCATE(nodesNumberOfFaces)
    ERRORSEXITS("DecompositionTopology_FacesCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE DecompositionTopology_FacesCalculate

  !
  !================================================================================================================================
  !

  !>Finalises the faces in the given decomposition topology.
  SUBROUTINE DecompositionTopology_FacesFinalise(decompositionFaces,err,error,*)

    !Argument variables
    TYPE(DecompositionFacesType), POINTER :: decompositionFaces !<A pointer to the decomposition faces to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: faceIdx

    ENTERS("DecompositionTopology_FacesFinalise",err,error,*999)

    IF(ASSOCIATED(decompositionFaces)) THEN
      DO faceIdx=1,SIZE(decompositionFaces%faces,1)
        CALL DecompositionFace_Finalise(decompositionFaces%faces(faceIdx),err,error,*999)
      ENDDO !faceIdx
      IF(ALLOCATED(decompositionFaces%faces)) DEALLOCATE(decompositionFaces%faces)
      DEALLOCATE(decompositionFaces)
    ENDIF

    EXITS("DecompositionTopology_FacesFinalise")
    RETURN
999 ERRORSEXITS("DecompositionTopology_FacesFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE DecompositionTopology_FacesFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the face data structures for a decomposition topology.
  SUBROUTINE DecompositionTopology_FacesInitialise(decompositionTopology,err,error,*)

    !Argument variables
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology !<A pointer to the decomposition topology to initialise the faces for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("DecompositionTopology_FacesInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(decompositionTopology)) CALL FlagError("Decomposition topology is not associated.",err,error,*998)
    IF(ASSOCIATED(decompositionTopology%faces)) &
      & CALL FlagError("Decomposition already has topology faces associated",err,error,*998)
      
    ALLOCATE(decompositionTopology%faces,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate topology faces",err,error,*999)
    decompositionTopology%faces%numberOfFaces=0
    decompositionTopology%faces%decompositionTopology=>decompositionTopology

    EXITS("DecompositionTopology_FacesInitialise")
    RETURN
999 CALL DecompositionTopology_FacesFinalise(decompositionTopology%faces,dummyErr,dummyError,*998)
998 ERRORSEXITS("DecompositionTopology_FacesInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE DecompositionTopology_FacesInitialise

  !
  !================================================================================================================================
  !

  !>Gets the decomposition type for a decomposition. \see OpenCMISS::Iron::cmfe_Decomposition_TypeGet
  SUBROUTINE Decomposition_TypeGet(decomposition,type,err,error,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer to the decomposition to get the type for
    INTEGER(INTG), INTENT(OUT) :: TYPE !<On return, the decomposition type for the specified decomposition \see MESH_ROUTINES_DecompositionTypes,MESH_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Decomposition_TypeGet",err,error,*999)

    CALL Decomposition_AssertIsFinished(decomposition,err,error,*999)
    
    type=decomposition%domainDecompositionType
 
    EXITS("Decomposition_TypeGet")
    RETURN
999 ERRORSEXITS("Decomposition_TypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_TypeGet

  !
  !================================================================================================================================
  !

  !>Sets/changes the decomposition type for a decomposition.  \see OpenCMISS::Iron::cmfe_Decomposition_TypeSet
  SUBROUTINE Decomposition_TypeSet(decomposition,type,err,error,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer to the decomposition to set the type for
    INTEGER(INTG), INTENT(IN) :: TYPE !<The decomposition type to set \see DecompositionRoutines_DecompositionTypes,DecompositionRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Decomposition_TypeSet",err,error,*999)

    CALL Decomposition_AssertNotFinished(decomposition,err,error,*999)
    
    SELECT CASE(TYPE)
    CASE(DECOMPOSITION_ALL_TYPE)
      decomposition%domainDecompositionType=DECOMPOSITION_ALL_TYPE
    CASE(DECOMPOSITION_CALCULATED_TYPE)
      decomposition%domainDecompositionType=DECOMPOSITION_CALCULATED_TYPE
    CASE(DECOMPOSITION_USER_DEFINED_TYPE)
      decomposition%domainDecompositionType=DECOMPOSITION_USER_DEFINED_TYPE
    CASE DEFAULT
      localError="Decomposition type "//TRIM(NumberToVString(type,"*",err,error))//" is not valid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Decomposition_TypeSet")
    RETURN
999 ERRORSEXITS("Decomposition_TypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_TypeSet

  !
  !================================================================================================================================
  !

  !>Sets/changes whether lines should be calculated in the the decomposition. \see OpenCMISS::Iron::cmfe_Decomposition_CalculateLinesSet
  SUBROUTINE Decomposition_CalculateLinesSet(decomposition,calculateLinesFlag,err,error,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer to the decomposition
    LOGICAL, INTENT(IN) :: calculateLinesFlag !<The boolean flag to determine whether the lines should be calculated or not
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("Decomposition_CalculateLinesSet",err,error,*999)

    CALL Decomposition_AssertNotFinished(decomposition,err,error,*999)
    
    decomposition%calculateLines=calculateLinesFlag
 
    EXITS("Decomposition_CalculateLinesSet")
    RETURN
999 ERRORSEXITS("Decomposition_CalculateLinesSet",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_CalculateLinesSet

  !
  !================================================================================================================================
  !

  !>Sets/changes whether faces should be calculated in the the decomposition. \see OpenCMISS::Iron::cmfe_Decomposition_CalculateFacesSet
  SUBROUTINE Decomposition_CalculateFacesSet(decomposition,calculateFacesFlag,err,error,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer to the decomposition
    LOGICAL, INTENT(IN) :: calculateFacesFlag !<The boolean flag to determine whether the faces should be calculated or not
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("Decomposition_CalculateFacesSet",err,error,*999)

    CALL Decomposition_AssertNotFinished(decomposition,err,error,*999)
    
    decomposition%calculateFaces=calculateFacesFlag

    EXITS("Decomposition_CalculateFacesSet")
    RETURN
999 ERRORSEXITS("Decomposition_CalculateFacesSet",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_CalculateFacesSet

  !
  !================================================================================================================================
  !

  !>Finalises the domain for a given decomposition and deallocates all memory.
  SUBROUTINE Domain_Finalise(domain,err,error,*)

    !Argument variables
    TYPE(DomainType), POINTER :: domain !<A pointer to the decomposition to finalise the domain for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("Domain_Finalise",err,error,*999)

    IF(ASSOCIATED(domain)) THEN
      IF(ALLOCATED(domain%nodeDomain)) DEALLOCATE(domain%nodeDomain)
      CALL Domain_MappingsFinalise(domain%mappings,err,error,*999)        
      CALL Domain_TopologyFinalise(domain%topology,err,error,*999)
      DEALLOCATE(domain)
    ENDIF
   
    EXITS("Domain_Finalise")
    RETURN
999 ERRORSEXITS("Domain_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE Domain_Finalise
  
  !
  !================================================================================================================================
  !

  !>Finalises the domains for a given decomposition and deallocates all memory.
  SUBROUTINE Decomposition_DomainsFinalise(decomposition,err,error,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer to the decomposition to finalise the domains for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx
    
    ENTERS("Decomposition_DomainsFinalise",err,error,*999)

    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)

    IF(ALLOCATED(decomposition%domain)) THEN
      DO componentIdx=1,SIZE(decomposition%domain,1)
        CALL Domain_Finalise(decomposition%domain(componentIdx)%ptr,err,error,*999)
      ENDDO !componentIdx
      DEALLOCATE(decomposition%domain)
    ENDIF
   
    EXITS("Decomposition_DomainsFinalise")
    RETURN
999 ERRORSEXITS("Decomposition_DomainsFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_DomainsFinalise
  
  !
  !================================================================================================================================
  !

  !>Allocates and initialises a domain.
  SUBROUTINE Domain_Initialise(domain,err,error,*)

    !Argument variables
    TYPE(DomainType), POINTER :: domain !<A pointer to the domain to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("Domain_Initialise",err,error,*998)

    IF(ASSOCIATED(domain)) CALL FlagError("Domain is already associated.",err,error,*998)
    
    ALLOCATE(domain,STAT=err)    
    IF(err/=0) CALL FlagError("Domain could not be allocated.",err,error,*999)
    NULLIFY(domain%decomposition)
    NULLIFY(domain%mesh)
    NULLIFY(domain%region)
    domain%meshComponentNumber=0
    domain%numberOfDimensions=0
    NULLIFY(domain%mappings)
    NULLIFY(domain%topology)
    CALL Domain_MappingsInitialise(domain,err,error,*999)
    CALL Domain_TopologyInitialise(domain,err,error,*999)
    
    EXITS("Domain_Initialise")
    RETURN
999 CALL Domain_Finalise(domain,dummyErr,dummyError,*998)
998 ERRORSEXITS("Domain_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE Domain_Initialise
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE Decomposition_DomainsInitialise(decomposition,err,error,*)

    !Argument variables
    TYPE(DecompositionType), POINTER :: decomposition !<A pointer to the decomposition to initialise the domains for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,dummyErr
    TYPE(MeshType), POINTER :: mesh
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("Decomposition_DomainsInitialise",err,error,*999)

    IF(.NOT.ASSOCIATED(decomposition)) CALL FlagError("Decomposition is not associated.",err,error,*999)
    IF(ALLOCATED(decomposition%domain)) CALL FlagError("Decomposition already has domain allocated.",err,error,*999)

    NULLIFY(mesh)
    CALL Decomposition_MeshGet(decomposition,mesh,err,error,*999)
    ALLOCATE(decomposition%domain(decomposition%numberOfComponents),STAT=err)
    IF(err/=0) CALL FlagError("Decomposition domain could not be allocated.",err,error,*999)
    DO componentIdx=1,decomposition%numberOfComponents
      NULLIFY(decomposition%domain(componentIdx)%ptr)
      CALL Domain_Initialise(decomposition%domain(componentIdx)%ptr,err,error,*999)
      decomposition%domain(componentIdx)%ptr%decomposition=>decomposition
      decomposition%domain(componentIdx)%ptr%mesh=>mesh
      decomposition%domain(componentIdx)%ptr%meshComponentNumber=componentIdx
      decomposition%domain(componentIdx)%ptr%region=>decomposition%region
      decomposition%domain(componentIdx)%ptr%numberOfDimensions=decomposition%numberOfDimensions
    ENDDO !componentIdx
    
    EXITS("Decomposition_DomainsInitialise")
    RETURN
999 CALL Decomposition_DomainsFinalise(decomposition,dummyErr,dummyError,*998)
998 ERRORSEXITS("Decomposition_DomainsInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE Decomposition_DomainsInitialise
  
  !
  !================================================================================================================================
  !

  !>Finalises the dofs mapping.
  SUBROUTINE DomainMappings_DofsFinalise(dofsMapping,err,error,*)

    !Argument variables
    TYPE(DomainMappingType), POINTER :: dofsMapping !<A pointer to the dofs mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("DomainMappings_DofsFinalise",err,error,*999)

    CALL DOMAIN_MAPPINGS_MAPPING_FINALISE(dofsMapping,err,error,*999)
 
    EXITS("DomainMappings_DofsFinalise")
    RETURN
999 ERRORSEXITS("DomainMappings_DofsFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE DomainMappings_DofsFinalise

  !
  !================================================================================================================================
  !

  !>Intialises the dofs mapping in the given domain mapping.
  SUBROUTINE DomainMappings_DofsInitialise(domainMappings,err,error,*)

    !Argument variables
    TYPE(DomainMappingsType), POINTER :: domainMappings !<A pointer to the domain mappings to initialise the dofs for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("DomainMappings_DofsInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(domainMappings)) CALL FlagError("Domain mappings is not associated.",err,error,*998)
    IF(ASSOCIATED(domainMappings%dofs)) CALL FlagError("Domain dofs mappings are already associated.",err,error,*998)
     
    ALLOCATE(domainMappings%dofs,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate domain mappings dofs.",err,error,*999)
    CALL DomainMappings_MappingInitialise(domainMappings%dofs,err,error,*999)
 
    EXITS("DomainMappings_DofsInitialise")
    RETURN
999 CALL DomainMappings_DofsFinalise(domainMappings%dofs,dummyErr,dummyError,*998)
998 ERRORSEXITS("DomainMappings_DofsInitialise",err,error)
    RETURN 1
   
  END SUBROUTINE DomainMappings_DofsInitialise

  !
  !================================================================================================================================
  !

  !>Calculates the local/global element mappings for a domain decomposition.
  SUBROUTINE DomainMappings_ElementsCalculate(domainMappings,err,error,*)

    !Argument variables
    TYPE(DomainMappingsType), POINTER :: domainMappings !<A pointer to the domain to calculate the element mappings for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,adjacentElementIdx,adjacentElement,domainNumber,domainIdx,elementIdx,localNodeIdx,nodeNumber, &
      & numberOfDomains,numberOfAdjacentElements,componentIdx
    INTEGER(INTG), ALLOCATABLE :: adjacentElements(:),domains(:),localElementNumbers(:)
    TYPE(BasisType), POINTER :: basis
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainMappingType), POINTER :: elementsMapping
    TYPE(LIST_TYPE), POINTER :: adjacentDomainsList
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: adjacentElementsList(:)
    TYPE(MeshType), POINTER :: mesh
    TYPE(MeshElementsType), POINTER :: meshElements
    TYPE(MeshNodesType), POINTER :: meshNodes
    TYPE(MeshTopologyType), POINTER :: meshTopology
    TYPE(VARYING_STRING) :: dummyError

    NULLIFY(elementsMapping)
    
    ENTERS("DomainMappings_ElementsCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(domainMappings)) CALL FlagError("Domain mappings is not associated.",err,error,*999)
    
    CALL DomainMappings_ElementsMappingGet(domainMappings,elementsMapping,err,error,*999)
    NULLIFY(domain)
    CALL DomainMappings_DomainGet(domainMappings,domain,err,error,*999)
    NULLIFY(decomposition)
    CALL Domain_DecompositionGet(domain,decomposition,err,error,*999)
    NULLIFY(mesh)
    CALL Decomposition_MeshGet(decomposition,mesh,err,error,*999)
    componentIdx=decomposition%meshComponentNumber
    NULLIFY(meshTopology)
    CALL Mesh_MeshTopologyGet(mesh,componentIdx,meshTopology,err,error,*999)
    NULLIFY(meshElements)
    CALL MeshTopology_MeshElementsGet(meshTopology,meshElements,err,error,*999)
    NULLIFY(meshNodes)
    CALL MeshTopology_MeshNodesGet(meshTopology,meshNodes,err,error,*999)
                            
    !Calculate the local and global numbers and set up the mappings                           
    ALLOCATE(elementsMapping%globalToLocalMap(mesh%numberOfElements),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate element mapping global to local map.",err,error,*999)
    elementsMapping%numberOfGlobal=meshElements%numberOfElements
    !Loop over the global elements and calculate local numbers
    ALLOCATE(localElementNumbers(0:decomposition%numberOfDomains-1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate local element numbers.",err,error,*999)
    localElementNumbers=0
    ALLOCATE(adjacentElementsList(0:decomposition%numberOfDomains-1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate adjacent elements list.",err,error,*999)
    DO domainIdx=0,decomposition%numberOfDomains-1
      NULLIFY(adjacentElementsList(domainIdx)%ptr)
      CALL List_CreateStart(adjacentElementsList(domainIdx)%ptr,err,error,*999)
      CALL List_DataTypeSet(adjacentElementsList(domainIdx)%ptr,LIST_INTG_TYPE,err,error,*999)
      CALL List_InitialSizeSet(adjacentElementsList(domainIdx)%ptr,MAX(INT(MESH%numberOfElements/2),1),err,error,*999)
      CALL List_CreateFinish(adjacentElementsList(domainIdx)%ptr,err,error,*999)
    ENDDO !domainIdx
    
    DO elementIdx=1,mesh%numberOfElements
      !Calculate the local numbers
      domainNumber=decomposition%elementDomain(elementIdx)
      localElementNumbers(domainNumber)=localElementNumbers(domainNumber)+1
      !Calculate the adjacent elements to the computation domains and the adjacent domain numbers themselves
      NULLIFY(basis)
      CALL MeshElements_BasisGet(meshElements,elementIdx,basis,err,error,*999)
      NULLIFY(adjacentDomainsList)
      CALL List_CreateStart(adjacentDomainsList,err,error,*999)
      CALL List_DataTypeSet(adjacentDomainsList,LIST_INTG_TYPE,err,error,*999)
      CALL List_InitialSizeSet(adjacentDomainsList,decomposition%numberOfDomains,err,error,*999)
      CALL List_CreateFinish(adjacentDomainsList,err,error,*999)
      CALL List_ItemAdd(adjacentDomainsList,domainNumber,err,error,*999)
      DO localNodeIdx=1,basis%numberOfNodes
        nodeNumber=meshElements%elements(elementIdx)%meshElementNodes(localNodeIdx)
        DO adjacentElementIdx=1,meshNodes%nodes(nodeNumber)%numberOfSurroundingElements
          adjacentElement=meshNodes%nodes(nodeNumber)%surroundingElements(adjacentElementIdx)
          IF(decomposition%elementDomain(adjacentElement)/=domainNumber) THEN
            CALL List_ItemAdd(adjacentElementsList(domainNumber)%ptr,adjacentElement,err,error,*999)
            CALL List_ItemAdd(adjacentDomainsList,decomposition%elementDomain(adjacentElement),err,error,*999)
          ENDIF
        ENDDO !adjacentElementIdx
      ENDDO !localNodeIdx
      CALL List_RemoveDuplicates(adjacentDomainsList,err,error,*999)
      CALL List_DetachAndDestroy(adjacentDomainsList,numberOfDomains,domains,err,error,*999)
      DEALLOCATE(domains)
      CALL DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE(elementsMapping%globalToLocalMap(elementIdx),err,error,*999)
      ALLOCATE(elementsMapping%globalToLocalMap(elementIdx)%localNumber(numberOfDomains),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate element global to local map local number.",err,error,*999)
      ALLOCATE(elementsMapping%globalToLocalMap(elementIdx)%domainNumber(numberOfDomains),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate element global to local map domain number.",err,error,*999)
      ALLOCATE(elementsMapping%globalToLocalMap(elementIdx)%localType(numberOfDomains),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate element global to local map local type.",err,error,*999)
      elementsMapping%globalToLocalMap(elementIdx)%numberOfDomains=1
      elementsMapping%globalToLocalMap(elementIdx)%localNumber(1)=localElementNumbers(domainNumber)
      elementsMapping%globalToLocalMap(elementIdx)%domainNumber(1)=decomposition%elementDomain(elementIdx) 
      IF(numberOfDomains==1) THEN
        !Element is an internal element
        elementsMapping%globalToLocalMap(elementIdx)%localType(1)=DOMAIN_LOCAL_INTERNAL
      ELSE
        !Element is on the boundary of computation domains
        elementsMapping%globalToLocalMap(elementIdx)%localType(1)=DOMAIN_LOCAL_BOUNDARY
      ENDIF
    ENDDO !elementIdx
    
    !Compute ghost element mappings
    DO domainIdx=0,decomposition%numberOfDomains-1
      CALL List_RemoveDuplicates(adjacentElementsList(domainIdx)%ptr,err,error,*999)
      CALL List_DetachAndDestroy(adjacentElementsList(domainIdx)%ptr,numberOfAdjacentElements,adjacentElements,err,error,*999)
      DO adjacentElementIdx=1,numberOfAdjacentElements
        adjacentElement=adjacentElements(adjacentElementIdx)
        localElementNumbers(domainIdx)=localElementNumbers(domainIdx)+1
        elementsMapping%globalToLocalMap(adjacentElement)%numberOfDomains= &
          & elementsMapping%globalToLocalMap(adjacentElement)%numberOfDomains+1
        elementsMapping%globalToLocalMap(adjacentElement)%localNumber(elementsMapping%globalToLocalMap(adjacentElement)% &
          & numberOfDomains)=localElementNumbers(domainIdx)
        elementsMapping%globalToLocalMap(adjacentElement)%domainNumber(elementsMapping%globalToLocalMap(adjacentElement)% &
          & numberOfDomains)=domainIdx
        elementsMapping%globalToLocalMap(adjacentElement)%localType(elementsMapping%globalToLocalMap(adjacentElement)% &
          & numberOfDomains)=DOMAIN_LOCAL_GHOST
      ENDDO !adjacentElementIdx
      IF(ALLOCATED(adjacentElements)) DEALLOCATE(adjacentElements)
    ENDDO !domainIdx
              
    DEALLOCATE(adjacentElementsList)
    DEALLOCATE(localElementNumbers)
    
    !Calculate element local to global maps from global to local map
    CALL DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE(elementsMapping,err,error,*999)
        
    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Element mappings :",err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Global to local map :",err,error,*999)
      DO elementIdx=1,MESH%numberOfElements
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Global element = ",elementIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of domains  = ", &
          & elementsMapping%globalToLocalMap(elementIdx)%numberOfDomains,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,elementsMapping%globalToLocalMap(elementIdx)% &
          & numberOfDomains,8,8,elementsMapping%globalToLocalMap(elementIdx)%localNumber, &
          & '("      Local number :",8(X,I7))','(20X,8(X,I7))',err,error,*999)      
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,elementsMapping%globalToLocalMap(elementIdx)% &
          & numberOfDomains,8,8,elementsMapping%globalToLocalMap(elementIdx)%domainNumber, &
          & '("      Domain number:",8(X,I7))','(20X,8(X,I7))',err,error,*999)      
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,elementsMapping%globalToLocalMap(elementIdx)% &
          & numberOfDomains,8,8,elementsMapping%globalToLocalMap(elementIdx)%localType, &
          & '("      Local type   :",8(X,I7))','(20X,8(X,I7))',err,error,*999)      
      ENDDO !elementIdx
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Local to global map :",err,error,*999)
      DO elementIdx=1,elementsMapping%totalNumberOfLocal
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Local element = ",elementIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Global element = ", &
          & elementsMapping%localToGlobalMap(elementIdx),err,error,*999)
      ENDDO !elementIdx
      IF(diagnostics2) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Internal elements :",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of internal elements = ", &
          & elementsMapping%numberOfInternal,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,elementsMapping%numberOfInternal,8,8, &
          & elementsMapping%domainList(elementsMapping%internalStart:elementsMapping%internalFinish), &
          & '("    Internal elements:",8(X,I7))','(22X,8(X,I7))',err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Boundary elements :",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of boundary elements = ", &
          & elementsMapping%numberOfBoundary,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,elementsMapping%numberOfBoundary,8,8, &
          & elementsMapping%domainList(elementsMapping%boundaryStart:elementsMapping%boundaryFinish), &
          & '("    Boundary elements:",8(X,I7))','(22X,8(X,I7))',err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Ghost elements :",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of ghost elements = ", &
          & elementsMapping%numberOfGhost,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,elementsMapping%numberOfGhost,8,8, &
          & elementsMapping%domainList(elementsMapping%ghostStart:elementsMapping%ghostFinish), &
          & '("    Ghost elements   :",8(X,I7))','(22X,8(X,I7))',err,error,*999)
      ENDIF
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Adjacent domains :",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of adjacent domains = ", &
        & elementsMapping%numberOfAdjacentDomains,err,error,*999)
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,elementsMapping%numberOfDomains+1,8,8, &
        & elementsMapping%adjacentDomainsPtr,'("    Adjacent domains ptr  :",8(X,I7))','(27X,8(X,I7))',err,error,*999)
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,elementsMapping%adjacentDomainsPtr( &
        & elementsMapping%numberOfDomains)-1,8,8,elementsMapping%adjacentDomainsList, &
        '("    Adjacent domains list :",8(X,I7))','(27X,8(X,I7))',err,error,*999)
      DO domainIdx=1,elementsMapping%numberOfAdjacentDomains
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Adjacent domain idx : ",domainIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Domain number = ", &
          & elementsMapping%adjacentDomains(domainIdx)%domainNumber,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of send ghosts    = ", &
          & elementsMapping%adjacentDomains(domainIdx)%numberOfSendGhosts,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,elementsMapping%adjacentDomains(domainIdx)% &
          & numberOfSendGhosts,6,6,elementsMapping%adjacentDomains(domainIdx)%localGhostSendIndices, &
          & '("      Local send ghost indicies       :",6(X,I7))','(39X,6(X,I7))',err,error,*999)      
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of recieve ghosts = ", &
          & elementsMapping%adjacentDomains(domainIdx)%numberOfReceiveGhosts,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,elementsMapping%adjacentDomains(domainIdx)% &
          & numberOfReceiveGhosts,6,6,elementsMapping%adjacentDomains(domainIdx)%localGhostReceiveIndices, &
          & '("      Local receive ghost indicies    :",6(X,I7))','(39X,6(X,I7))',err,error,*999)              
      ENDDO !domainIdx
    ENDIF
    
    EXITS("DomainMappings_ElementsCalculate")
    RETURN
999 IF(ALLOCATED(domains)) DEALLOCATE(domains)
    IF(ALLOCATED(adjacentElements)) DEALLOCATE(adjacentElements)    
    CALL DomainMappings_ElementsFinalise(elementsMapping,dummyErr,dummyError,*998)
998 ERRORSEXITS("DomainMappings_ElementsCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE DomainMappings_ElementsCalculate
  
  !
  !================================================================================================================================
  !

  !>Finalises the mappings in the given domain. 
  SUBROUTINE Domain_MappingsFinalise(domainMappings,err,error,*)

    !Argument variables
    TYPE(DomainMappingsType), POINTER :: domainMappings !<A pointer to the domain to finalise the mappings for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Domain_MappingsFinalise",err,error,*999)

    IF(ASSOCIATED(domainMappings)) THEN
      CALL DomainMappings_ElementsFinalise(domainMappings%elements,err,error,*999)
      CALL DomainMappings_NodesFinalise(domainMappings%nodes,err,error,*999)
      CALL DomainMappings_DofsFinalise(domainMappings%dofs,err,error,*999)
      DEALLOCATE(domainMappings)
    ENDIF
 
    EXITS("Domain_MappingsFinalise")
    RETURN
999 ERRORSEXITS("Domain_MappingsFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE Domain_MappingsFinalise

  !
  !================================================================================================================================
  !

  !>Finalises the element mapping in the given domain mapping. 
  SUBROUTINE DomainMappings_ElementsFinalise(domainElementsMapping,err,error,*)

    !Argument variables
    TYPE(DomainMappingType), POINTER :: domainElementsMapping !<A pointer to the domain mappings to finalise the elements for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("DomainMappings_ElementsFinalise",err,error,*999)

    CALL DomainMapping_Finalise(domainElementsMapping,err,error,*999)
 
    EXITS("DomainMappings_ElementsFinalise")
    RETURN
999 ERRORSEXITS("DomainMappings_ElementsFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE DomainMappings_ElementsFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the element mapping in the given domain mapping.
  SUBROUTINE DomainMappings_ElementsInitialise(domainMappings,err,error,*)

    !Argument variables
    TYPE(DomainMappingsType), POINTER :: domainMappings !<A pointer to the domain mappings to initialise the elements for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("DomainMappings_ElementsInitialise",err,error,*999)

    IF(.NOT.ASSOCIATED(domainMappings)) CALL FlagError("Domain mappings is not associated.",err,error,*999)
    IF(ASSOCIATED(domainMappings%elements)) CALL FlagError("Domain elements mappings are already associated.",err,error,*999)
    
    ALLOCATE(domainMappings%elements,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate domain mappings elements.",err,error,*999)
    CALL DomainMappings_MappingInitialise(domainMappings%elements,err,error,*999)
 
    EXITS("DomainMappings_ElementsInitialise")
    RETURN
999 CALL DomainMappings_ElementsFinalise(domainMappings%elements,dummyErr,dummyError,*998)
998 ERRORSEXITS("DomainMappings_ElementsInitialise",err,error)
    RETURN 1
   
  END SUBROUTINE DomainMappings_ElementsInitialise
  
  !
  !================================================================================================================================
  !

  !>Initialises the mappings for a domain decomposition. 
  SUBROUTINE Domain_MappingsInitialise(domain,err,error,*)

    !Argument variables
    TYPE(DomainType), POINTER :: domain !<A pointer to the domain to initialise the mappings for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("Domain_MappingsInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(domain)) CALL FlagError("Domain is not associated.",err,error,*998)    
    IF(ASSOCIATED(domain%mappings)) CALL FlagError("Domain already has mappings associated.",err,error,*998)
     
    ALLOCATE(domain%mappings,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate domain mappings.",err,error,*999)
    domain%mappings%domain=>domain
    NULLIFY(domain%mappings%elements)
    NULLIFY(domain%mappings%nodes)
    NULLIFY(domain%mappings%dofs)
    !Calculate the node and element mappings
    CALL DomainMappings_ElementsInitialise(domain%mappings,err,error,*999)
    CALL DomainMappings_NodesInitialise(domain%mappings,err,error,*999)
    CALL DomainMappings_DofsInitialise(domain%mappings,err,error,*999)
    CALL DomainMappings_ElementsCalculate(domain%mappings,err,error,*999)
    CALL DomainMappings_NodesDofsCalculate(domain%mappings,err,error,*999)
   
    EXITS("Domain_MappingsInitialise")
    RETURN
999 CALL Domain_MappingsFinalise(domain%mappings,dummyErr,dummyError,*998)
998 ERRORSEXITS("Domain_MappingsInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE Domain_MappingsInitialise
  
  !
  !================================================================================================================================
  !

  !>Calculates the local/global node and dof mappings for a domain decomposition.
  SUBROUTINE DomainMappings_NodesDofsCalculate(domainMappings,err,error,*)

    !Argument variables
    TYPE(DomainMappingsType), POINTER :: domainMappings !<A pointer to the domain to calculate the node dofs for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,adjacentElementIdx,computationNodeIdx,ghostNodeIdx,adjacentElement,ghostNode, &
      & numberOfNodesPerDomain,domainIdx,domainIdx2,domainNumber,nodeIdx,derivativeIdx,versionIdx,dofIdx,numberOfDomains, &
      & maxNumberDomains,numberOfGhostNodes,myGroupComputationNodeNumber,numberOfGroupComputationNodes,componentIdx
    INTEGER(INTG), ALLOCATABLE :: localNodeNumbers(:),localDofNumbers(:),nodeCount(:),numberInternalNodes(:), &
      & numberBoundaryNodes(:)
    INTEGER(INTG), ALLOCATABLE :: domains(:),allDomains(:),ghostNodes(:)
    LOGICAL :: boundaryDomain
    TYPE(LIST_TYPE), POINTER :: adjacentDomainsList,allAdjacentDomainsList
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: ghostNodesList(:)
    TYPE(MeshType), POINTER :: mesh
    TYPE(MeshTopologyType), POINTER :: meshTopology
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainMappingType), POINTER :: elementsMapping
    TYPE(DomainMappingType), POINTER :: nodesMapping
    TYPE(DomainMappingType), POINTER :: dofsMapping
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("DomainMappings_NodesDofsCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(domainMappings)) CALL FlagError("Domain mappings is not associated.",err,error,*999)
    
    NULLIFY(domain)
    CALL DomainMappings_DomainGet(domainMappings,domain,err,error,*999)
    NULLIFY(nodesMapping)
    CALL DomainMappings_NodesMappingGet(domainMappings,nodesMapping,err,error,*999)
    NULLIFY(elementsMapping)
    CALL DomainMappings_ElementsMappingGet(domainMappings,elementsMapping,err,error,*999)
    NULLIFY(dofsMapping)
    CALL DomainMappings_DofsMappingGet(domainMappings,dofsMapping,err,error,*999)
    NULLIFY(decomposition)
    CALL Domain_DecompositionGet(domain,decomposition,err,error,*999)
    NULLIFY(mesh)
    CALL Decomposition_MeshGet(decomposition,mesh,err,error,*999)
    componentIdx=decomposition%meshComponentNumber
    NULLIFY(meshTopology)
    CALL Mesh_MeshTopologyGet(mesh,componentIdx,meshTopology,err,error,*999)
 
    CALL WorkGroup_NumberOfGroupNodesGet(nodesMapping%workGroup,numberOfGroupComputationNodes,err,error,*999)
    CALL WorkGroup_GroupNodeNumberGet(nodesMapping%workGroup,myGroupComputationNodeNumber,err,error,*999)
                   
    !Calculate the local and global numbers and set up the mappings
    ALLOCATE(nodesMapping%globalToLocalMap(meshTopology%nodes%numberOfNodes),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate node mapping global to local map.",err,error,*999)
    nodesMapping%numberOfGlobal=meshTopology%nodes%numberOfNodes
    ALLOCATE(localNodeNumbers(0:decomposition%numberOfDomains-1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate local node numbers.",err,error,*999)
    localNodeNumbers=0
    ALLOCATE(dofsMapping%globalToLocalMap(meshTopology%dofs%numberOfDofs),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate dofs mapping global to local map.",err,error,*999)
    dofsMapping%numberOfGlobal=meshTopology%dofs%numberOfDofs
    ALLOCATE(localDofNumbers(0:decomposition%numberOfDomains-1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate local dof numbers.",err,error,*999)
    localDofNumbers=0
    ALLOCATE(ghostNodesList(0:decomposition%numberOfDomains-1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate ghost nodes list.",err,error,*999)
    DO domainIdx=0,decomposition%numberOfDomains-1
      NULLIFY(ghostNodesList(domainIdx)%ptr)
      CALL List_CreateStart(ghostNodesList(domainIdx)%ptr,err,error,*999)
      CALL List_DataTypeSet(ghostNodesList(domainIdx)%ptr,LIST_INTG_TYPE,err,error,*999)
      CALL List_InitialSizeSet(ghostNodesList(domainIdx)%ptr,INT(meshTopology%nodes%numberOfNodes/2),err,error,*999)
      CALL List_CreateFinish(ghostNodesList(domainIdx)%ptr,err,error,*999)
    ENDDO !domainIdx
    ALLOCATE(numberInternalNodes(0:decomposition%numberOfDomains-1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate number of internal nodes.",err,error,*999)
    numberInternalNodes=0
    ALLOCATE(numberBoundaryNodes(0:decomposition%numberOfDomains-1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate number of boundary nodes.",err,error,*999)
    numberBoundaryNodes=0
    
    !For the first pass just determine the internal and boundary nodes
    DO nodeIdx=1,meshTopology%nodes%numberOfNodes
      NULLIFY(adjacentDomainsList)
      CALL List_CreateStart(adjacentDomainsList,err,error,*999)
      CALL List_DataTypeSet(adjacentDomainsList,LIST_INTG_TYPE,err,error,*999)
      CALL List_InitialSizeSet(adjacentDomainsList,decomposition%numberOfDomains,err,error,*999)
      CALL List_CreateFinish(adjacentDomainsList,err,error,*999)
      NULLIFY(allAdjacentDomainsList)
      CALL List_CreateStart(allAdjacentDomainsList,err,error,*999)
      CALL List_DataTypeSet(allAdjacentDomainsList,LIST_INTG_TYPE,err,error,*999)
      CALL List_InitialSizeSet(allAdjacentDomainsList,decomposition%numberOfDomains,err,error,*999)
      CALL List_CreateFinish(allAdjacentDomainsList,err,error,*999)
      DO adjacentElementIdx=1,meshTopology%nodes%nodes(nodeIdx)%numberOfSurroundingElements
        adjacentElement=meshTopology%nodes%nodes(nodeIdx)%surroundingElements(adjacentElementIdx)
        domainNumber=decomposition%elementDomain(adjacentElement)
        CALL List_ItemAdd(adjacentDomainsList,domainNumber,err,error,*999)
        DO domainIdx=1,elementsMapping%globalToLocalMap(adjacentElement)%numberOfDomains
          CALL List_ItemAdd(allAdjacentDomainsList,elementsMapping%globalToLocalMap(adjacentElement)%domainNumber(domainIdx), &
            & err,error,*999)
        ENDDO !domainIdx
      ENDDO !adjacentElementIdx
      CALL List_RemoveDuplicates(adjacentDomainsList,err,error,*999)
      CALL List_DetachAndDestroy(adjacentDomainsList,numberOfDomains,domains,err,error,*999)
      CALL List_RemoveDuplicates(allAdjacentDomainsList,err,error,*999)
      CALL List_DetachAndDestroy(allAdjacentDomainsList,maxNumberDomains,allDomains,err,error,*999)
      IF(numberOfDomains/=maxNumberDomains) THEN !Ghost node
        DO domainIdx=1,maxNumberDomains
          domainNumber=allDomains(domainIdx)
          boundaryDomain=.FALSE.
          DO domainIdx2=1,numberOfDomains
            IF(domainNumber==domains(domainIdx2)) THEN
              boundaryDomain=.TRUE.
              EXIT
            ENDIF
          ENDDO !domainIdx2
          IF(.NOT.boundaryDomain) CALL List_ItemAdd(ghostNodesList(domainNumber)%ptr,nodeIdx,err,error,*999)
        ENDDO !domainIdx
      ENDIF
      ALLOCATE(nodesMapping%globalToLocalMap(nodeIdx)%localNumber(maxNumberDomains),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate node global to local map local number.",err,error,*999)
      ALLOCATE(nodesMapping%globalToLocalMap(nodeIdx)%domainNumber(maxNumberDomains),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate node global to local map domain number.",err,error,*999)
      ALLOCATE(nodesMapping%globalToLocalMap(nodeIdx)%localType(maxNumberDomains),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate node global to local map local type.",err,error,*999)
      DO derivativeIdx=1,meshTopology%nodes%nodes(nodeIdx)%numberOfDerivatives
        DO versionIdx=1,meshTopology%nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%numberOfVersions
          dofIdx=meshTopology%nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%dofIndex(versionIdx)
          ALLOCATE(dofsMapping%globalToLocalMap(dofIdx)%localNumber(maxNumberDomains),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate dof global to local map local number.",err,error,*999)
          ALLOCATE(dofsMapping%globalToLocalMap(dofIdx)%domainNumber(maxNumberDomains),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate dof global to local map domain number.",err,error,*999)
          ALLOCATE(dofsMapping%globalToLocalMap(dofIdx)%localType(maxNumberDomains),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate dof global to local map local type.",err,error,*999)
        ENDDO !versionIdx
      ENDDO !derivativeIdx
      IF(numberOfDomains==1) THEN
        !Node is an internal node
        domainNumber=domains(1)
        numberInternalNodes(domainNumber)=numberInternalNodes(domainNumber)+1
        !localNodeNumbers(domainNumber)=localNodeNumbers(domainNumber)+1
        nodesMapping%globalToLocalMap(nodeIdx)%numberOfDomains=1
        !nodesMapping%globalToLocalMap(nodeIdx)%localNumber(1)=localNodeNumbers(domains(1))
        nodesMapping%globalToLocalMap(nodeIdx)%localNumber(1)=-1
        nodesMapping%globalToLocalMap(nodeIdx)%domainNumber(1)=domains(1) 
        nodesMapping%globalToLocalMap(nodeIdx)%localType(1)=DOMAIN_LOCAL_INTERNAL
        DO derivativeIdx=1,meshTopology%nodes%nodes(nodeIdx)%numberOfDerivatives
          DO versionIdx=1,meshTopology%nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%numberOfVersions
            dofIdx=meshTopology%nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%dofIndex(versionIdx)
            dofsMapping%globalToLocalMap(dofIdx)%numberOfDomains=1
            dofsMapping%globalToLocalMap(dofIdx)%localNumber(1)=-1
            dofsMapping%globalToLocalMap(dofIdx)%domainNumber(1)=domainNumber
            dofsMapping%globalToLocalMap(dofIdx)%localType(1)=DOMAIN_LOCAL_INTERNAL
          ENDDO !versionIdx
        ENDDO !derivativeIdx
      ELSE
        !Node is on the boundary of computation domains
        nodesMapping%globalToLocalMap(nodeIdx)%numberOfDomains=numberOfDomains
        DO derivativeIdx=1,meshTopology%nodes%nodes(nodeIdx)%numberOfDerivatives
          DO versionIdx=1,meshTopology%nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%numberOfVersions
            dofIdx=meshTopology%nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%dofIndex(versionIdx)
            dofsMapping%globalToLocalMap(dofIdx)%numberOfDomains=numberOfDomains
          ENDDO !versionIdx
        ENDDO !derivativeIdx
        DO domainIdx=1,numberOfDomains
          domainNumber=domains(domainIdx)
          !localNodeNumbers(domainNumber)=localNodeNumbers(domainNumber)+1
          !nodesMapping%globalToLocalMap(nodeIdx)%localNumber(domainIdx)=localNodeNumbers(domainNumber)
          nodesMapping%globalToLocalMap(nodeIdx)%localNumber(domainIdx)=-1
          nodesMapping%globalToLocalMap(nodeIdx)%domainNumber(domainIdx)=domainNumber
          nodesMapping%globalToLocalMap(nodeIdx)%localType(domainIdx)=DOMAIN_LOCAL_BOUNDARY
          DO derivativeIdx=1,meshTopology%nodes%nodes(nodeIdx)%numberOfDerivatives
            DO versionIdx=1,meshTopology%nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%numberOfVersions
              dofIdx=meshTopology%nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%dofIndex(versionIdx)
              dofsMapping%globalToLocalMap(dofIdx)%localNumber(domainIdx)=-1
              dofsMapping%globalToLocalMap(dofIdx)%domainNumber(domainIdx)=domainNumber
              dofsMapping%globalToLocalMap(dofIdx)%localType(domainIdx)=DOMAIN_LOCAL_BOUNDARY
            ENDDO !versionIdx
          ENDDO !derivativeIdx
        ENDDO !domainIdx
      ENDIF
      DEALLOCATE(domains) 
      DEALLOCATE(allDomains)
    ENDDO !nodeIdx
    
    !For the second pass assign boundary nodes to one domain on the boundary and set local node numbers.
    numberOfNodesPerDomain=FLOOR(REAL(meshTopology%nodes%numberOfNodes,DP)/REAL(decomposition%numberOfDomains,DP))
    ALLOCATE(domain%nodeDomain(meshTopology%nodes%numberOfNodes),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate node domain",err,error,*999)
    domain%nodeDomain=-1
    DO nodeIdx=1,meshTopology%nodes%numberOfNodes
      IF(nodesMapping%globalToLocalMap(nodeIdx)%numberOfDomains==1) THEN !Internal node
        domainNumber=nodesMapping%globalToLocalMap(nodeIdx)%domainNumber(1)
        domain%nodeDomain(nodeIdx)=domainNumber
        localNodeNumbers(domainNumber)=localNodeNumbers(domainNumber)+1
        nodesMapping%globalToLocalMap(nodeIdx)%localNumber(1)=localNodeNumbers(domainNumber)
        DO derivativeIdx=1,meshTopology%nodes%nodes(nodeIdx)%numberOfDerivatives
          DO versionIdx=1,meshTopology%nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%numberOfVersions
            dofIdx=meshTopology%nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%dofIndex(versionIdx)
            localDofNumbers(domainNumber)=localDofNumbers(domainNumber)+1
            dofsMapping%globalToLocalMap(dofIdx)%localNumber(1)=localDofNumbers(domainNumber)
          ENDDO !versionIdx
        ENDDO !derivativeIdx
      ELSE !Boundary node
        numberOfDomains=nodesMapping%globalToLocalMap(nodeIdx)%numberOfDomains
        DO domainIdx=1,numberOfDomains
          domainNumber=nodesMapping%globalToLocalMap(nodeIdx)%domainNumber(domainIdx)
          IF(domain%nodeDomain(nodeIdx)<0) THEN
            IF((numberInternalNodes(domainNumber)+numberBoundaryNodes(domainNumber)<numberOfNodesPerDomain).OR. &
              & (domainIdx==nodesMapping%globalToLocalMap(nodeIdx)%numberOfDomains)) THEN
              !Allocate the node to this domain
              domain%nodeDomain(nodeIdx)=domainNumber
              numberBoundaryNodes(domainNumber)=numberBoundaryNodes(domainNumber)+1
              localNodeNumbers(domainNumber)=localNodeNumbers(domainNumber)+1
              !Reset the boundary information to be in the first domain index. The remaining domain indicies will
              !be overwritten when the ghost nodes are calculated below. 
              nodesMapping%globalToLocalMap(nodeIdx)%numberOfDomains=1
              nodesMapping%globalToLocalMap(nodeIdx)%localNumber(1)=localNodeNumbers(domainNumber)
              nodesMapping%globalToLocalMap(nodeIdx)%domainNumber(1)=domainNumber
              nodesMapping%globalToLocalMap(nodeIdx)%localType(1)=DOMAIN_LOCAL_BOUNDARY
              DO derivativeIdx=1,meshTopology%nodes%nodes(nodeIdx)%numberOfDerivatives
                DO versionIdx=1,meshTopology%nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%numberOfVersions
                  dofIdx=meshTopology%nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%dofIndex(versionIdx)
                  localDofNumbers(domainNumber)=localDofNumbers(domainNumber)+1
                  dofsMapping%globalToLocalMap(dofIdx)%numberOfDomains=1
                  dofsMapping%globalToLocalMap(dofIdx)%localNumber(1)=localDofNumbers(domainNumber)
                  dofsMapping%globalToLocalMap(dofIdx)%domainNumber(1)=domainNumber
                  dofsMapping%globalToLocalMap(dofIdx)%localType(1)=DOMAIN_LOCAL_BOUNDARY
                ENDDO !versionIdx
              ENDDO !derivativeIdx
            ELSE
              !The node as already been assigned to a domain so it must be a ghost node in this domain
              CALL List_ItemAdd(ghostNodesList(domainNumber)%ptr,nodeIdx,err,error,*999)
            ENDIF
          ELSE
            !The node as already been assigned to a domain so it must be a ghost node in this domain
            CALL List_ItemAdd(ghostNodesList(domainNumber)%ptr,nodeIdx,err,error,*999)
          ENDIF
        ENDDO !domainIdx
      ENDIF
    ENDDO !nodeIdx
    DEALLOCATE(numberInternalNodes)
                  
    !Calculate ghost node and dof mappings
    DO domainIdx=0,decomposition%numberOfDomains-1
      CALL List_RemoveDuplicates(ghostNodesList(domainIdx)%ptr,err,error,*999)
      CALL List_DetachAndDestroy(ghostNodesList(domainIdx)%ptr,numberOfGhostNodes,ghostNodes,err,error,*999)
      DO ghostNodeIdx=1,numberOfGhostNodes
        ghostNode=ghostNodes(ghostNodeIdx)
        localNodeNumbers(domainIdx)=localNodeNumbers(domainIdx)+1
        nodesMapping%globalToLocalMap(ghostNode)%numberOfDomains=nodesMapping%globalToLocalMap(ghostNode)%numberOfDomains+1
        nodesMapping%globalToLocalMap(ghostNode)%localNumber(nodesMapping%globalToLocalMap(ghostNode)%numberOfDomains)= &
          & localNodeNumbers(domainIdx)
        nodesMapping%globalToLocalMap(ghostNode)%domainNumber(nodesMapping%globalToLocalMap(ghostNode)%numberOfDomains)=domainIdx
        nodesMapping%globalToLocalMap(ghostNode)%localType(nodesMapping%globalToLocalMap(ghostNode)%numberOfDomains)= &
          & DOMAIN_LOCAL_GHOST
        DO derivativeIdx=1,meshTopology%nodes%nodes(ghostNode)%numberOfDerivatives
          DO versionIdx=1,meshTopology%nodes%nodes(ghostNode)%derivatives(derivativeIdx)%numberOfVersions
            dofIdx=meshTopology%nodes%nodes(ghostNode)%derivatives(derivativeIdx)%dofIndex(versionIdx)
            localDofNumbers(domainIdx)=localDofNumbers(domainIdx)+1
            dofsMapping%globalToLocalMap(dofIdx)%numberOfDomains=dofsMapping%globalToLocalMap(dofIdx)%numberOfDomains+1
            dofsMapping%globalToLocalMap(dofIdx)%localNumber(dofsMapping%globalToLocalMap(dofIdx)%numberOfDomains)= &
              & localDofNumbers(domainIdx)
            dofsMapping%globalToLocalMap(dofIdx)%domainNumber(dofsMapping%globalToLocalMap(dofIdx)%numberOfDomains)=domainIdx
            dofsMapping%globalToLocalMap(dofIdx)%localType(dofsMapping%globalToLocalMap(dofIdx)%numberOfDomains)= &
              & DOMAIN_LOCAL_GHOST
          ENDDO !versionIdx
        ENDDO !derivativeIdx
      ENDDO !ghostNodeIdx
      DEALLOCATE(ghostNodes)
    ENDDO !domainIdx
                  
    !Check decomposition and check that each domain has a node in it.
    ALLOCATE(nodeCount(0:numberOfGroupComputationNodes-1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate node count.",err,error,*999)
    nodeCount=0
    DO nodeIdx=1,meshTopology%nodes%numberOfNodes
      computationNodeIdx=domain%nodeDomain(nodeIdx)
      IF(computationNodeIdx>=0.AND.computationNodeIdx<numberOfGroupComputationNodes) THEN
        nodeCount(computationNodeIdx)=nodeCount(computationNodeIdx)+1
      ELSE
        localError="The computation node number of "// &
          & TRIM(NumberToVString(computationNodeIdx,"*",err,error))// &
          & " for node number "//TRIM(NumberToVString(nodeIdx,"*",err,error))// &
          & " is invalid. The computation node number must be between 0 and "// &
          & TRIM(NumberToVString(numberOfGroupComputationNodes-1,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDDO !nodeIdx
    DO computationNodeIdx=0,numberOfGroupComputationNodes-1
      IF(nodeCount(computationNodeIdx)==0) THEN
        localError="Invalid decomposition. There are no nodes in computation node "// &
          & TRIM(NumberToVString(computationNodeIdx,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDDO !computationNodeIdx
    DEALLOCATE(nodeCount)
    
    DEALLOCATE(ghostNodesList)
    DEALLOCATE(localNodeNumbers)
                  
    !Calculate node and dof local to global maps from global to local map
    CALL DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE(nodesMapping,err,error,*999)
    CALL DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE(dofsMapping,err,error,*999)
        
    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Node decomposition :",err,error,*999)
      DO nodeIdx=1,meshTopology%nodes%numberOfNodes
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Node = ",nodeIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Domain = ",domain%nodeDomain(nodeIdx),err,error,*999)
      ENDDO !nodeIdx
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Node mappings :",err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Global to local map :",err,error,*999)
      DO nodeIdx=1,meshTopology%nodes%numberOfNodes
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Global node = ",nodeIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of domains  = ", &
          & nodesMapping%globalToLocalMap(nodeIdx)%numberOfDomains,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,nodesMapping%globalToLocalMap(nodeIdx)% &
          & numberOfDomains,8,8,nodesMapping%globalToLocalMap(nodeIdx)%localNumber, &
          & '("      Local number :",8(X,I7))','(20X,8(X,I7))',err,error,*999)      
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,nodesMapping%globalToLocalMap(nodeIdx)% &
            & numberOfDomains,8,8,nodesMapping%globalToLocalMap(nodeIdx)%domainNumber, &
            & '("      Domain number:",8(X,I7))','(20X,8(X,I7))',err,error,*999)      
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,nodesMapping%globalToLocalMap(nodeIdx)% &
          & numberOfDomains,8,8,nodesMapping%globalToLocalMap(nodeIdx)%localType, &
          & '("      Local type   :",8(X,I7))','(20X,8(X,I7))',err,error,*999)      
      ENDDO !nodeIdx     
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Local to global map :",err,error,*999)
      DO nodeIdx=1,nodesMapping%totalNumberOfLocal
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Local node = ",nodeIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Global node = ", &
          & nodesMapping%localToGlobalMap(nodeIdx),err,error,*999)
      ENDDO !nodeIdx
      IF(diagnostics2) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Internal nodes :",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of internal nodes = ", &
          & nodesMapping%numberOfInternal,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,nodesMapping%numberOfInternal,8,8, &
          & nodesMapping%domainList(nodesMapping%internalStart:nodesMapping%internalFinish), &
          & '("    Internal nodes:",8(X,I7))','(19X,8(X,I7))',err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Boundary nodes :",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of boundary nodes = ", &
          & nodesMapping%numberOfBoundary,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,nodesMapping%numberOfBoundary,8,8, &
          & nodesMapping%domainList(nodesMapping%boundaryStart:nodesMapping%boundaryFinish), &
          & '("    Boundary nodes:",8(X,I7))','(19X,8(X,I7))',err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Ghost nodes :",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of ghost nodes = ", &
          & nodesMapping%numberOfGhost,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,nodesMapping%numberOfGhost,8,8, &
          & nodesMapping%domainList(nodesMapping%ghostStart:nodesMapping%ghostFinish), &
          & '("    Ghost nodes   :",8(X,I7))','(19X,8(X,I7))',err,error,*999)
      ENDIF
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Adjacent domains :",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of adjacent domains = ", &
        & nodesMapping%numberOfAdjacentDomains,err,error,*999)
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,nodesMapping%numberOfDomains+1,8,8, &
        & nodesMapping%adjacentDomainsPtr,'("    Adjacent domains ptr  :",8(X,I7))','(27X,8(X,I7))',err,error,*999)
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,nodesMapping%adjacentDomainsPtr( &
        & nodesMapping%numberOfDomains)-1,8,8,nodesMapping%adjacentDomainsList, &
        '("    Adjacent domains list :",8(X,I7))','(27X,8(X,I7))',err,error,*999)
      DO domainIdx=1,nodesMapping%numberOfAdjacentDomains
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Adjacent domain idx : ",domainIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Domain number = ", &
          & nodesMapping%adjacentDomains(domainIdx)%domainNumber,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of send ghosts    = ", &
          & nodesMapping%adjacentDomains(domainIdx)%numberOfSendGhosts,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,nodesMapping%adjacentDomains(domainIdx)% &
          & numberOfSendGhosts,6,6,nodesMapping%adjacentDomains(domainIdx)%localGhostSendIndices, &
          & '("      Local send ghost indicies       :",6(X,I7))','(39X,6(X,I7))',err,error,*999)      
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of recieve ghosts = ", &
          & nodesMapping%adjacentDomains(domainIdx)%numberOfReceiveGhosts,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,nodesMapping%adjacentDomains(domainIdx)% &
          & numberOfReceiveGhosts,6,6,nodesMapping%adjacentDomains(domainIdx)%localGhostReceiveIndices, &
          & '("      Local receive ghost indicies    :",6(X,I7))','(39X,6(X,I7))',err,error,*999)              
      ENDDO !domainIdx
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Dofs mappings :",err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Global to local map :",err,error,*999)
      DO dofIdx=1,meshTopology%dofs%numberOfDofs
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Global dof = ",dofIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of domains  = ", &
          & dofsMapping%globalToLocalMap(dofIdx)%numberOfDomains,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,dofsMapping%globalToLocalMap(dofIdx)% &
          & numberOfDomains,8,8,dofsMapping%globalToLocalMap(dofIdx)%localNumber, &
          & '("      Local number :",8(X,I7))','(20X,8(X,I7))',err,error,*999)      
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,dofsMapping%globalToLocalMap(dofIdx)% &
            & numberOfDomains,8,8,dofsMapping%globalToLocalMap(dofIdx)%domainNumber, &
            & '("      Domain number:",8(X,I7))','(20X,8(X,I7))',err,error,*999)      
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,dofsMapping%globalToLocalMap(dofIdx)% &
          & numberOfDomains,8,8,dofsMapping%globalToLocalMap(dofIdx)%localType, &
          & '("      Local type   :",8(X,I7))','(20X,8(X,I7))',err,error,*999)      
      ENDDO !dofIdx   
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Local to global map :",err,error,*999)
      DO dofIdx=1,dofsMapping%totalNumberOfLocal
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Local dof = ",dofIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Global dof = ", &
          & dofsMapping%localToGlobalMap(dofIdx),err,error,*999)
      ENDDO !nodeIdx
      IF(DIAGNOSTICS2) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Internal dofs :",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of internal dofs = ", &
          & dofsMapping%numberOfInternal,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,dofsMapping%numberOfInternal,8,8, &
          & dofsMapping%domainList(dofsMapping%internalStart:dofsMapping%internalFinish), &
          & '("    Internal dofs:",8(X,I7))','(18X,8(X,I7))',err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Boundary dofs :",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of boundary dofs = ", &
          & dofsMapping%numberOfBoundary,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,dofsMapping%numberOfBoundary,8,8, &
          & dofsMapping%domainList(dofsMapping%boundaryStart:dofsMapping%boundaryFinish), &
          & '("    Boundary dofs:",8(X,I7))','(18X,8(X,I7))',err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Ghost dofs :",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of ghost dofs = ", &
          & dofsMapping%numberOfGhost,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,dofsMapping%numberOfGhost,8,8, &
          & dofsMapping%domainList(dofsMapping%ghostStart:dofsMapping%ghostFinish), &
          & '("    Ghost dofs   :",8(X,I7))','(18X,8(X,I7))',err,error,*999)
      ENDIF
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Adjacent domains :",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of adjacent domains = ", &
        & dofsMapping%numberOfAdjacentDomains,err,error,*999)
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,dofsMapping%numberOfDomains+1,8,8, &
        & dofsMapping%adjacentDomainsPtr,'("    Adjacent domains ptr  :",8(X,I7))','(27X,8(X,I7))',err,error,*999)
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,dofsMapping%adjacentDomainsPtr( &
        & dofsMapping%numberOfDomains)-1,8,8,dofsMapping%adjacentDomainsList, &
        '("    Adjacent domains list :",8(X,I7))','(27X,8(X,I7))',err,error,*999)
      DO domainIdx=1,dofsMapping%numberOfAdjacentDomains
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Adjacent domain idx : ",domainIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Domain number = ", &
          & dofsMapping%adjacentDomains(domainIdx)%domainNumber,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of send ghosts    = ", &
          & dofsMapping%adjacentDomains(domainIdx)%numberOfSendGhosts,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,dofsMapping%adjacentDomains(domainIdx)% &
          & numberOfSendGhosts,6,6,dofsMapping%adjacentDomains(domainIdx)%localGhostSendIndices, &
          & '("      Local send ghost indicies       :",6(X,I7))','(39X,6(X,I7))',err,error,*999)      
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of recieve ghosts = ", &
          & dofsMapping%adjacentDomains(domainIdx)%numberOfReceiveGhosts,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,dofsMapping%adjacentDomains(domainIdx)% &
          & numberOfReceiveGhosts,6,6,dofsMapping%adjacentDomains(domainIdx)%localGhostReceiveIndices, &
          & '("      Local receive ghost indicies    :",6(X,I7))','(39X,6(X,I7))',err,error,*999)              
      ENDDO !domainIdx
    ENDIF
    
    EXITS("DomainMappings_NodesDofsCalculate")
    RETURN
999 IF(ALLOCATED(domains)) DEALLOCATE(domains)
    IF(ALLOCATED(allDomains)) DEALLOCATE(allDomains)
    IF(ALLOCATED(ghostNodes)) DEALLOCATE(ghostNodes)
    IF(ALLOCATED(numberInternalNodes)) DEALLOCATE(numberInternalNodes)
    IF(ALLOCATED(numberBoundaryNodes)) DEALLOCATE(numberBoundaryNodes)
    IF(ASSOCIATED(nodesMapping)) CALL DomainMappings_NodesFinalise(nodesMapping,dummyErr,dummyError,*998)
998 IF(ASSOCIATED(dofsMapping)) CALL DomainMappings_DofsFinalise(dofsMapping,dummyErr,dummyError,*997)
997 ERRORSEXITS("DomainMappings_NodesDofsCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE DomainMappings_NodesDofsCalculate
  
  !
  !================================================================================================================================
  !

  !>Finalises the node mapping in the given domain mappings.
  SUBROUTINE DomainMappings_NodesFinalise(domainNodesMapping,err,error,*)

    !Argument variables
    TYPE(DomainMappingType), POINTER :: domainNodesMapping !<A pointer to the nodes domain mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("DomainMappings_NodesFinalise",err,error,*999)

    CALL DomainMapping_Finalise(domainNodesMapping,err,error,*999)
 
    EXITS("DomainMappings_NodesFinalise")
    RETURN
999 ERRORSEXITS("DomainMappings_NodesFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE DomainMappings_NodesFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the node mapping in the given domain mapping.
  SUBROUTINE DomainMappings_NodesInitialise(domainMappings,err,error,*)

    !Argument variables
    TYPE(DomainMappingsType), POINTER :: domainMappings !<A pointer to the domain mappings to initialise the nodes for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("DomainMappings_NodesInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(domainMappings)) CALL FlagError("Domain mappings is not associated.",err,error,*998)
    IF(ASSOCIATED(domainMappings%nodes)) &
        & CALL FlagError("Domain nodes mappings are already associated.",err,error,*998)
      
    ALLOCATE(domainMappings%nodes,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate domain mappings nodes.",err,error,*999)
    CALL DomainMappings_MappingInitialise(domainMappings%nodes,err,error,*999)
 
    EXITS("DomainMappings_NodesInitialise")
    RETURN
999 CALL DomainMappings_NodesFinalise(domainMappings%nodes,dummyErr,dummyError,*998)
998 ERRORSEXITS("DomainMappings_NodesInitialise",err,error)
    RETURN 1
   
  END SUBROUTINE DomainMappings_NodesInitialise

  !
  !================================================================================================================================
  !

  !>Calculates the domain topology.
  SUBROUTINE DomainTopology_Calculate(domainTopology,err,error,*)

   !Argument variables
    TYPE(DomainTopologyType), POINTER :: domainTopology !<A pointer to the domain topology to calculate.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementIdx,nodeIdx
    TYPE(BasisType), POINTER :: basis
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainNodesType), POINTER :: domainNodes

    ENTERS("DomainTopology_Calculate",err,error,*999)

    IF(.NOT.ASSOCIATED(domainTopology)) CALL FlagError("Domain topology is not associated.",err,error,*999)

    NULLIFY(domainElements)
    CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
    NULLIFY(domainNodes)
    CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
    
    !Find maximum number of element parameters for all elements in the domain toplogy
    domainElements%maximumNumberOfElementParameters=-1
    DO elementIdx=1,domainElements%totalNumberOfElements
      NULLIFY(basis)
      CALL DomainElements_BasisGet(domainElements,elementIdx,basis,err,error,*999)
      IF(basis%numberOfElementParameters>domainElements%maximumNumberOfElementParameters) &
        & domainElements%maximumNumberOfElementParameters=basis%numberOfElementParameters
    ENDDO !elementIdx
    !Find maximum number of derivatives for all nodes in the domain toplogy
    domainNodes%maximumNumberOfDerivatives=-1
    DO nodeIdx=1,domainNodes%totalNumberOfNodes
      IF(domainNodes%nodes(nodeIdx)%numberOfDerivatives>domainNodes%maximumNumberOfDerivatives) &
        & domainNodes%maximumNumberOfDerivatives=domainNodes%nodes(nodeIdx)%numberOfDerivatives
    ENDDO !nodeIdx
    !Calculate the elements surrounding the nodes in the domain topology
    CALL DomainTopology_NodesSurroundingElementsCalculate(domainTopology,err,error,*999)
    
    EXITS("DomainTopology_Calculate")
    RETURN
999 ERRORSEXITS("DomainTopology_Calculate",err,error)
    RETURN 1
    
  END SUBROUTINE DomainTopology_Calculate
  
  !
  !================================================================================================================================
  !

  !>Initialises the local domain topology from the mesh topology.
  SUBROUTINE DomainTopology_InitialiseFromMesh(domainTopology,err,error,*)

    !Argument variables
    TYPE(DomainTopologyType), POINTER :: domainTopology !<A pointer to the domain topology to initialise from the mesh topology for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,derivativeIdx,derivativeIdx2,dofIdx,elementIdx,globalElement,globalNode,insertStatus, &
      & localNode,localNodeIdx,nodeIdx,versionIdx
    LOGICAL :: found
    TYPE(BasisType), POINTER :: basis
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainDOFsType), POINTER :: domainDofs
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainMappingType), POINTER :: dofsMapping
    TYPE(DomainMappingType), POINTER :: elementsMapping
    TYPE(DomainMappingType), POINTER :: nodesMapping
    TYPE(DomainMappingsType), POINTER :: domainMappings
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(MeshType), POINTER :: mesh
    TYPE(MeshElementsType), POINTER :: meshElements
    TYPE(MeshNodesType), POINTER :: meshNodes
    TYPE(MeshTopologyType), POINTER :: meshTopology
 
    ENTERS("DomainTopology_InitialiseFromMesh",err,error,*999)
    
    IF(.NOT.ASSOCIATED(domainTopology)) CALL FlagError("Domain topology is not associated.",err,error,*999)

    NULLIFY(domain)
    CALL DomainTopology_DomainGet(domainTopology,domain,err,error,*999)
    NULLIFY(domainMappings)
    CALL Domain_DomainMappingsGet(domain,domainMappings,err,error,*999)
    NULLIFY(decomposition)
    CALL Domain_DecompositionGet(domain,decomposition,err,error,*999)
    NULLIFY(mesh)
    CALL Decomposition_MeshGet(decomposition,mesh,err,error,*999)
    componentIdx=decomposition%meshComponentNumber
    NULLIFY(meshTopology)
    CALL Mesh_MeshTopologyGet(mesh,componentIdx,meshTopology,err,error,*999)
    NULLIFY(meshElements)
    CALL MeshTopology_MeshElementsGet(meshTopology,meshElements,err,error,*999)
    NULLIFY(domainElements)
    CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
    NULLIFY(meshNodes)
    CALL MeshTopology_MeshNodesGet(meshTopology,meshNodes,err,error,*999)
    NULLIFY(domainNodes)
    CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
    NULLIFY(domainDofs)
    CALL DomainTopology_DomainDofsGet(domainTopology,domainDofs,err,error,*999)
    NULLIFY(elementsMapping)
    CALL DomainMappings_ElementsMappingGet(domainMappings,elementsMapping,err,error,*999)
    NULLIFY(nodesMapping)
    CALL DomainMappings_NodesMappingGet(domainMappings,nodesMapping,err,error,*999)
    NULLIFY(dofsMapping)
    CALL DomainMappings_DofsMappingGet(domainMappings,dofsMapping,err,error,*999)
    
    ALLOCATE(domainElements%elements(elementsMapping%totalNumberOfLocal),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate domain elements elements.",err,error,*999)
    domainElements%numberOfElements=elementsMapping%numberOfLocal
    domainElements%totalNumberOfElements=elementsMapping%totalNumberOfLocal
    domainElements%numberOfGlobalElements=elementsMapping%numberOfGlobal
    ALLOCATE(domainNodes%nodes(nodesMapping%totalNumberOfLocal),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate domain nodes nodes.",err,error,*999)
    domainNodes%numberOfNodes=nodesMapping%numberOfLocal
    domainNodes%totalNumberOfNodes=nodesMapping%totalNumberOfLocal
    domainNodes%numberOfGlobalNodes=nodesMapping%numberOfGlobal
    ALLOCATE(domainDofs%dofIndex(3,dofsMapping%totalNumberOfLocal),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate domain dofs dof index.",err,error,*999)
    domainDofs%numberOfDofs=dofsMapping%numberOfLocal
    domainDofs%totalNumberOfDofs=dofsMapping%totalNumberOfLocal
    domainDofs%numberOfGlobalDofs=dofsMapping%numberOfGlobal
    !Loop over the domain nodes and calculate the parameters from the mesh nodes
    CALL Tree_CreateStart(domainNodes%nodesTree,err,error,*999)
    CALL Tree_InsertTypeSet(domainNodes%nodesTree,TREE_NO_DUPLICATES_ALLOWED,err,error,*999)
    CALL Tree_CreateFinish(domainNodes%nodesTree,err,error,*999)
    dofIdx=0
    DO localNodeIdx=1,domainNodes%totalNumberOfNodes
      CALL DomainNode_Initialise(domainNodes%nodes(localNodeIdx),err,error,*999)
      globalNode=nodesMapping%localToGlobalMap(localNodeIdx)
      domainNodes%nodes(localNodeIdx)%localNumber=localNodeIdx
      domainNodes%nodes(localNodeIdx)%meshNumber=globalNode
      domainNodes%nodes(localNodeIdx)%globalNumber=meshNodes%nodes(globalNode)%globalNumber
      domainNodes%nodes(localNodeIdx)%userNumber=meshNodes%nodes(globalNode)%userNumber
      CALL Tree_ItemInsert(domainNodes%nodesTree,domainNodes%nodes(localNodeIdx)%userNumber,localNodeIdx,insertStatus, &
        & err,error,*999)
      domainNodes%nodes(localNodeIdx)%numberOfSurroundingElements=0
      domainNodes%nodes(localNodeIdx)%numberOfDerivatives=meshNodes%nodes(globalNode)%numberOfDerivatives
      ALLOCATE(domainNodes%nodes(localNodeIdx)%derivatives(meshNodes%nodes(globalNode)%numberOfDerivatives),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate domain node derivatives.",err,error,*999)
      DO derivativeIdx=1,domainNodes%nodes(localNodeIdx)%numberOfDerivatives
        CALL DomainNodeDerivative_Initialise(domainNodes%nodes(localNodeIdx)%derivatives(derivativeIdx),err,error,*999)
        domainNodes%nodes(localNodeIdx)%derivatives(derivativeIdx)%globalDerivativeIndex= & 
          & meshNodes%nodes(globalNode)%derivatives(derivativeIdx)%globalDerivativeIndex
        domainNodes%nodes(localNodeIdx)%derivatives(derivativeIdx)%partialDerivativeIndex= &
          & meshNodes%nodes(globalNode)%derivatives(derivativeIdx)%partialDerivativeIndex
        domainNodes%nodes(localNodeIdx)%derivatives(derivativeIdx)%numberOfVersions= & 
          & meshNodes%nodes(globalNode)%derivatives(derivativeIdx)%numberOfVersions
        ALLOCATE(domainNodes%nodes(localNodeIdx)%derivatives(derivativeIdx)%userVersionNumbers( &
          & meshNodes%nodes(globalNode)%derivatives(derivativeIdx)%numberOfVersions),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate node derivative version numbers.",err,error,*999)
        domainNodes%nodes(localNodeIdx)%derivatives(derivativeIdx)%userVersionNumbers(1: &
          & meshNodes%nodes(globalNode)%derivatives(derivativeIdx)%numberOfVersions)= &
          & meshNodes%nodes(globalNode)%derivatives(derivativeIdx)%userVersionNumbers(1: &
          & meshNodes%nodes(globalNode)%derivatives(derivativeIdx)%numberOfVersions)
        ALLOCATE(domainNodes%nodes(localNodeIdx)%derivatives(derivativeIdx)%dofIndex( &
          & meshNodes%nodes(globalNode)%derivatives(derivativeIdx)%numberOfVersions),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate node dervative versions dof index.",err,error,*999)
        DO versionIdx=1,domainNodes%nodes(localNodeIdx)%derivatives(derivativeIdx)%numberOfVersions
          dofIdx=dofIdx+1
          domainNodes%nodes(localNodeIdx)%derivatives(derivativeIdx)%dofIndex(versionIdx)=dofIdx
          domainDofs%dofIndex(1,dofIdx)=versionIdx
          domainDofs%dofIndex(2,dofIdx)=derivativeIdx
          domainDofs%dofIndex(3,dofIdx)=localNodeIdx
        ENDDO !versionIdx
      ENDDO !derivativeIdx
      domainNodes%nodes(localNodeIdx)%boundaryNode=meshNodes%nodes(globalNode)%boundaryNode
    ENDDO !localNodeIdx
    !Loop over the domain elements and renumber from the mesh elements
    DO elementIdx=1,domainElements%totalNumberOfElements
      CALL DomainElement_Initialise(domainElements%elements(elementIdx),err,error,*999)
      globalElement=elementsMapping%localToGlobalMap(elementIdx)
      NULLIFY(basis)
      CALL MeshElements_BasisGet(meshElements,globalElement,basis,err,error,*999)
      domainElements%elements(elementIdx)%number=elementIdx
      domainElements%elements(elementIdx)%basis=>basis
      ALLOCATE(domainElements%elements(elementIdx)%elementNodes(basis%numberOfNodes),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate domain elements element nodes.",err,error,*999)
      ALLOCATE(domainElements%elements(elementIdx)%elementDerivatives(basis%maximumNumberOfDerivatives,basis%numberOfNodes), &
        & STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate domain elements element derivatives.",err,error,*999)
      ALLOCATE(domainElements%elements(elementIdx)%elementVersions(basis%maximumNumberOfDerivatives,basis%numberOfNodes),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate domain elements element versions.",err,error,*999)
      DO localNodeIdx=1,basis%numberOfNodes
        globalNode=meshElements%elements(globalElement)%meshElementNodes(localNodeIdx)
        localNode=nodesMapping%globalToLocalMap(globalNode)%localNumber(1)
        domainElements%elements(elementIdx)%elementNodes(localNodeIdx)=localNode
        DO derivativeIdx=1,basis%numberOfDerivatives(localNodeIdx)
          !Find equivalent node derivative by matching partial derivative index
          !/todo Take a look at this - is it needed any more?
          found=.FALSE.
          DO derivativeIdx2=1,domainNodes%nodes(localNode)%numberOfDerivatives
            IF(domainNodes%nodes(localNode)%derivatives(derivativeIdx2)%partialDerivativeIndex == &
              & basis%partialDerivativeIndex(derivativeIdx,localNodeIdx)) THEN
              found=.TRUE.
              EXIT
            ENDIF
          ENDDO !derivativeIdx2
          IF(found) THEN
            domainElements%elements(elementIdx)%elementDerivatives(derivativeIdx,localNodeIdx)=derivativeIdx2
            domainElements%elements(elementIdx)%elementVersions(derivativeIdx,localNodeIdx) = & 
              & meshElements%elements(globalElement)%userElementNodeVersions(derivativeIdx,localNodeIdx)
          ELSE
            CALL FlagError("Could not find equivalent node derivative",err,error,*999)
          ENDIF
        ENDDO !derivativeIdx
      ENDDO !localNodeIdx
    ENDDO !elementIdx                       
   
    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Initialised domain topology :",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Total number of domain nodes = ",domainNodes%totalNumberOfNodes, &
        & err,error,*999)
      DO nodeIdx=1,domainNodes%totalNumberOfNodes
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Node number = ",domainNodes%nodes(nodeIdx)%localNumber, &
          & err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Node mesh number = ",domainNodes%nodes(nodeIdx)%meshNumber, &
          & err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Node global number = ",domainNodes%nodes(nodeIdx)%globalNumber, &
          & err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Node user number = ",domainNodes%nodes(nodeIdx)%userNumber, &
          & err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Number of derivatives = ", &
          & domainNodes%nodes(nodeIdx)%numberOfDerivatives,err,error,*999)
        DO derivativeIdx=1,domainNodes%nodes(nodeIdx)%numberOfDerivatives
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Node local derivative number = ",derivativeIdx,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Global derivative index = ", &
            & domainNodes%nodes(nodeIdx)%derivatives(derivativeIdx)%globalDerivativeIndex,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Partial derivative index = ", &
            & domainNodes%nodes(nodeIdx)%derivatives(derivativeIdx)%partialDerivativeIndex,err,error,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1, &
            & domainNodes%nodes(nodeIdx)%derivatives(derivativeIdx)%numberOfVersions,4,4, &
            & domainNodes%nodes(nodeIdx)%derivatives(derivativeIdx)%dofIndex, &
            & '("        Degree-of-freedom index(versionIdx)  :",4(X,I9))','(36X,4(X,I9))',err,error,*999)
        ENDDO !derivativeIdx
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Boundary node = ", &
          & domainNodes%nodes(nodeIdx)%boundaryNode,err,error,*999)
     ENDDO !nodeIdx
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Total Number of domain dofs = ",domainDofs%totalNumberOfDofs, &
        & err,error,*999)
      DO dofIdx=1,domainDofs%totalNumberOfDofs
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Dof number = ",dofIdx,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,3,3, &
          & domainDofs%dofIndex(:,dofIdx),'("    Degree-of-freedom index :",3(X,I9))','(29X,3(X,I9))', &
          & err,error,*999)
      ENDDO !dofIdx
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Total number of domain elements = ", &
        & domainElements%totalNumberOfElements,err,error,*999)
      DO elementIdx=1,domainElements%totalNumberOfElements
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Element number = ",domainElements%elements(elementIdx)%number, &
          & err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Basis user number = ", &
          & domainElements%elements(elementIdx)%basis%userNumber,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of local nodes = ", &
          & domainElements%elements(elementIdx)%basis%numberOfNodes,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,domainElements%elements(elementIdx)%basis%numberOfNodes,8,8, &
          & domainElements%elements(elementIdx)%elementNodes,'("    Element nodes(nn) :",8(X,I9))','(23X,8(X,I9))', &
          & err,error,*999)
        DO localNodeIdx=1,domainElements%elements(elementIdx)%basis%numberOfNodes
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Local node number : ",localNodeIdx,err,error,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,domainElements%elements(elementIdx)%basis% &
            & numberOfDerivatives(localNodeIdx),8,8,domainElements%elements(elementIdx)%elementDerivatives(:,localNodeIdx), &
            & '("        Element derivatives :",8(X,I2))','(29X,8(X,I2))',err,error,*999)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,domainElements%elements(elementIdx)%basis% &
            & numberOfDerivatives(localNodeIdx),8,8,domainElements%elements(elementIdx)%elementVersions(:,localNodeIdx), &
            & '("        Element versions    :",8(X,I2))','(29X,8(X,I2))',err,error,*999)
        ENDDO !localNodeIdx
      ENDDO !elementIdx
    ENDIF
    
    EXITS("DomainTopology_InitialiseFromMesh")
    RETURN
999 ERRORSEXITS("DomainTopology_InitialiseFromMesh",err,error)
    RETURN 1
    
  END SUBROUTINE DomainTopology_InitialiseFromMesh
  
  !
  !================================================================================================================================
  !

  !>Finalises the dofs in the given domain topology. 
  SUBROUTINE DomainTopology_DOFsFinalise(domainDofs,err,error,*)

    !Argument variables
    TYPE(DomainDOFsType), POINTER :: domainDofs !<A pointer to the domain dofs to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DomainTopology_DOFsFinalise",err,error,*999)

    IF(ASSOCIATED(domainDofs)) THEN
      IF(ALLOCATED(domainDofs%dofIndex)) DEALLOCATE(domainDofs%dofIndex)
      DEALLOCATE(domainDofs)
    ENDIF

    EXITS("DomainTopology_DOFsFinalise")
    RETURN
999 ERRORSEXITS("DomainTopology_DOFsFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE DomainTopology_DOFsFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the dofs data structures for a domain topology. \todo finalise on exit
  SUBROUTINE DomainTopology_DOFsInitialise(domainTopology,err,error,*)

    !Argument variables
    TYPE(DomainTopologyType), POINTER :: domainTopology !<A pointer to the domain topology to initialise the dofs for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("DomainTopology_DOFsInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(domainTopology)) CALL FlagError("Domain topology is not associated.",err,error,*998)
    IF(ASSOCIATED(domainTopology%dofs)) CALL FlagError("Domain topology already has dofs associated.",err,error,*998)
    
    ALLOCATE(domainTopology%dofs,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate domain topology dofs.",err,error,*999)
    domainTopology%dofs%numberOfDofs=0
    domainTopology%dofs%totalNumberOfDofs=0
    domainTopology%dofs%numberOfGlobalDofs=0
    domainTopology%dofs%domainTopology=>domainTopology
    
    EXITS("DomainTopology_DOFsInitialise")
    RETURN
999 CALL DomainTopology_DOFsFinalise(domainTopology%dofs,dummyErr,dummyError,*998)
998 ERRORSEXITS("DomainTopology_DOFsInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE DomainTopology_DOFsInitialise
  
  !
  !================================================================================================================================
  !

  !>Finalises the given domain element.
  SUBROUTINE DomainElement_Finalise(domainElement,err,error,*)

    !Argument variables
    TYPE(DomainElementType) :: domainElement !<The domain element to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DomainElement_Finalise",err,error,*999)

    IF(ALLOCATED(domainElement%elementNodes)) DEALLOCATE(domainElement%elementNodes)
    IF(ALLOCATED(domainElement%elementDerivatives)) DEALLOCATE(domainElement%elementDerivatives)
    IF(ALLOCATED(domainElement%elementVersions)) DEALLOCATE(domainElement%elementVersions)
 
    EXITS("DomainElement_Finalise")
    RETURN
999 ERRORSEXITS("DomainElement_Finalise",err,error)
    RETURN 1
   
  END SUBROUTINE DomainElement_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the given domain element.
  SUBROUTINE DomainElement_Initialise(domainElement,err,error,*)

    !Argument variables
    TYPE(DomainElementType) :: domainElement !<The domain element to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DomainElement_Initialise",err,error,*999)

    domainElement%number=0
    NULLIFY(domainElement%basis)
  
    EXITS("DomainElement_Initialise")
    RETURN
999 ERRORSEXITS("DomainElement_Initialise",err,error)
    RETURN 1
   
  END SUBROUTINE DomainElement_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises the elements in the given domain topology. 
  SUBROUTINE DomainTopology_ElementsFinalise(domainElements,err,error,*)

    !Argument variables
    TYPE(DomainElementsType), POINTER :: domainElements !<A pointer to the domain topology to finalise the elements for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementIdx

    ENTERS("DomainTopology_ElementsFinalise",err,error,*999)

    IF(ASSOCIATED(domainElements)) THEN
      DO elementIdx=1,domainElements%totalNumberOfElements
        CALL DomainElement_Finalise(domainElements%elements(elementIdx),err,error,*999)
      ENDDO !elementIdx
      IF(ALLOCATED(domainElements%elements)) DEALLOCATE(domainElements%elements)
      DEALLOCATE(domainElements)
    ENDIF
    
    EXITS("DomainTopology_ElementsFinalise")
    RETURN
999 ERRORSEXITS("DomainTopology_ElementsFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE DomainTopology_ElementsFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the element data structures for a domain topology. 
  SUBROUTINE DomainTopology_ElementsInitialise(domainTopology,err,error,*)

    !Argument variables
    TYPE(DomainTopologyType), POINTER :: domainTopology !<A pointer to the domain topology to initialise the elements for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("DomainTopology_ElementsInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(domainTopology)) CALL FlagError("Domain topology is not associated.",err,error,*998)
    IF(ASSOCIATED(domainTopology%elements)) CALL FlagError("Domain already has topology elements associated.",err,error,*998)
     
    ALLOCATE(domainTopology%elements,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate domain topology elements.",err,error,*999)
    domainTopology%elements%numberOfElements=0
    domainTopology%elements%totalNumberOfElements=0
    domainTopology%elements%numberOfGlobalElements=0
    domainTopology%elements%domainTopology=>domainTopology
    domainTopology%elements%maximumNumberOfElementParameters=0
    
    EXITS("DomainTopology_ElementsInitialise")
    RETURN
999 CALL DomainTopology_ElementsFinalise(domainTopology%elements,dummyErr,dummyError,*998)
998 ERRORSEXITS("DomainTopology_ElementsInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE DomainTopology_ElementsInitialise  
  
  !
  !================================================================================================================================
  !

  !>Finalises the topology in the given domain. 
  SUBROUTINE Domain_TopologyFinalise(domainTopology,err,error,*)

   !Argument variables
    TYPE(DomainTopologyType), POINTER :: domainTopology !<A pointer to the domain to finalise the topology for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Domain_TopologyFinalise",err,error,*999)

    IF(ASSOCIATED(domainTopology)) THEN
      CALL DomainTopology_NodesFinalise(domainTopology%nodes,err,error,*999)
      CALL DomainTopology_DOFsFinalise(domainTopology%dofs,err,error,*999)
      CALL DomainTopology_ElementsFinalise(domainTopology%elements,err,error,*999)
      CALL DomainTopology_LinesFinalise(domainTopology%lines,err,error,*999)
      CALL DomainTopology_FacesFinalise(domainTopology%faces,err,error,*999)
      DEALLOCATE(domainTopology)
    ENDIF
 
    EXITS("Domain_TopologyFinalise")
    RETURN
999 ERRORSEXITS("Domain_TopologyFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE Domain_TopologyFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the topology for a given domain. \todo finalise on error
  SUBROUTINE Domain_TopologyInitialise(domain,err,error,*)

    !Argument variables
    TYPE(DomainType), POINTER :: domain !A pointer to the domain to initialise the topology for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("Domain_TopologyInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(domain)) CALL FlagError("Domain is not associated.",err,error,*998)
    IF(ASSOCIATED(domain%topology)) CALL FlagError("Domain already has topology associated.",err,error,*998)
    
    !Allocate domain topology
    ALLOCATE(domain%topology,STAT=err)
    IF(err/=0) CALL FlagError("Domain topology could not be allocated",err,error,*999)
    domain%topology%domain=>domain
    NULLIFY(domain%topology%elements)
    NULLIFY(domain%topology%nodes)
    NULLIFY(domain%topology%dofs)
    NULLIFY(domain%topology%lines)
    NULLIFY(domain%topology%faces)
    !Initialise the topology components
    CALL DomainTopology_ElementsInitialise(domain%topology,err,error,*999)
    CALL DomainTopology_NodesInitialise(domain%topology,err,error,*999)
    CALL DomainTopology_DOFsInitialise(domain%topology,err,error,*999)
    CALL DomainTopology_LinesInitialise(domain%topology,err,error,*999)
    CALL DomainTopology_FacesInitialise(domain%topology,err,error,*999)
    !Initialise the domain topology from the domain mappings and the mesh it came from
    CALL DomainTopology_InitialiseFromMesh(domain%topology,err,error,*999)
    !Calculate the topological information.
    CALL DomainTopology_Calculate(domain%topology,err,error,*999)
    
    EXITS("Domain_TopologyInitialise")
    RETURN
999 CALL Domain_TopologyFinalise(domain%topology,dummyErr,dummyError,*998)
998 ERRORSEXITS("Domain_TopologyInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE Domain_TopologyInitialise
  
  !
  !================================================================================================================================
  !

  !>Finalises a line in the given domain topology and deallocates all memory.
  SUBROUTINE DomainLine_Finalise(domainLine,err,error,*)

    !Argument variables
    TYPE(DomainLineType) :: domainLine !<The domain line to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DomainLine_Finalise",err,error,*999)

    domainLine%number=0
    NULLIFY(domainLine%basis)
    IF(ALLOCATED(domainLine%nodesInLine)) DEALLOCATE(domainLine%nodesInLine)
    IF(ALLOCATED(domainLine%derivativesInLine)) DEALLOCATE(domainLine%derivativesInLine)
 
    EXITS("DomainLine_Finalise")
    RETURN
999 ERRORSEXITS("DomainLine_Finalise",err,error)
    RETURN 1
   
  END SUBROUTINE DomainLine_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the line data structure for a domain topology line.
  SUBROUTINE DomainLine_Initialise(domainLine,err,error,*)

    !Argument variables
    TYPE(DomainLineType) :: domainLine !<The domain line to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DomainLine_Initialise",err,error,*999)

    domainLine%number=0
    NULLIFY(domainLine%basis)
    domainLine%boundaryLine=.FALSE.
    
    EXITS("DomainLine_Initialise")
    RETURN
999 ERRORSEXITS("DomainLine_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE DomainLine_Initialise
  
  !
  !================================================================================================================================
  !

  !>Finalises the lines in the given domain topology. 
  SUBROUTINE DomainTopology_LinesFinalise(domainLines,err,error,*)

    !Argument variables
    TYPE(DomainLinesType), POINTER :: domainLines !<A pointer to the domain topology to finalise the lines for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: lineIdx
    
    ENTERS("DomainTopology_LinesFinalise",err,error,*999)

    IF(ASSOCIATED(domainLines)) THEN
      DO lineIdx=1,SIZE(domainLines%lines,1)
        CALL DomainLine_Finalise(domainLines%lines(lineIdx),err,error,*999)
      ENDDO !lineIdx
      IF(ALLOCATED(domainLines%lines)) DEALLOCATE(domainLines%lines)
      DEALLOCATE(domainLines)
    ENDIF

    EXITS("DomainTopology_LinesFinalise")
    RETURN
999 ERRORSEXITS("DomainTopology_LinesFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE DomainTopology_LinesFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the line data structures for a domain topology. 
  SUBROUTINE DomainTopology_LinesInitialise(domainTopology,err,error,*)

    !Argument variables
    TYPE(DomainTopologyType), POINTER :: domainTopology !<A pointer to the domain topology to initialise the lines for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("DomainTopology_LinesInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(domainTopology)) CALL FlagError("Domain topology is not associated.",err,error,*998)
    IF(ASSOCIATED(domainTopology%lines)) CALL FlagError("Decomposition already has topology lines associated.",err,error,*998)
     
    ALLOCATE(domainTopology%lines,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate topology lines.",err,error,*999)
    domainTopology%lines%numberOfLines=0
    domainTopology%lines%domainTopology=>domainTopology
    
    EXITS("DomainTopology_LinesInitialise")
    RETURN
999 CALL DomainTopology_LinesFinalise(domainTopology%lines,dummyErr,dummyError,*998)
998 ERRORSEXITS("DomainTopology_LinesInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE DomainTopology_LinesInitialise
  
  !
  !================================================================================================================================
  !
  !>Finalises a face in the given domain and deallocates all memory.
  SUBROUTINE DomainFace_Finalise(domainFace,err,error,*)

    !Argument variables
    TYPE(DomainFaceType) :: domainFace !<The domain face to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DomainFace_Finalise",err,error,*999)

    domainFace%number=0
    NULLIFY(domainFace%basis)
    IF(ALLOCATED(domainFace%nodesInFace)) DEALLOCATE(domainFace%nodesInFace)
    IF(ALLOCATED(domainFace%derivativesInFace)) DEALLOCATE(domainFace%derivativesInFace)
 
    EXITS("DomainFace_Finalise")
    RETURN
999 ERRORSEXITS("DomainFace_Finalise",err,error)
    RETURN 1
   
  END SUBROUTINE DomainFace_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the face data structure for a domain topology face.
  SUBROUTINE DomainFace_Initialise(domainFace,err,error,*)

    !Argument variables
    TYPE(DomainFaceType) :: domainFace !<The domain face to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DomainFace_Initialise",err,error,*999)

    domainFace%number=0
    NULLIFY(domainFace%basis)
    domainFace%boundaryFace=.FALSE.
    
    EXITS("DomainFace_Initialise")
    RETURN
999 ERRORSEXITS("DomainFace_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE DomainFace_Initialise
  
  !
  !================================================================================================================================
  !

  !>Finalises the faces in the given domain topology.
  SUBROUTINE DomainTopology_FacesFinalise(domainFaces,err,error,*)

    !Argument variables
    TYPE(DomainFacesType), POINTER :: domainFaces !<A pointer to the domain faces to finalise 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: faceIdx
    
    ENTERS("DomainTopology_FacesFinalise",err,error,*999)

    IF(ASSOCIATED(domainFaces)) THEN
      DO faceIdx=1,SIZE(domainFaces%faces,1)
        CALL DomainFace_Finalise(domainFaces%faces(faceIdx),err,error,*999)
      ENDDO !faceIdx
      IF(ALLOCATED(domainFaces%faces)) DEALLOCATE(domainFaces%faces)
      DEALLOCATE(domainFaces)
    ENDIF
  
    EXITS("DomainTopology_FacesFinalise")
    RETURN
999 ERRORSEXITS("DomainTopology_FacesFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE DomainTopology_FacesFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the face data structures for a domain topology. 
  SUBROUTINE DomainTopology_FacesInitialise(domainTopology,err,error,*)

    !Argument variables
    TYPE(DomainTopologyType), POINTER :: domainTopology !<A pointer to the domain topology to initialise the faces for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("DomainTopology_FacesInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(domainTopology)) CALL FlagError("Domain topology is not associated.",err,error,*998)    
    IF(ASSOCIATED(domainTopology%faces)) CALL FlagError("Domain topology already has faces associated.",err,error,*998)
     
    ALLOCATE(domainTopology%faces,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate domain topology faces.",err,error,*999)
    domainTopology%faces%numberOfFaces=0
    domainTopology%faces%domainTopology=>domainTopology
   
    EXITS("DomainTopology_FacesInitialise")
    RETURN
999 CALL DomainTopology_FacesFinalise(domainTopology%faces,dummyErr,dummyError,*998)
998 ERRORSEXITS("DomainTopology_FacesInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE DomainTopology_FacesInitialise
  
  !
  !================================================================================================================================
  !

  !>Finalises the given domain node derivative and deallocates all memory.
  SUBROUTINE DomainNodeDerivative_Finalise(domainNodeDerivative,err,error,*)

    !Argument variables
    TYPE(DomainNodeDerivativeType), INTENT(INOUT) :: domainNodeDerivative !<The domain node to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DomainNodeDerivative_Finalise",err,error,*999)

    IF(ALLOCATED(domainNodeDerivative%userVersionNumbers)) DEALLOCATE(domainNodeDerivative%userVersionNumbers)
    IF(ALLOCATED(domainNodeDerivative%dofIndex)) DEALLOCATE(domainNodeDerivative%dofIndex)

    EXITS("DomainNodeDerivative_Finalise")
    RETURN
999 ERRORSEXITS("DomainNodeDerivative_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE DomainNodeDerivative_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the given domain node derivative.
  SUBROUTINE DomainNodeDerivative_Initialise(domainNodeDerivative,err,error,*)

    !Argument variables
    TYPE(DomainNodeDerivativeType), INTENT(OUT) :: domainNodeDerivative !<The domain node derivative to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DomainNodeDerivative_Initialise",err,error,*999)

    domainNodeDerivative%numberOfVersions=0
    domainNodeDerivative%globalDerivativeIndex=0
    domainNodeDerivative%partialDerivativeIndex=0

    EXITS("DomainNodeDerivative_Initialise")
    RETURN
999 ERRORSEXITS("DomainNodeDerivative_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE DomainNodeDerivative_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises the given domain node and deallocates all memory.
  SUBROUTINE DomainNode_Finalise(domainNode,err,error,*)

    !Argument variables
    TYPE(DomainNodeType), INTENT(INOUT) :: domainNode !<The domain node to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: derivativeIdx

    ENTERS("DomainNode_Finalise",err,error,*999)

    IF(ALLOCATED(domainNode%derivatives)) THEN
      DO derivativeIdx=1,SIZE(domainNode%derivatives,1)
        CALL DomainNodeDerivative_Finalise(domainNode%derivatives(derivativeIdx),err,error,*999)
      ENDDO !derivativeIdx
      DEALLOCATE(domainNode%derivatives)
    ENDIF
    IF(ALLOCATED(domainNode%surroundingElements)) DEALLOCATE(domainNode%surroundingElements)
    IF(ALLOCATED(domainNode%nodeLines)) DEALLOCATE(domainNode%nodeLines)
 
    EXITS("DomainNode_Finalise")
    RETURN
999 ERRORSEXITS("DomainNode_Finalise",err,error)
    RETURN 1
   
  END SUBROUTINE DomainNode_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the given domain node.
  SUBROUTINE DomainNode_Initialise(domainNode,err,error,*)

    !Argument variables
    TYPE(DomainNodeType) :: domainNode !<The domain node to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DomainNode_Initialise",err,error,*999)

    domainNode%localNumber=0
    domainNode%meshNumber=0
    domainNode%globalNumber=0
    domainNode%userNumber=0
    domainNode%numberOfSurroundingElements=0
    domainNode%numberOfNodeLines=0
    domainNode%boundaryNode=.FALSE.
    
    EXITS("DomainNode_Initialise")
    RETURN
999 ERRORSEXITS("DomainNode_Initialise",err,error)
    RETURN 1
   
  END SUBROUTINE DomainNode_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises the nodes in a domain. 
  SUBROUTINE DomainTopology_NodesFinalise(domainNodes,err,error,*)

    !Argument variables
    TYPE(DomainNodesType), POINTER :: domainNodes !<A pointer to the domain nodes to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: nodeIdx

    ENTERS("DomainTopology_NodesFinalise",err,error,*999)

    IF(ASSOCIATED(domainNodes)) THEN
      DO nodeIdx=1,SIZE(domainNodes%nodes,1)
        CALL DomainNode_Finalise(domainNodes%nodes(nodeIdx),err,error,*999)
      ENDDO !nodeIdx
      IF(ALLOCATED(domainNodes%nodes)) DEALLOCATE(domainNodes%nodes)
      IF(ASSOCIATED(domainNodes%nodesTree)) CALL Tree_Destroy(domainNodes%nodesTree,err,error,*999)
      DEALLOCATE(domainNodes)
    ENDIF
 
    EXITS("DomainTopology_NodesFinalise")
    RETURN
999 ERRORSEXITS("DomainTopology_NodesFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE DomainTopology_NodesFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the nodes data structures for a domain topology.
  SUBROUTINE DomainTopology_NodesInitialise(domainTopology,err,error,*)

    !Argument variables
    TYPE(DomainTopologyType), POINTER :: domainTopology !<A pointer to the domain topology to initialise the nodes for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("DomainTopology_NodesInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(domainTopology)) CALL FlagError("Domain topology is not associated.",err,error,*998)
    IF(ASSOCIATED(domainTopology%nodes)) CALL FlagError("Domain topology already has nodes associated.",err,error,*998)
      
    ALLOCATE(domainTopology%nodes,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate domain topology nodes.",err,error,*999)
    domainTopology%nodes%numberOfNodes=0
    domainTopology%nodes%totalNumberOfNodes=0
    domainTopology%nodes%numberOfGlobalNodes=0
    domainTopology%nodes%maximumNumberOfDerivatives=0
    domainTopology%nodes%domainTopology=>domainTopology
    NULLIFY(domainTopology%nodes%nodesTree)
   
    EXITS("DomainTopology_NodesInitialise")
    RETURN
999 CALL DomainTopology_NodesFinalise(domainTopology%nodes,dummyErr,dummyError,*998)
998 ERRORSEXITS("DomainTopology_NodesInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE DomainTopology_NodesInitialise
  
  !
  !================================================================================================================================
  !

  !>Calculates the element numbers surrounding a node for a domain.
  SUBROUTINE DomainTopology_NodesSurroundingElementsCalculate(domainTopology,err,error,*)

    !Argument variables
    TYPE(DomainTopologyType), POINTER :: domainTopology !<A pointer to the domain topology to calculate the elements surrounding the nodes for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementIdx,elementNumber,insertPosition,localNodeIdx,nodeIdx,nodeNumber,surroundingElementNumber
    INTEGER(INTG), ALLOCATABLE :: newSurroundingElements(:)
    LOGICAL :: foundElement
    TYPE(BasisType), POINTER :: basis
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainNodesType), POINTER :: domainNodes

    ENTERS("DomainTopology_NodesSurroundingElementsCalculate",err,error,*999)
    
    IF(.NOT.ASSOCIATED(domainTopology)) CALL FlagError("Domain topology is not associated.",err,error,*999)

    NULLIFY(domainElements)
    CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
    NULLIFY(domainNodes)
    CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
    IF(.NOT.ALLOCATED(domainNodes%nodes)) CALL FlagError("Domain nodes nodes is not allocated.",err,error,*999)
    
    DO nodeIdx=1,domainNodes%totalNumberOfNodes
      domainNodes%nodes(nodeIdx)%numberOfSurroundingElements=0
      IF(ALLOCATED(domainNodes%nodes(nodeIdx)%surroundingElements)) DEALLOCATE(domainNodes%nodes(nodeIdx)%surroundingElements)
    ENDDO !nodeIdx
    DO elementIdx=1,domainElements%totalNumberOfElements
      NULLIFY(basis)
      CALL DomainElements_BasisGet(domainElements,elementIdx,basis,err,error,*999)
      DO localNodeIdx=1,basis%numberOfNodes
        nodeNumber=domainElements%elements(elementIdx)%elementNodes(localNodeIdx)
        foundElement=.FALSE.
        elementNumber=1
        insertPosition=1
        DO WHILE(elementNumber<=domainNodes%nodes(nodeNumber)%numberOfSurroundingElements.AND..NOT.foundElement)
          surroundingElementNumber=domainNodes%nodes(nodeNumber)%surroundingElements(elementNumber)
          IF(surroundingElementNumber==elementIdx) THEN
            foundElement=.TRUE.
          ENDIF
          elementNumber=elementNumber+1
          IF(elementIdx>=surroundingElementNumber) THEN
            insertPosition=elementNumber
          ENDIF
        ENDDO
        IF(.NOT.foundElement) THEN
          !Insert element into surrounding elements
          ALLOCATE(newSurroundingElements(domainNodes%nodes(nodeNumber)%numberOfSurroundingElements+1),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate new surrounding elements",err,error,*999)
          IF(ALLOCATED(domainNodes%nodes(nodeNumber)%surroundingElements)) THEN
            newSurroundingElements(1:insertPosition-1)=domainNodes%nodes(nodeNumber)%surroundingElements(1:insertPosition-1)
            newSurroundingElements(insertPosition)=elementIdx
            newSurroundingElements(insertPosition+1:domainNodes%nodes(nodeNumber)%numberOfSurroundingElements+1)= &
              & domainNodes%nodes(nodeNumber)%surroundingElements(insertPosition:domainNodes%nodes(nodeNumber)% &
              & numberOfSurroundingElements)
            DEALLOCATE(domainNodes%nodes(nodeNumber)%surroundingElements)
          ELSE
            newSurroundingElements(1)=elementIdx
          ENDIF
          CALL MOVE_ALLOC(newSurroundingElements,domainNodes%nodes(nodeNumber)%surroundingElements)
          domainNodes%nodes(nodeNumber)%numberOfSurroundingElements=domainNodes%nodes(nodeNumber)%numberOfSurroundingElements+1
        ENDIF
      ENDDO !localNodeIdx
    ENDDO !elementIdx

    EXITS("DomainTopology_NodesSurroundingElementsCalculate")
    RETURN
999 IF(ALLOCATED(newSurroundingElements)) DEALLOCATE(newSurroundingElements)
    ERRORS("DomainTopology_NodesSurroundingElementsCalculate",err,error)
    EXITS("DomainTopology_NodesSurroundingElementsCalculate")
    RETURN 1
    
  END SUBROUTINE DomainTopology_NodesSurroundingElementsCalculate
  
  !
  !================================================================================================================================
  !

END MODULE DecompositionRoutines
