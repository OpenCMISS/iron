!> \file
!> \author Chris Bradley
!> \brief This module handles all node routines.
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

!> This module handles all node routines.
MODULE NodeRoutines

  USE BaseRoutines
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE Kinds
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

  INTERFACE NODE_CHECK_EXISTS
    MODULE PROCEDURE Node_CheckExists
  END INTERFACE NODE_CHECK_EXISTS

  INTERFACE NODES_CREATE_START
    MODULE PROCEDURE Nodes_CreateStartRegion
    MODULE PROCEDURE Nodes_CreateStartInterface
  END INTERFACE NODES_CREATE_START

  !>Starts the process of creating nodes for an interface or region
  INTERFACE Nodes_CreateStart
    MODULE PROCEDURE Nodes_CreateStartRegion
    MODULE PROCEDURE Nodes_CreateStartInterface
  END INTERFACE Nodes_CreateStart

  INTERFACE NODES_CREATE_FINISH
    MODULE PROCEDURE Nodes_CreateFinish
  END INTERFACE NODES_CREATE_FINISH

  !>Initialises nodes for an interface or region
  INTERFACE Nodes_Initialise
    MODULE PROCEDURE Nodes_InitialiseRegion
    MODULE PROCEDURE Nodes_InitialiseInterface
  END INTERFACE Nodes_Initialise

  INTERFACE NODES_LABEL_GET
    MODULE PROCEDURE Nodes_LabelGetC
    MODULE PROCEDURE Nodes_LabelGetVS
  END INTERFACE NODES_LABEL_GET

  !>Gets the label for a node identified by a given global number.
  INTERFACE Nodes_LabelGet
    MODULE PROCEDURE Nodes_LabelGetC
    MODULE PROCEDURE Nodes_LabelGetVS
  END INTERFACE Nodes_LabelGet

   INTERFACE NODES_LABEL_SET
    MODULE PROCEDURE Nodes_LabelSetC
    MODULE PROCEDURE Nodes_LabelSetVS
  END INTERFACE NODES_LABEL_SET

  !>Changes/sets the label for a node identified by a given global number.
  INTERFACE Nodes_LabelSet
    MODULE PROCEDURE Nodes_LabelSetC
    MODULE PROCEDURE Nodes_LabelSetVS
  END INTERFACE Nodes_LabelSet

  INTERFACE NODES_NUMBER_OF_NODES_GET
    MODULE PROCEDURE Nodes_NumberOfNodesGet
  END INTERFACE NODES_NUMBER_OF_NODES_GET

  INTERFACE NODES_USER_NUMBER_GET
    MODULE PROCEDURE Nodes_UserNumberGet
  END INTERFACE NODES_USER_NUMBER_GET

  INTERFACE NODES_USER_NUMBER_SET
    MODULE PROCEDURE Nodes_UserNumberSet
  END INTERFACE NODES_USER_NUMBER_SET

  PUBLIC NODE_CHECK_EXISTS

  PUBLIC Node_CheckExists

  PUBLIC NODES_CREATE_FINISH,NODES_CREATE_START

  PUBLIC Nodes_CreateFinish,Nodes_CreateStart

  PUBLIC Nodes_Destroy

  PUBLIC NODES_LABEL_GET,NODES_LABEL_SET

  PUBLIC Nodes_LabelGet,Nodes_LabelSet

  PUBLIC NODES_NUMBER_OF_NODES_GET

  PUBLIC NODES_NumberOfNodesGet
  
  PUBLIC NODES_USER_NUMBER_GET,NODES_USER_NUMBER_SET

  PUBLIC Nodes_UserNumberGet,Nodes_UserNumberSet
  
  PUBLIC Nodes_UserNumbersAllSet


CONTAINS

  !
  !================================================================================================================================
  !

  !>Checks that a user node number is defined on the specified region.
  SUBROUTINE Node_CheckExists(nodes,userNumber,nodeExists,globalNumber,err,error,*)

    !Argument variables
    TYPE(NodesType), POINTER :: nodes !<A pointer to the nodes to check
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user node number to check if it exists
    LOGICAL, INTENT(OUT) :: nodeExists !<On exit, is .TRUE. if the node user number exists in the region, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: globalNumber !<On exit, if the node exists the global number corresponding to the user node number. If the node does not exist then global number will be 0.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(TREE_NODE_TYPE), POINTER :: treeNode
   
    ENTERS("Node_CheckExists",err,error,*999)

    IF(.NOT.ASSOCIATED(nodes)) CALL FlagError("Nodes is not associated.",err,error,*999)

    nodeExists=.FALSE.
    globalNumber=0
    NULLIFY(treeNode)
    CALL Tree_Search(nodes%nodesTree,userNumber,treeNode,err,error,*999)
    IF(ASSOCIATED(treeNode)) THEN
      CALL Tree_NodeValueGet(nodes%nodesTree,treeNode,globalNumber,err,error,*999)
      nodeExists=.TRUE.
    ENDIF

    EXITS("Node_CheckExists")
    RETURN
999 ERRORSEXITS("Node_CheckExists",err,error)
    RETURN 1
    
  END SUBROUTINE Node_CheckExists
  
  !
  !================================================================================================================================
  !

  !>Finalises a node and deallocates all memory
  SUBROUTINE Node_Finalise(node,err,error,*)
    
    !Argument variables
    TYPE(NodeType), INTENT(OUT) :: node !<The node to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Node_Finalise",err,error,*999)

    node%globalNumber=0
    node%userNumber=0
    node%label=""
    
    EXITS("Node_Finalise")
    RETURN
999 ERRORSEXITS("Node_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE Node_Finalise
  
  !
  !================================================================================================================================
  !

  !>Initialises a node 
  SUBROUTINE Node_Initialise(node,err,error,*)
    
    !Argument variables
    TYPE(NodeType), INTENT(OUT) :: node !<The node to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Node_Initialise",err,error,*999)

    node%globalNumber=0
    node%userNumber=0
    node%label=""
    
    EXITS("Node_Initialise")
    RETURN
999 ERRORSEXITS("Node_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE Node_Initialise
  
  !
  !================================================================================================================================
  !

  !>Finishes the process of creating nodes in the region. \see OpenCMISS::Iron::cmfe_Nodes_CreateFinish
  SUBROUTINE Nodes_CreateFinish(nodes,err,error,*)

    !Argument variables
    TYPE(NodesType), POINTER :: nodes !<A pointer to the nodes to be finished
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: nodeIdx
    
    ENTERS("Nodes_CreateFinish",err,error,*999)

    IF(.NOT.ASSOCIATED(nodes)) CALL FlagError("Nodes is not associated.",err,error,*999)
    IF(nodes%nodesFinished) CALL FlagError("Nodes have already been finished.",err,error,*999)
    
    nodes%nodesFinished=.TRUE.
    
    IF(diagnostics1) THEN
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Number of nodes = ",nodes%numberOfNodes,err,error,*999)
      DO nodeIdx=1,nodes%numberOfNodes
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Node = ",nodeIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Global number    = ",nodes%nodes(nodeIdx)%globalNumber,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    User number      = ",nodes%nodes(nodeIdx)%userNumber,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Label            = ",nodes%nodes(nodeIdx)%label,err,error,*999)
      ENDDO !nodeIdx
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"User->Global number tree:",err,error,*999)
      CALL Tree_Output(DIAGNOSTIC_OUTPUT_TYPE,nodes%nodesTree,err,error,*999)
    ENDIF

    EXITS("Nodes_CreateFinish")
    RETURN
999 ERRORSEXITS("Nodes_CreateFinish",err,error)
    RETURN 1
   
  END SUBROUTINE Nodes_CreateFinish
    
  !
  !================================================================================================================================
  !

  !>Starts the process of creating generic nodes.
  SUBROUTINE Nodes_CreateStartGeneric(nodes,numberOfNodes,err,error,*)

    !Argument variables
    TYPE(NodesType), POINTER :: nodes !<The nodes pointer
    INTEGER(INTG), INTENT(IN) :: numberOfNodes !<The number of nodes to create
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: insertStatus,nodeIdx
    TYPE(VARYING_STRING) :: localError

    ENTERS("Nodes_CreateStartGeneric",err,error,*999)

    IF(.NOT.ASSOCIATED(nodes)) CALL FlagError("Nodes is not associated.",err,error,*999)
    
    IF(numberOfNodes<=0) THEN
      localError="The specified number of nodes of "//TRIM(NumberToVString(numberOfNodes,"*",err,error))// &
        & " is invalid. The number of nodes must be > 0."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    ALLOCATE(nodes%nodes(numberOfNodes),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate nodes nodes.",err,error,*999)
    nodes%numberOfNodes=numberOfNodes
    CALL Tree_CreateStart(nodes%nodesTree,err,error,*999)
    CALL Tree_InsertTypeSet(nodes%nodesTree,TREE_NO_DUPLICATES_ALLOWED,err,error,*999)
    CALL Tree_CreateFinish(nodes%nodesTree,err,error,*999)
    !Set default node numbers
    DO nodeIdx=1,nodes%numberOfNodes
      CALL Node_Initialise(nodes%nodes(nodeIdx),err,error,*999)
      nodes%nodes(nodeIdx)%globalNumber=nodeIdx
      nodes%nodes(nodeIdx)%userNumber=nodeIdx
      CALL Tree_ItemInsert(nodes%nodesTree,nodeIdx,nodeIdx,insertStatus,err,error,*999)
    ENDDO !nodeIdx
    
    EXITS("Nodes_CreateStartGeneric")
    RETURN
999 ERRORSEXITS("Nodes_CreateStartGeneric",err,error)
    RETURN 1
    
  END SUBROUTINE Nodes_CreateStartGeneric

  !
  !================================================================================================================================
  !

  !>Starts the process of creating nodes in an interface. \see OpenCMISS::Iron::cmfe_Nodes_CreateStart
  SUBROUTINE Nodes_CreateStartInterface(interface,numberOfNodes,nodes,err,error,*)

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: interface !<A pointer to the interface in which to create the nodes
    INTEGER(INTG), INTENT(IN) :: numberOfNodes !<The number of nodes to create
    TYPE(NodesType), POINTER :: nodes !<On exit, a pointer to the created nodes. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("Nodes_CreateStartInterface",err,error,*998)

    IF(.NOT.ASSOCIATED(interface)) CALL FlagError("Interface is not associated.",err,error,*998)
    IF(ASSOCIATED(interface%nodes)) CALL FlagError("Interface already has nodes associated.",err,error,*998)
    IF(ASSOCIATED(nodes)) CALL FlagError("Nodes is already associated.",err,error,*999)
        
    !Initialise the nodes for the interface
    CALL Nodes_Initialise(interface,err,error,*999)
    !Create the nodes 
    CALL Nodes_CreateStartGeneric(interface%nodes,numberOfNodes,err,error,*999)
    !Return the pointer        
    nodes=>interface%nodes
    
    EXITS("Nodes_CreateStartInterface")
    RETURN
999 CALL Nodes_Finalise(interface%nodes,dummyErr,dummyError,*998)    
998 ERRORSEXITS("Nodes_CreateStartInterface",err,error)
    RETURN 1
   
  END SUBROUTINE Nodes_CreateStartInterface

  !
  !================================================================================================================================
  !

  !>Starts the process of creating nodes in the region.  \see OpenCMISS::Iron::cmfe_Nodes_CreateStart
  SUBROUTINE Nodes_CreateStartRegion(region,numberOfNodes,nodes,err,error,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: region !<A pointer to the region in which to create the nodes
    INTEGER(INTG), INTENT(IN) :: numberOfNodes !<The number of nodes to create
    TYPE(NodesType), POINTER :: nodes !<On exit, a pointer to the created nodes
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("Nodes_CreateStartRegion",err,error,*998)

    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated.",err,error,*998)
    IF(ASSOCIATED(region%nodes)) CALL FlagError("Region already has nodes associated.",err,error,*998)
    IF(ASSOCIATED(nodes)) CALL FlagError("Nodes is already associated.",err,error,*998)
    
    !Initialise the nodes for the region    
    CALL Nodes_Initialise(region,err,error,*999)
    !Create the generic nodes
    CALL Nodes_CreateStartGeneric(region%nodes,numberOfNodes,err,error,*999)
    !Return the pointer        
    nodes=>region%nodes
   
    EXITS("Nodes_CreateStartRegion")
    RETURN
999 CALL Nodes_Finalise(region%nodes,dummyErr,dummyError,*998)    
998 ERRORSEXITS("Nodes_CreateStartRegion",err,error)
    RETURN 1
   
  END SUBROUTINE Nodes_CreateStartRegion

  !
  !================================================================================================================================
  !

  !>Destroys nodes. \see OpenCMISS::Iron::cmfe_Nodes_Destroy
  SUBROUTINE Nodes_Destroy(nodes,err,error,*)

    !Argument variables
    TYPE(NodesType), POINTER :: nodes !<A pointer to the nodes to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("Nodes_Destroy",err,error,*999)

    IF(.NOT.ASSOCIATED(nodes)) CALL FlagError("Nodes is not associated.",err,error,*999)
    
    IF(ASSOCIATED(nodes%region)) THEN
      NULLIFY(nodes%region%nodes)
    ELSE
      IF(ASSOCIATED(nodes%interface)) THEN
        NULLIFY(nodes%interface%nodes)
      ELSE
        CALL FlagError("Nodes region and interface are not associated.",err,error,*999)
      ENDIF
    ENDIF
    
    CALL Nodes_Finalise(nodes,err,error,*999)
    
    EXITS("Nodes_Destroy")
    RETURN
999 ERRORSEXITS("Nodes_Destroy",err,error)
    RETURN 1
   
  END SUBROUTINE Nodes_Destroy
    
  !
  !===============================================================================================================================
  !

  !>Finalises the nodes and deallocates any memory. 
  SUBROUTINE Nodes_Finalise(nodes,err,error,*)

    !Argument variables
    TYPE(NodesType), POINTER :: nodes !<A pointer to the nodes to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: nodeIdx

    ENTERS("Nodes_Finalise",err,error,*999)

    IF(ASSOCIATED(nodes)) THEN
      IF(ALLOCATED(nodes%nodes)) THEN
        DO nodeIdx=1,SIZE(nodes%nodes,1)
          CALL Node_Finalise(nodes%nodes(nodeIdx),err,error,*999)
        ENDDO !nodeIdx
        DEALLOCATE(nodes%nodes)
      ENDIF
      IF(ASSOCIATED(nodes%nodesTree)) CALL Tree_Destroy(nodes%nodesTree,err,error,*999)
      DEALLOCATE(nodes)
    ENDIF
    
    EXITS("Nodes_Finalise")
    RETURN
999 ERRORSEXITS("Nodes_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE Nodes_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the nodes.
  SUBROUTINE Nodes_InitialiseGeneric(nodes,err,error,*)

    !Argument variables
    TYPE(NodesType), POINTER :: nodes !<A pointer to the nodes to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Nodes_InitialiseGeneric",err,error,*999)

    IF(.NOT.ASSOCIATED(nodes)) CALL FlagError("Nodes is not associated.",err,error,*999)
    
    NULLIFY(nodes%region)
    NULLIFY(nodes%interface)
    nodes%nodesFinished=.FALSE.
    nodes%numberOfNodes=0
    NULLIFY(nodes%nodesTree)
   
    EXITS("Nodes_InitialiseGeneric")
    RETURN
999 ERRORSEXITS("Nodes_InitialiseGeneric",err,error)
    RETURN 1
    
  END SUBROUTINE Nodes_InitialiseGeneric

  !
  !================================================================================================================================
  !

  !>Initialises the nodes in a given interface.
  SUBROUTINE Nodes_InitialiseInterface(INTERFACE,err,error,*)

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: interface !<A pointer to the interface to initialise the nodes for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("Nodes_InitialiseInterface",err,error,*998)

    IF(.NOT.ASSOCIATED(interface)) CALL FlagError("Interface is not associated.",err,error,*998)
    IF(ASSOCIATED(interface%nodes)) CALL FlagError("Interface already has associated nodes.",err,error,*998)
    
    ALLOCATE(interface%nodes,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate interface nodes.",err,error,*999)
    CALL Nodes_InitialiseGeneric(interface%nodes,err,error,*999)
    interface%nodes%interface=>interface
    
    EXITS("Nodes_InitialiseInterface")
    RETURN
999 CALL Nodes_Finalise(interface%nodes,dummyErr,dummyError,*998)
998 ERRORSEXITS("Nodes_InitialiseInterface",err,error)
    RETURN 1
    
  END SUBROUTINE Nodes_InitialiseInterface

  !
  !================================================================================================================================
  !

  !>Initialises the nodes in a given region.
  SUBROUTINE Nodes_InitialiseRegion(region,err,error,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: region !<A pointer to the region to initialise the nodes for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("Nodes_InitialiseRegion",err,error,*998)

    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated.",err,error,*998)
    IF(ASSOCIATED(region%nodes)) CALL FlagError("Region has associated nodes.",err,error,*998)
    
    ALLOCATE(region%nodes,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate region nodes.",err,error,*999)
    CALL Nodes_InitialiseGeneric(region%nodes,err,error,*999)
    region%nodes%region=>region
          
    EXITS("Nodes_InitialiseRegion")
    RETURN
999 CALL Nodes_Finalise(region%nodes,dummyErr,dummyError,*998)
998 ERRORSEXITS("Nodes_InitialiseRegion",err,error)
    RETURN 1
    
  END SUBROUTINE Nodes_InitialiseRegion

  !
  !================================================================================================================================
  !

  !>Gets the character label for a node identified by a given global number. \see OpenCMISS::Iron::cmfe_Nodes_LabelGet
  SUBROUTINE Nodes_LabelGetC(nodes,globalNumber,label,err,error,*)

    !Argument variables
    TYPE(NodesType), POINTER :: nodes !<A pointer to the nodes to get the label for
    INTEGER(INTG), INTENT(IN) :: globalNumber !<The global number to get the label for
    CHARACTER(LEN=*), INTENT(OUT) :: label !<On exit, the label of the specified global node
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER :: cLength,vsLength
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Nodes_LabelGetC",err,error,*999)

    IF(.NOT.ASSOCIATED(nodes)) CALL FlagError("Nodes is not associated.",err,error,*999)
    IF(.NOT.nodes%nodesFinished) CALL FlagError("Nodes have not been finished.",err,error,*999)
    IF(globalNumber<1.OR.globalNumber>nodes%numberOfNodes) THEN
      localError="The specified global node number of "//TRIM(NumberToVString(globalNumber,"*",err,error))// &
        & " is invalid. The global node number should be between 1 and "// &
        & TRIM(NumberToVString(nodes%numberOfNodes,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    cLength=LEN(label)
    vsLength=LEN_TRIM(nodes%nodes(globalNumber)%label)
    IF(cLength>vsLength) THEN
      label=CHAR(LEN_TRIM(nodes%nodes(globalNumber)%label))
    ELSE
      label=CHAR(nodes%nodes(globalNumber)%label,cLength)
    ENDIF
    
    EXITS("Nodes_LabelGetC")
    RETURN
999 ERRORSEXITS("Nodes_LabelGetC",err,error)    
    RETURN 1
   
  END SUBROUTINE Nodes_LabelGetC
        
  !
  !================================================================================================================================
  !

  !>Gets the varying string label for a node identified by a given global number. \see OpenCMISS::Iron::cmfe_Nodes_LabelGet
  SUBROUTINE Nodes_LabelGetVS(nodes,globalNumber,label,err,error,*)

    !Argument variables
    TYPE(NodesType), POINTER :: nodes !<A pointer to the nodes to get the label for
    INTEGER(INTG), INTENT(IN) :: globalNumber !<The global number to get the label for
    TYPE(VARYING_STRING), INTENT(OUT) :: label !<On exit, the label of the specified global node
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Nodes_LabelGetVS",err,error,*999)

    IF(.NOT.ASSOCIATED(nodes)) CALL FlagError("Nodes is not associated.",err,error,*999)
    IF(.NOT.nodes%nodesFinished) CALL FlagError("Nodes have not been finished.",err,error,*999)
    IF(globalNumber<1.OR.globalNumber>nodes%numberOfNodes) THEN
      localError="The specified global node number of "//TRIM(NumberToVString(globalNumber,"*",err,error))// &
        & " is invalid. The global node number should be between 1 and "// &
        & TRIM(NumberToVString(nodes%numberOfNodes,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    label=nodes%nodes(globalNumber)%label
    
    EXITS("Nodes_LabelGetVS")
    RETURN
999 ERRORSEXITS("Nodes_LabelGetVS",err,error)    
    RETURN 1
   
  END SUBROUTINE Nodes_LabelGetVS
        
  !
  !================================================================================================================================
  !

  !>Changes/sets the character label for a node identified by a given global number. \see OpenCMISS::Iron::cmfe_Nodes_LabelSet
  SUBROUTINE Nodes_LabelSetC(nodes,globalNumber,label,err,error,*)

    !Argument variables
    TYPE(NodesType), POINTER :: nodes !<A pointer to the nodes to set the label for
    INTEGER(INTG), INTENT(IN) :: globalNumber !<The global number to set the label for
    CHARACTER(LEN=*), INTENT(IN) :: label !<The label to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Nodes_LabelSetC",err,error,*999)

    IF(.NOT.ASSOCIATED(nodes)) CALL FlagError("Nodes is not associated.",err,error,*999)
    IF(globalNumber<1.OR.globalNumber>nodes%numberOfNodes) THEN
      localError="The specified global node number of "//TRIM(NumberToVString(globalNumber,"*",err,error))// &
        & " is invalid. The global node number should be between 1 and "// &
        & TRIM(NumberToVString(nodes%numberOfNodes,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    nodes%nodes(globalNumber)%label=label
    
    EXITS("Nodes_LabelSetC")
    RETURN
999 ERRORSEXITS("Nodes_LabelSetC",err,error)    
    RETURN 1
   
  END SUBROUTINE Nodes_LabelSetC
        
  !
  !================================================================================================================================
  !

  !>Changes/sets the varying string label for a node identified by a given global number. \see OpenCMISS::Iron::cmfe_Nodes_LabelSet
  SUBROUTINE Nodes_LabelSetVS(nodes,globalNumber,label,err,error,*)

    !Argument variables
    TYPE(NodesType), POINTER :: nodes !<A pointer to the nodes to set the label for
    INTEGER(INTG), INTENT(IN) :: globalNumber !<The global number to set the label for
    TYPE(VARYING_STRING), INTENT(IN) :: label !<The label to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Nodes_LabelSetVS",err,error,*999)

    IF(.NOT.ASSOCIATED(nodes)) CALL FlagError("Nodes is not associated.",err,error,*999)
    IF(globalNumber<1.OR.globalNumber>nodes%numberOfNodes) THEN
      localError="The specified global node number of "//TRIM(NumberToVString(globalNumber,"*",err,error))// &
        & " is invalid. The global node number should be between 1 and "// &
        & TRIM(NumberToVString(nodes%numberOfNodes,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    nodes%nodes(globalNumber)%label=label
    
    EXITS("Nodes_LabelSetVS")
    RETURN
999 ERRORSEXITS("Nodes_LabelSetVS",err,error)    
    RETURN 1
   
  END SUBROUTINE Nodes_LabelSetVS
        
  !
  !================================================================================================================================
  !

  !>Returns the number of nodes. \see OpenCMISS::Iron::cmfe_Nodes_NumberOfNodesGet
  SUBROUTINE Nodes_NumberOfNodesGet(nodes,numberOfNodes,err,error,*)

    !Argument variables
    TYPE(NodesType), POINTER :: nodes !<A pointer to the nodes to get the number of nodes for
    INTEGER(INTG), INTENT(OUT) :: numberOfNodes !<On return, the number of nodes
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("Nodes_NumberOfNodesGet",err,error,*999)

    IF(.NOT.ASSOCIATED(nodes)) CALL FlagError("Nodes is not associated.",err,error,*999)
    IF(.NOT.nodes%nodesFinished) CALL FlagError("Nodes have not been finished.",err,error,*999)
    
    numberOfNodes=nodes%numberOfNodes
    
    EXITS("Nodes_NumberOfNodesGet")
    RETURN
999 ERRORSEXITS("Nodes_NumberOfNodesGet",err,error)    
    RETURN 1
   
  END SUBROUTINE Nodes_NumberOfNodesGet
        
  !
  !================================================================================================================================
  !

  !>Gets the user number for a node identified by a given global number. \see OpenCMISS::Iron::cmfe_Nodes_UserNumberGet
  SUBROUTINE Nodes_UserNumberGet(nodes,globalNumber,userNumber,err,error,*)

    !Argument variables
    TYPE(NodesType), POINTER :: nodes !<A pointer to the nodes to get the number for
    INTEGER(INTG), INTENT(IN) :: globalNumber !<The global number to get the user number for
    INTEGER(INTG), INTENT(OUT) :: userNumber !<On exit, the user number of the specified global node
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Nodes_UserNumberGet",err,error,*999)

    IF(.NOT.ASSOCIATED(nodes)) CALL FlagError("Nodes is not associated.",err,error,*999)
    IF(.NOT.nodes%nodesFinished) CALL FlagError("Nodes have not been finished.",err,error,*999)
    IF(globalNumber<1.OR.globalNumber>nodes%numberOfNodes) THEN
      localError="The specified global node number of "//TRIM(NumberToVString(globalNumber,"*",err,error))// &
        & " is invalid. The global node number should be between 1 and "// &
        & TRIM(NumberToVString(nodes%numberOfNodes,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    userNumber=nodes%nodes(globalNumber)%userNumber
    
    EXITS("Nodes_UserNumberGet")
    RETURN
999 ERRORSEXITS("Nodes_UserNumberGet",err,error)    
    RETURN 1
   
  END SUBROUTINE Nodes_UserNumberGet
        
  !
  !================================================================================================================================
  !

  !>Changes/sets the user number for a node identified by a given global number. \see OpenCMISS::Iron::cmfe_Nodes_UserNumberSet
  SUBROUTINE Nodes_UserNumberSet(nodes,globalNumber,userNumber,err,error,*)

    !Argument variables
    TYPE(NodesType), POINTER :: nodes !<A pointer to the nodes to set the number for
    INTEGER(INTG), INTENT(IN) :: globalNumber !<The global number to set the user number for
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: insertStatus,oldGlobalNumber
    LOGICAL :: nodeExists
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Nodes_UserNumberSet",err,error,*999)

    IF(ASSOCIATED(nodes)) CALL FlagError("Nodes is not associated.",err,error,*999)
    IF(nodes%nodesFinished) CALL FlagError("Nodes have been finished.",err,error,*999)
    IF(globalNumber<1.OR.globalNumber>nodes%numberOfNodes) THEN
      localError="The specified global node number of "//TRIM(NumberToVString(globalNumber,"*",err,error))// &
        & " is invalid. The global node number should be between 1 and "// &
        & TRIM(NumberToVString(nodes%numberOfNodes,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    !Check the node user number is not already used
    CALL Node_CheckExists(nodes,userNumber,nodeExists,oldGlobalNumber,err,error,*999)
    IF(nodeExists) THEN
      IF(oldGlobalNumber/=globalNumber) THEN
        localError="The specified node user number of "//TRIM(NumberToVString(userNumber,"*",err,error))// &
          & " is already used by global node number "//TRIM(NumberToVString(oldGlobalNumber,"*",err,error))// &
          & ". User node numbers must be unique."
        CALL FlagError(localError,err,error,*999)
      ENDIF      
    ELSE
      CALL Tree_ItemDelete(nodes%nodesTree,nodes%nodes(globalNumber)%userNumber,err,error,*999)
      CALL Tree_ItemInsert(nodes%nodesTree,userNumber,globalNumber,insertStatus,err,error,*999)
      IF(insertStatus/=TREE_NODE_INSERT_SUCESSFUL) CALL FlagError("Unsucessful nodes tree insert.",err,error,*999)
      nodes%nodes(globalNumber)%userNumber=userNumber
    ENDIF
    
    EXITS("Nodes_UserNumberSet")
    RETURN
999 ERRORSEXITS("Nodes_UserNumberSet",err,error)    
    RETURN 1
   
  END SUBROUTINE Nodes_UserNumberSet

  !
  !================================================================================================================================
  !

  !>Changes/sets the user numbers for all nodes. \see OpenCMISS::Iron::cmfe_Nodes_UserNumbersAllSet
  SUBROUTINE Nodes_UserNumbersAllSet(nodes,userNumbers,err,error,*)

    !Argument variables
    TYPE(NodesType), POINTER :: nodes !<A pointer to the nodes to set the numbers for
    INTEGER(INTG), INTENT(IN) :: userNumbers(:) !<The user numbers to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: nodeIdx,insertStatus
    TYPE(TREE_TYPE), POINTER :: newNodesTree
    TYPE(VARYING_STRING) :: localError

    NULLIFY(newNodesTree)
    
    ENTERS("Nodes_UserNumbersAllSet",err,error,*999)

    IF(.NOT.ASSOCIATED(nodes)) CALL FlagError("Nodes is not associated.",err,error,*999)
    IF(nodes%nodesFinished) CALL FlagError("Nodes have been finished.",err,error,*999)
    IF(SIZE(userNumbers,1)/=nodes%numberOfNodes) THEN
      localError="The number of specified node user numbers of "// &
        TRIM(NumberToVstring(SIZE(userNumbers,1),"*",err,error))// &
        " does not match number of nodes of "// &
        TRIM(NumberToVstring(nodes%numberOfNodes,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    !Check the users numbers to ensure that there are no duplicates                    
    CALL Tree_CreateStart(newNodesTree,err,error,*999)
    CALL Tree_InsertTypeSet(newNodesTree,TREE_NO_DUPLICATES_ALLOWED,err,error,*999)
    CALL Tree_CreateFinish(newNodesTree,err,error,*999)
    DO nodeIdx=1,nodes%numberOfNodes
      CALL Tree_ItemInsert(newNodesTree,userNumbers(nodeIdx),nodeIdx,insertStatus,err,error,*999)
      IF(insertStatus/=TREE_NODE_INSERT_SUCESSFUL) THEN
        localError="The specified user number of "//TRIM(NumberToVstring(userNumbers(nodeIdx),"*",err,error))// &
          & " for global node number "//TRIM(NumberToVstring(nodeIdx,"*",err,error))// &
          & " is a duplicate. The user node numbers must be unique."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDDO !nodeIdx
    CALL Tree_Destroy(nodes%nodesTree,err,error,*999)
    nodes%nodesTree=>newNodesTree
    NULLIFY(newNodesTree)
    DO nodeIdx=1,nodes%numberOfNodes
      nodes%nodes(nodeIdx)%userNumber=userNumbers(nodeIdx)
    ENDDO !nodesIdx
    
    EXITS("Nodes_UserNumbersAllSet")
    RETURN
999 IF(ASSOCIATED(newNodesTree)) CALL Tree_Destroy(newNodesTree,err,error,*998)
998 ERRORSEXITS("Nodes_UserNumbersAllSet",err,error)    
    RETURN 1
   
  END SUBROUTINE Nodes_UserNumbersAllSet
  
  !
  !================================================================================================================================
  !

END MODULE NodeRoutines
