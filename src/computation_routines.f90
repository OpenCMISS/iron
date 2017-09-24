!> \file
!> \author Chris Bradley
!> \brief This module contains all computation routines.
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

!> This module contains all computation routines.

MODULE ComputationRoutines

  USE BaseRoutines
  USE CmissMPI
  USE CmissPetsc
  USE Constants
  USE Kinds
#ifndef NOMPIMOD
  USE MPI
#endif
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE Strings
 
#include "macros.h"

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif

  PRIVATE
  
  !Module parameters  

  !Module types
  
  !>pointer type to ComputationWorkGroupType
  TYPE :: ComputationWorkGroupPtrType
    TYPE(ComputationWorkGroupType), POINTER :: ptr
  END TYPE ComputationWorkGroupPtrType
  
  !>Contains information on logical working groups
  TYPE :: ComputationWorkGroupType
    LOGICAL :: workGroupFinished !<Is .TRUE. if the work group has been finished. .FALSE. if not. 
    TYPE(ComputationWorkGroupType), POINTER:: parent !<Parent of this working groups
    INTEGER(INTG) :: numberOfComputationNodes !<size of the total compurational nodes belonging to this group
    INTEGER(INTG) :: numberOfAvailableRanks !<The number of available ranks for this work group.
    INTEGER(INTG), ALLOCATABLE :: availableRanks(:) !<The list of available ranks for this work group.
    INTEGER(INTG) :: numberOfSubGroups !<The number of sub work groups
    TYPE(ComputationWorkGroupPtrType), ALLOCATABLE:: subGroups(:) !<subGroups(subGroupIdx). A pointer to the subGroupIdx'th sub work group.
    TYPE(ComputationEnvironmentType), POINTER :: computationEnvironment !<pointer to the actual working environment
    LOGICAL :: computationEnvironmentFinished !<Is .TURE. if the actual working environment has been generated, .FALSE. if not
    INTEGER(INTG) :: mpiGroupCommunicator !<The MPI communicator for this work group
    INTEGER(INTG) :: mpiGroup !<The MPI communicator for this work group
    INTEGER(INTG) :: myComputationNodeNumber !<The index of the running process
  END TYPE ComputationWorkGroupType

  !>Contains information on a cache heirarchy
  TYPE ComputationCacheType
    INTEGER(INTG) :: numberOfLevels !<The number of levels in the cache hierarchy
    INTEGER(INTG), ALLOCATABLE :: size(:) !<size(levelIdx). The size of the levelIdx'th cache level.
  END TYPE ComputationCacheType

  !>Contains information on a computation node containing a number of processors
  TYPE ComputationNodeType
    INTEGER(INTG) :: numberOfProcessors !<The number of processors for this computation node
    INTEGER(INTG) :: rank !<The MPI rank of this computation node in the world communicator
    !TYPE(ComputationCacheType) :: CACHE 
    INTEGER(INTG) :: nodeNameLength !<The length of the name of the computation node
    CHARACTER(LEN=MPI_MAX_PROCESSOR_NAME) :: nodeName !<The name of the computation node
  END TYPE ComputationNodeType

  !>Contains information on the MPI type to transfer information about a computation node
  TYPE MPIComputationNodeType
    INTEGER(INTG) :: mpiType !<The MPI data type
    INTEGER(INTG) :: numberOfBlocks !<The number of blocks in the MPI data type. This will be equal to 4.
    INTEGER(INTG) :: blockLengths(4) !<The length of each block.
    INTEGER(INTG) :: types(4) !<The data types of each block.
    INTEGER(MPI_ADDRESS_KIND) :: displacements(4) !<The address displacements to each block.
  END TYPE MPIComputationNodeType

  !>Contains information on the computation environment the program is running in.
  TYPE ComputationEnvironmentType
    INTEGER(INTG) :: mpiVersion !<The version of MPI that we are running with
    INTEGER(INTG) :: mpiSubVersion !<The sub-version of MPI that we are running with
    INTEGER(INTG) :: mpiWorldCommunicator !<The MPI world communicator for OpenCMISS
    INTEGER(INTG) :: mpiCommWorld !<The clone of the MPI world communicator for OpenCMISS
    INTEGER(INTG) :: numberOfWorldComputationNodes !<The number of computation nodes in the world communicator
    INTEGER(INTG) :: myWorldComputationNodeNumber !<The rank of the running process in the world communicator
    TYPE(ComputationNodeType), ALLOCATABLE :: computationNodes(:) !<computationNodes(node_idx). Contains information on the node_idx'th computation node.
    TYPE(MPIComputationNodeType) :: mpiComputationNode !<The MPI data type information to transfer the computation node information.
    TYPE(ComputationWorkGroupType), POINTER :: worldWorkGroup !<A pointer to the work group corresponding to the world communicator
  END TYPE ComputationEnvironmentType

  !Module variables

  LOGICAL, SAVE :: cmissMPIInitialised !<Is .TRUE. if OpenCMISS has initialised MPI
  TYPE(ComputationEnvironmentType), TARGET :: computationEnvironment !<The computation environment the program is running in.

  !Interfaces
  
  ! Access specifiers for subroutines and interfaces(if any)
  PUBLIC ComputationEnvironmentType
  
  PUBLIC ComputationNodeType

  PUBLIC ComputationWorkGroupType
  
  PUBLIC Computation_Initialise,Computation_Finalise

  PUBLIC ComputationEnvironment_WorldCommunicatorGet,ComputationEnvironment_WorldCommunicatorSet
  
  PUBLIC ComputationEnvironment_NodeNumberGet,ComputationEnvironment_NumberOfNodesGet

  PUBLIC Computation_WorkGroupCreateFinish,Computation_WorkGroupCreateStart
  
  PUBLIC Computation_WorkGroupSubGroupAdd
  
  PUBLIC computationEnvironment

CONTAINS

  
  !
  !================================================================================================================================
  !

  !>Finalise a work group and deallocate all memory
  RECURSIVE SUBROUTINE Computation_WorkGroupFinalise(workGroup,err,error,*)

    !Argument Variables
    TYPE(ComputationWorkGroupType),POINTER :: workGroup !<A pointer to the work group to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: subGroupIdx

    ENTERS("Computation_WorkGroupFinalise",err,error,*999)

    IF(ASSOCIATED(workGroup)) THEN
    
      IF(ALLOCATED(workGroup%availableRanks)) DEALLOCATE(workGroup%availableRanks)
      IF(ALLOCATED(workGroup%subGroups)) THEN
        DO subGroupIdx=1,SIZE(workGroup%subGroups,1)
          CALL Computation_WorkGroupFinalise(workGroup%subGroups(subGroupIdx)%ptr,err,error,*999)
        ENDDO !subGroupIdx
        DEALLOCATE(workGroup%subGroups)
      ENDIF
      DEALLOCATE(workGroup)
    ENDIF
    
    EXITS("Computation_WorkGroupFinalise")
    RETURN
999 ERRORSEXITS("Computation_WorkGroupFinalise",err,error)    
    RETURN 1
    
  END SUBROUTINE Computation_WorkGroupFinalise
  
  !
  !================================================================================================================================
  !

  !>Add the work sub-group to the parent group based on the computation requirements (called by user)
  SUBROUTINE Computation_WorkGroupInitialise(workGroup,err,error,*)

    !Argument Variables
    TYPE(ComputationWorkGroupType),POINTER, INTENT(OUT) :: workGroup !<A pointer to the work group to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG):: dummyErr
    TYPE(ComputationWorkgroupPtrType) newWorkGroup
    TYPE(ComputationWorkGroupType),POINTER ::  tmpParentWorkGroup
    TYPE(ComputationWorkgroupPtrType), ALLOCATABLE :: subGroups(:)
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("Computation_WorkGroupInitialise",err,error,*998)

    IF(ASSOCIATED(workGroup)) CALL FlagError("Work group is already associated.",err,error,*998)
    
    ALLOCATE(workGroup,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate work group.",err,error,*999)
    workGroup%workGroupFinished=.FALSE.
    NULLIFY(workGroup%parent)
    workGroup%numberOfComputationNodes=0
    workGroup%numberOfSubGroups=0
    NULLIFY(workGroup%computationEnvironment)
    workGroup%mpiGroupCommunicator=MPI_COMM_NULL
    workGroup%mpiGroup=MPI_GROUP_NULL
    workGroup%myComputationNodeNumber=0
    
    EXITS("Computation_WorkGroupInitialise")
    RETURN
999 CALL Computation_WorkGroupFinalise(workGroup,dummyErr,dummyError,*998)
998 ERRORSEXITS("Computation_WorkGroupIntialise",err,error)    
    RETURN 1
    
  END SUBROUTINE Computation_WorkGroupInitialise
  
  !
  !================================================================================================================================
  !

  !>Add the work sub-group to the parent group based on the computation requirements (called by user)
  SUBROUTINE Computation_WorkGroupSubGroupAdd(parentWorkGroup,numberOfComputationNodes,subWorkGroup,err,error,*)

    !Argument Variables
    TYPE(ComputationWorkGroupType),POINTER, INTENT(INOUT) :: parentWorkGroup
    TYPE(ComputationWorkGroupType),POINTER, INTENT(INOUT) :: subWorkGroup
    INTEGER(INTG),INTENT(IN) :: numberOfComputationNodes
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ComputationWorkgroupPtrType) newWorkGroup
    TYPE(ComputationWorkGroupType),POINTER ::  tmpParentWorkGroup
    TYPE(ComputationWorkgroupPtrType), ALLOCATABLE :: subGroups(:)
    INTEGER(INTG):: I

    ENTERS("Computation_WorkGroupSubGroupAdd",err,error,*999)

    ALLOCATE(newWorkGroup%ptr)
    newWorkGroup%ptr%numberOfComputationNodes = numberOfComputationNodes
    newWorkGroup%ptr%numberOfSubGroups = 0

    IF(ASSOCIATED(parentWorkGroup)) THEN 
      ALLOCATE(subGroups(parentWorkGroup%numberOfSubGroups+1))
      DO I=1,parentWorkGroup%numberOfSubGroups
        subGroups(I)%ptr=>parentWorkGroup%subGroups(I)%ptr
      ENDDO
      !subGroups(1:parentWorkGroup%numberOfSubGroups)=>parentWorkGroup%subGroups(:)
      
      IF(ALLOCATED(parentWorkGroup%subGroups)) THEN 
        DEALLOCATE(parentWorkGroup%subGroups)
      ENDIF
      subGroups(1+parentWorkGroup%numberOfSubGroups)%ptr=>newWorkGroup%ptr
      ALLOCATE(parentWorkGroup%subGroups(SIZE(subGroups,1)))
      DO I=1,SIZE(subGroups,1)
        parentWorkGroup%subGroups(I)%ptr => subGroups(I)%ptr
      ENDDO
      !parentWorkGroup%subGroups(:) => subGroups(:)
      
      DEALLOCATE(subGroups)
      parentWorkGroup%numberOfSubGroups = 1+parentWorkGroup%numberOfSubGroups
      newWorkGroup%ptr%PARENT => parentWorkGroup
      tmpParentWorkGroup => parentWorkGroup 
      DO WHILE(ASSOCIATED(tmpParentWorkGroup)) !Update the computation number of its ancestors
        tmpParentWorkGroup%numberOfComputationNodes = tmpParentWorkGroup%numberOfComputationNodes &
          & + newWorkGroup%ptr%numberOfComputationNodes
        tmpParentWorkGroup => tmpParentWorkGroup%PARENT
      ENDDO
    ELSE !Top level group
      CALL FlagError('parentWorkGroup is not associated, call COMPUTATION_WORK_GROUP_CREATE_START first',&
      & err,error,*999)
    ENDIF
    subWorkGroup => newWorkGroup%ptr

    EXITS("Computation_WorkGroupSubGroupAdd")
    RETURN
999 ERRORSEXITS("Computation_WorkGroupSubGroupAdd",err,error)
    RETURN 1
    
  END SUBROUTINE Computation_WorkGroupSubGroupAdd

  !
  !================================================================================================================================
  !

  !>Start the creation of a work group
  SUBROUTINE Computation_WorkGroupCreateStart(parentWorkGroup,numberOfComputationNodes,workGroup,err,error,*)

    !Argument Variables
    TYPE(ComputationWorkGroupType), POINTER, INTENT(INOUT) :: parentWorkGroup !<The parent work group to create the work group as a sub group for.
    INTEGER(INTG),INTENT(IN) :: numberOfComputationNodes !<The number of nodes to be created in the work group. 
    TYPE(ComputationWorkGroupType), POINTER, INTENT(OUT) :: workGroup !<On exit, the created work group. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Computation_WorkGroupCreateStart",err,error,*999)

    IF(.NOT.ASSOCIATED(parentWorkGroup)) CALL FlagError('Parent work group is not associated.',err,error,*999)
    IF(ASSOCIATED(workGroup)) CALL FlagError("Work group is already associated.",err,error,*999)
    IF(numberOfComputationNodes<1.OR.numberOfComputationNodes>parentWorkGroup%numberOfAvailableRanks) THEN
      localError="The requested number of computation nodes is invalid. The number of computation nodes must be > 0 and <= "// &
        & TRIM(NumberToVString(parentWorkGroup%numberOfAvailableRanks,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    CALL Computation_WorkGroupInitialise(workGroup,err,error,*999)
    
    EXITS("Computation_WorkGroupCreateStart")
    RETURN
999 ERRORSEXITS("Computation_WorkGroupCreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE Computation_WorkGroupCreateStart

  !
  !================================================================================================================================
  !

  !>Generate computation environment for current level work group tree and all it's subgroups recursively 
  RECURSIVE SUBROUTINE Computation_WorkGroupGenerateCompEnviron(workGroup,availableRankList,err,error,*)

    !Argument Variables
    TYPE(ComputationWorkGroupType),POINTER, INTENT(INOUT) :: workGroup
    INTEGER(INTG), ALLOCATABLE, INTENT(INOUT) :: availableRankList(:)
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: rankIdx,mpiIError,rank,originalGroup,newGroup,subGroupIdx
    INTEGER(INTG), ALLOCATABLE :: newAvailableRankList(:)
    
    ENTERS("Computation_WorkGroupGenerateCompEnviron",err,error,*999)
    
    ALLOCATE(workGroup%computationEnvironment)

    !Set size of computation nodes in this communicator
    workGroup%computationEnvironment%numberOfWorldComputationNodes = workGroup%numberOfComputationNodes
    
    !Determine my processes rank
    CALL MPI_COMM_RANK(computationEnvironment%mpiWorldCommunicator,rank,mpiIError)
    CALL MPI_ERROR_CHECK("MPI_COMM_RANK",mpiIError,err,error,*999)
    workGroup%computationEnvironment%myWorldComputationNodeNumber=rank
    
    !Fill in the information for every computation node in this group
    ALLOCATE(workGroup%computationEnvironment%computationNodes(workGroup%computationEnvironment%numberOfWorldComputationNodes))
    IF(SIZE(availableRankList,1)-workGroup%computationEnvironment%numberOfWorldComputationNodes < 0) THEN
      CALL FlagError("NOT ENOUGH RANKS",err,error,*999)
    ENDIF
    DO rankIdx=1,workGroup%computationEnvironment%numberOfWorldComputationNodes,1
      workGroup%computationEnvironment%computationNodes(rankIdx) = &
        & computationEnvironment%computationNodes(availableRankList(rankIdx))
    ENDDO !rankIdx

    !Create a communicator
    CALL MPI_COMM_GROUP(MPI_COMM_WORLD,originalGroup,mpiIError)
    CALL MPI_ERROR_CHECK("MPI_COMM_RANK",mpiIError,err,error,*999)    
    CALL MPI_GROUP_INCL(originalGroup,rankIdx-1,availableRankList(1:rankIdx-1),newGroup,mpiIError)  !Choose the first I-1 ranks
    CALL MPI_ERROR_CHECK("MPI_COMM_RANK",mpiIError,err,error,*999)    
    CALL MPI_COMM_CREATE(MPI_COMM_WORLD,newGroup,workGroup%computationEnvironment%mpiWorldCommunicator,mpiIError)
    CALL MPI_ERROR_CHECK("MPI_COMM_RANK",mpiIError,err,error,*999)    
    CALL MPI_GROUP_FREE(originalGroup,mpiIError)
    CALL MPI_ERROR_CHECK("MPI_COMM_RANK",mpiIError,err,error,*999)    
    CALL MPI_GROUP_FREE(newGroup,mpiIError) 
    CALL MPI_ERROR_CHECK("MPI_COMM_RANK",mpiIError,err,error,*999)    

    !Shrink the availableRankList
    ALLOCATE(newAvailableRankList(SIZE(availableRankList,1)-workGroup%computationEnvironment%numberOfWorldComputationNodes))  
    newAvailableRankList(1:SIZE(newAvailableRankList)) = availableRankList(rankIdx:SIZE(availableRankList,1))
    DEALLOCATE(availableRankList)
    ALLOCATE(availableRankList(SIZE(newAvailableRankList,1)))
    availableRankList(:) = newAvailableRankList(:)

    workGroup%computationEnvironmentFinished = .TRUE.

    !Recursively do this to all its subgroups
    DO subGroupIdx=1,workGroup%numberOfSubGroups,1
      CALL Computation_WorkGroupGenerateCompEnviron(workGroup%subGroups(subGroupIdx)%ptr,&
        & availableRankList,err,error,*999)      
    ENDDO !subGroupIdx

    EXITS("Computation_WorkGroupGenerateCompEnviron")
    RETURN
999 ERRORSEXITS("Computation_WorkGroupGenerateCompEnviron",err,error)
    RETURN 1
    
  END SUBROUTINE Computation_WorkGroupGenerateCompEnviron

  !
  !================================================================================================================================
  !

  !>Generate the hierarchy computation environment based on work group tree
  SUBROUTINE Computation_WorkGroupCreateFinish(worldWorkGroup,err,error,*)

    !Argument Variables
    TYPE(ComputationWorkGroupType),POINTER,INTENT(INOUT) :: worldWorkGroup
    INTEGER(INTG),INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING),INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG),ALLOCATABLE:: availableRankList(:)
    INTEGER(INTG) :: rankIdx,subGroupIdx

    ENTERS("Computation_WorkGroupCreateFinish",err,error,*999)

    !Set the computation environment of the world work group to be the global computation environment
    !(the default communicator in OpenCMISS)
    worldWorkGroup%computationEnvironment => computationEnvironment 
    worldWorkGroup%computationEnvironmentFinished = .TRUE.

    !generate the communicators for subgroups if any
    ALLOCATE(availableRankList(worldWorkGroup%computationEnvironment%numberOfWorldComputationNodes))
    DO rankIdx=0,SIZE(availableRankList,1)-1
      availableRankList(rankIdx+1) = rankIdx
    ENDDO !rankIdx
    DO subGroupIdx=1,worldWorkGroup%numberOfSubGroups,1
      CALL Computation_WorkGroupGenerateCompEnviron(worldWorkGroup%subGroups(subGroupIdx)%ptr,availableRankList, &
        & err,error,*999)
    ENDDO !subGroupIdx

    EXITS("Computation_WorkGroupCreateFinish")
    RETURN
999 ERRORSEXITS("Computation_WorkGroupCreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE Computation_WorkGroupCreateFinish

  !
  !================================================================================================================================
  !

  !>Finalises the computation node data structures and deallocates all memory.
  SUBROUTINE Computation_ComputationNodeFinalise(computationNode,err,error,*)
  
    !Argument Variables
    TYPE(ComputationNodeType),INTENT(INOUT) :: computationNode !<The computation node to finalise
    INTEGER(INTG),INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING),INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Computation_ComputationNodeFinalise",err,error,*999)

    computationNode%numberOfProcessors=0
    computationNode%rank=-1
    computationNode%nodeNameLength=0
    computationNode%nodeName=""    

    EXITS("Computation_ComputationNodeFinalise")
    RETURN
999 ERRORS("Computation_ComputationNodeFinalise",err,error)
    EXITS("Computation_ComputationNodeFinalise")
    RETURN 1
    
  END SUBROUTINE Computation_ComputationNodeFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the computation node data structures.
  SUBROUTINE Computation_ComputationNodeInitialise(computationNode,rank,err,error,*)
  
    !Argument Variables
    TYPE(ComputationNodeType), INTENT(OUT) :: computationNode !<The computation node to initialise
    INTEGER(INTG), INTENT(IN) :: rank !<The MPI rank of the computation node
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: mpiIError

    ENTERS("Computation_ComputationNodeInitialise",err,error,*999)

    computationNode%numberOfProcessors=1
    computationNode%rank=rank
    CALL MPI_GET_PROCESSOR_NAME(computationNode%nodeName,computationNode%nodeNameLength,mpiIError)
    CALL MPI_ERROR_CHECK("MPI_GET_PROCESSOR_NAME",mpiIError,err,error,*999)
    
    EXITS("Computation_ComputationNodeInitialise")
    RETURN
999 ERRORS("Computation_ComputationNodeInitialise",err,error)
    EXITS("Computation_ComputationNodeInitialise")
    RETURN 1
    
  END SUBROUTINE Computation_ComputationNodeInitialise

  !
  !================================================================================================================================
  !

  !>Finalises the MPI computation node type data structure and deallocates all memory.
  SUBROUTINE Computation_MPIComputationNodeFinalise(mpiComputationNode,err,error,*)
  
    !Argument Variables
    TYPE(MPIComputationNodeType) :: mpiComputationNode !<The mpi computation node to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: blockIdx,mpiIError

    ENTERS("Computation_MPIComputationNodeFinalise",err,error,*999)

    DO blockIdx=1,mpiComputationNode%numberOfBlocks
      mpiComputationNode%types(blockIdx)=0
      mpiComputationNode%blockLengths(blockIdx)=0
      mpiComputationNode%displacements(blockIdx)=0
    ENDDO !blockIdx
    mpiComputationNode%numberOfBlocks=0

    IF(mpiComputationNode%mpiType/=MPI_DATATYPE_NULL) THEN
      CALL MPI_TYPE_FREE(mpiComputationNode%mpiType,mpiIError)
      CALL MPI_ERROR_CHECK("MPI_TYPE_FREE",mpiIError,err,error,*999)
    ENDIF

    EXITS("Computation_MPIComputationNodeFinalise")
    RETURN
999 ERRORS("Computation_MPIComputationNodeFinalise",err,error)
    EXITS("Computation_MPIComputationNodeFinalise")
    RETURN 1
    
  END SUBROUTINE Computation_MPIComputationNodeFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the data structure containing the MPI type information for the ComputationNodeType.
  SUBROUTINE Computation_MPIComputationNodeInitialise(computationEnvironment,rank,err,error,*)
  
    !Argument Variables
    TYPE(ComputationEnvironmentType), INTENT(INOUT) :: computationEnvironment !<A pointer to the computation environment to initialise the MPI computation node type for.
    INTEGER(INTG), INTENT(IN) :: rank !<The rank for the initailisation
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: blockIdx,dummyErr,mpiIError
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("Computation_MPIComputationNodeInitialise",err,error,*998)

    !IF(.NOT.ASSOCIATED(computationEnvironment)) CALL FlagError("Computation environment is not associated.",err,error,*998)
    IF(.NOT.ALLOCATED(computationEnvironment%computationNodes)) &
      & CALL FlagError("Computation environment computation nodes is not allocated.",err,error,*999)
    IF(rank<0.OR.rank>computationEnvironment%numberOfWorldComputationNodes) THEN
      localError="The specified rank of "//TRIM(NumberToVString(rank,"*",err,error))// &
        & " is invalid. The rank should be >= 0 and <= "// &
        & TRIM(NumberToVString(computationEnvironment%numberOfWorldComputationNodes,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    computationEnvironment%mpiComputationNode%mpiType=MPI_DATATYPE_NULL
    
    computationEnvironment%mpiComputationNode%numberOfBlocks=4
    computationEnvironment%mpiComputationNode%types=[MPI_INTEGER,MPI_INTEGER,MPI_INTEGER,MPI_CHARACTER]
    computationEnvironment%mpiComputationNode%blockLengths=[1,1,1,MPI_MAX_PROCESSOR_NAME]
		
    CALL MPI_GET_ADDRESS(computationEnvironment%computationNodes(rank)%numberOfProcessors, &
      & computationEnvironment%mpiComputationNode%displacements(1),mpiIError)
    CALL MPI_ERROR_CHECK("MPI_GET_ADDRESS",mpiIError,err,error,*999)
    CALL MPI_GET_ADDRESS(computationEnvironment%computationNodes(rank)%rank, &
      & computationEnvironment%mpiComputationNode%displacements(2),mpiIError)
    CALL MPI_ERROR_CHECK("MPI_GET_ADDRESS",mpiIError,err,error,*999)
    CALL MPI_GET_ADDRESS(computationEnvironment%computationNodes(rank)%nodeNameLength, &
      & computationEnvironment%mpiComputationNode%displacements(3),mpiIError)
    CALL MPI_ERROR_CHECK("MPI_GET_ADDRESS",mpiIError,err,error,*999)
    !CPB 19/02/07 AIX compiler complains about the type of the first parameter i.e., the previous 3 have been integers
    !and this one is not so cast the type.
    CALL MPI_GET_ADDRESS(computationEnvironment%computationNodes(rank)%nodeName, &
      & computationEnvironment%mpiComputationNode%displacements(4),mpiIError)
    CALL MPI_ERROR_CHECK("MPI_GET_ADDRESS",mpiIError,err,error,*999)

    DO blockIdx=4,1,-1
      computationEnvironment%mpiComputationNode%displacements(blockIdx)= &
        & computationEnvironment%mpiComputationNode%displacements(blockIdx)- &
        & computationEnvironment%mpiComputationNode%displacements(1)
    ENDDO !blockIdx

    CALL MPI_TYPE_CREATE_STRUCT(computationEnvironment%mpiComputationNode%numberOfBlocks, &
      & computationEnvironment%mpiComputationNode%blockLengths, &
      & computationEnvironment%mpiComputationNode%displacements, &
      & computationEnvironment%mpiComputationNode%types, &
      & computationEnvironment%mpiComputationNode%mpiType,mpiIError)
    CALL MPI_ERROR_CHECK("MPI_TYPE_CREATE_STRUCT",mpiIError,err,error,*999)

    CALL MPI_TYPE_COMMIT(computationEnvironment%mpiComputationNode%mpiType,mpiIError)
    CALL MPI_ERROR_CHECK("MPI_TYPE_COMMIT",mpiIError,err,error,*999)
    
    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"MPI Computation Node Type Data:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  MPI type = ", &
        & computationEnvironment%mpiComputationNode%mpiType,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of blocks  = ", &
        & computationEnvironment%mpiComputationNode%numberOfBlocks,err,error,*999)
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,computationEnvironment%mpiComputationNode%numberOfBlocks,4,4, &
        & computationEnvironment%mpiComputationNode%types,'("  Block types =",4(X,I15))','(15X,4(X,I15))',err,error,*999)
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,computationEnvironment%mpiComputationNode%numberOfBlocks,8,8, &
        & computationEnvironment%mpiComputationNode%blockLengths,'("  Block lengths =",8(X,I5))','(17X,8(X,I5))',err,error,*999)
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,computationEnvironment%mpiComputationNode%numberOfBlocks,8,8, &
        & computationEnvironment%mpiComputationNode%displacements,'("  Displacements =",8(X,I5))','(17X,8(X,I5))',err,error,*999)
    ENDIF

    EXITS("Computation_MPIComputationNodeInitialise")
    RETURN
999 CALL Computation_MPIComputationNodeFinalise(computationEnvironment%mpiComputationNode,dummyErr,dummyError,*998)
998 ERRORS("Computation_MPIComputationNodeInitialise",err,error)
    EXITS("Computation_MPIComputationNodeInitialise")
    RETURN 1
    
  END SUBROUTINE Computation_MPIComputationNodeInitialise

  !
  !=================================================================================================================================
  !

  !>Finalises the computation environment data structures and deallocates all memory.
  SUBROUTINE Computation_ComputationEnvironmentFinalise(computationEnvironment,err,error,*)

    !Argument Variables
    TYPE(ComputationEnvironmentType) :: computationEnvironment !<A pointer to the computation environment to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: computationNodeIdx,mpiIError
    LOGICAL :: mpiFinalised

    ENTERS("Computation_ComputationEnvironmentFinalise",err,error,*999)

    !IF(ASSOCIATED(computationEnvironment)) THEN
    IF(ALLOCATED(computationEnvironment%computationNodes)) THEN      
      DO computationNodeIdx=0,computationEnvironment%numberOfWorldComputationNodes-1
        CALL Computation_ComputationNodeFinalise(computationEnvironment%computationNodes(computationNodeIdx), &
          & err,error,*999)
      ENDDO !computationNodeIdx
      DEALLOCATE(computationEnvironment%computationNodes)
    ENDIF
    computationEnvironment%numberOfWorldComputationNodes=0
    
    CALL Computation_MPIComputationNodeFinalise(computationEnvironment%mpiComputationNode,err,error,*999)
    
    IF(computationEnvironment%mpiCommWorld /= MPI_COMM_NULL) THEN
      CALL MPI_COMM_FREE(computationEnvironment%mpiCommWorld,mpiIError)
      CALL MPI_ERROR_CHECK("MPI_COMM_FREE",mpiIError,err,error,*999)
    ENDIF
    
    !DEALLOCATE(computationEnvironment)
    !ENDIF

    EXITS("Computation_ComputationEnvironmentFinalise")
    RETURN
999 ERRORS("Computation_ComputationEnvironmentFinalise",err,error)
    EXITS("Computation_ComputationEnvironmentFinalise")
    RETURN 1
    
  END SUBROUTINE Computation_ComputationEnvironmentFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the computation environment data structures.
  SUBROUTINE Computation_ComputationEnvironmentInitialise(computationEnvironment,err,error,*)

    !Argument Variables
    TYPE(ComputationEnvironmentType) :: computationEnvironment !<A pointer to the comptuation environment to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: computationNodeIdx,dummyErr,mpiIError,rank
    LOGICAL :: mpiInitialised
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("Computation_ComputationEnvironmentInitialise",err,error,*999)

    !IF(ASSOCIATED(computationEnvironment)) CALL FlagError("Computation environment is already associated.",err,error,*999)

    !ALLOCATE(computationEnvironment,STAT=err)
    !IF(err/=0) CALL FlagError("Computation environment could not be allocated.",err,error,*999)

    
    computationEnvironment%mpiVersion=0
    computationEnvironment%mpiSubversion=0
    computationEnvironment%mpiWorldCommunicator=MPI_COMM_NULL    
    computationEnvironment%mpiCommWorld=MPI_COMM_NULL
    computationEnvironment%numberOfWorldComputationNodes=0
    computationEnvironment%myWorldComputationNodeNumber=-1
    NULLIFY(computationEnvironment%worldWorkGroup)

    !Get the version of MPI that we are running with
    CALL MPI_GET_VERSION(computationEnvironment%mpiVersion,computationEnvironment%mpiSubversion,mpiIError)
    CALL MPI_ERROR_CHECK("MPI_GET_VERSION",mpiIError,err,error,*999)
    
    !Create a (private) communicator for OpenCMISS as a duplicate MPI_COMM_WORLD
    CALL MPI_COMM_DUP(MPI_COMM_WORLD,computationEnvironment%mpiCommWorld,mpiIError)
    CALL MPI_ERROR_CHECK("MPI_COMM_DUP",mpiIError,err,error,*999)
    !Set the default MPI world communicator to be the duplicate of MPI_COMM_WORLD
    computationEnvironment%mpiWorldCommunicator=computationEnvironment%mpiCommWorld
    
    !Determine the number of ranks/computation nodes we have in our world computation environment
    CALL MPI_COMM_SIZE(computationEnvironment%mpiWorldCommunicator,computationEnvironment%numberOfWorldComputationNodes,mpiIError)
    CALL MPI_ERROR_CHECK("MPI_COMM_SIZE",mpiIError,err,error,*999)

    !Allocate the computation node data structures
    ALLOCATE(computationEnvironment%computationNodes(0:computationEnvironment%numberOfWorldComputationNodes-1),STAT=ERR)
    IF(ERR /=0) CALL FlagError("Could not allocate computational environment computation nodes.",err,error,*999)

    !Determine my processes rank in the world communicator
    CALL MPI_COMM_RANK(computationEnvironment%mpiWorldCommunicator,rank,mpiIError)
    CALL MPI_ERROR_CHECK("MPI_COMM_RANK",mpiIError,err,error,*999)
    computationEnvironment%myWorldComputationNodeNumber=rank
    
#ifdef TAUPROF
    CALL TAU_PROFILE_SET_NODE(rank)
#endif
    
    !Create the MPI type information for the ComputationNodeType
    CALL Computation_MPIComputationNodeInitialise(computationEnvironment,rank,err,error,*999)
    !Fill in all the computation node data structures for this rank at the root position (will be changed later with an
    !allgather call)
    CALL Computation_ComputationNodeInitialise(computationEnvironment%computationNodes(0),rank,err,error,*999)

    !Now transfer all the computation node information to the other computation nodes so that each rank has all the
    !information.
    CALL MPI_ALLGATHER(MPI_IN_PLACE,1,computationEnvironment%mpiComputationNode%mpiType, &
      & computationEnvironment%computationNodes(0),1,computationEnvironment%mpiComputationNode%mpiType, &
      & computationEnvironment%mpiWorldCommunicator,mpiIError)
    CALL MPI_ERROR_CHECK("MPI_ALLGATHER",mpiIError,err,error,*999)
    
    !Setup the world work group.
    !Initialise
    CALL Computation_WorkGroupInitialise(computationEnvironment%worldWorkGroup,err,error,*999)
    !Set the world work group to have all the ranks in the world communicator
    
    IF(diagnostics1) THEN
      !Just let the master node write out this information
      IF(RANK==0) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Computation environment:",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  MPI version = ", &
          & computationEnvironment%mpiVersion,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  MPI subversion = ", &
          & computationEnvironment%mpiSubversion,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of world computation nodes = ", &
          & computationEnvironment%numberOfWorldComputationNodes,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  My world computation node number = ", &
          & computationEnvironment%myWorldComputationNodeNumber,err,error,*999)
        IF(diagnostics2) THEN
          DO computationNodeIdx=0,computationEnvironment%numberOfWorldComputationNodes-1
            CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Computation Node:",err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of Processors = ", &
              & computationEnvironment%computationNodes(computationNodeIdx)%numberOfProcessors,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    MPI rank = ", &
              & computationEnvironment%computationNodes(computationNodeIdx)%rank,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Node Name = ", &
              & computationEnvironment%computationNodes(computationNodeIdx)%nodeName,err,error,*999)
          ENDDO !computationNodeIdx
        ENDIF
      ENDIF
    ENDIF
    
    EXITS("Computation_ComputationEnvironmentInitialise")    
    RETURN
999 CALL Computation_ComputationEnvironmentFinalise(computationEnvironment,dummyErr,dummyError,*998)
998 ERRORS("Computation_ComputationEnvironmentInitialise",err,error)
    EXITS("Computation_ComputationEnvironmentInitialise")
    RETURN 1
    
  END SUBROUTINE Computation_ComputationEnvironmentInitialise

  !
  !=================================================================================================================================
  !

  !>Finalises the computation data structures and deallocates all memory.
  SUBROUTINE Computation_Finalise(err,error,*)

    !Argument Variables
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: mpiIError
    LOGICAL :: mpiFinalised

    ENTERS("Computation_Finalise",err,error,*999)

    CALL Computation_ComputationEnvironmentFinalise(computationEnvironment,err,error,*999)
    
    !Finalise PETSc
    !Call this after MPI_COMM_FREE as PETSc routines are called when some MPI comm attributes are freed.
    !CALL Petsc_LogView(PETSC_COMM_WORLD,"OpenCMISSTest.petsc",err,error,*999)
    CALL Petsc_Finalise(err,error,*999)

    IF(cmissMPIInitialised) THEN
      !Check if MPI has been finalised
      CALL MPI_FINALIZED(mpiFinalised,mpiIError)
      CALL MPI_ERROR_CHECK("MPI_FINALIZED",mpiIError,err,error,*999)
      IF(.NOT.mpiFinalised) THEN
        CALL MPI_FINALIZE(mpiIError)
        CALL MPI_ERROR_CHECK("MPI_FINALIZE",mpiIError,err,error,*999)
      ENDIF
    ENDIF

    EXITS("Computation_Finalise")
    RETURN
999 ERRORSEXITS("Computation_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE Computation_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the computation data structures.
  SUBROUTINE Computation_Initialise(err,error,*)

    !Argument Variables
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: computationNodeIdx,dummyErr,mpiIError,rank
    LOGICAL :: mpiInitialised
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("Computation_Initialise",err,error,*999)
    
    !Check if MPI has been initialised
    cmissMPIInitialised=.FALSE.
    CALL MPI_INITIALIZED(mpiInitialised,mpiIError)
    CALL MPI_ERROR_CHECK("MPI_INITIALIZED",mpiIError,err,error,*999)
    IF(.NOT.mpiInitialised) THEN
      !Initialise the MPI environment
      CALL MPI_INIT(mpiIError)
      CALL MPI_ERROR_CHECK("MPI_INIT",mpiIError,err,error,*999)
      cmissMPIInitialised=.TRUE.
    ENDIF

    CALL Computation_ComputationEnvironmentInitialise(computationEnvironment,err,error,*999)

    !Initialise node numbers in base routines.
    CALL ComputationNodeNumbersSet(computationEnvironment%myWorldComputationNodeNumber,computationEnvironment% &
      & numberOfWorldComputationNodes,err,error,*999)
    
    !Initialise PETSc
    CALL Petsc_Initialise(PETSC_NULL_CHARACTER,err,error,*999)
    
    IF(diagnostics1) THEN
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"OpenCMISS MPI initialised = ",cmissMPIInitialised,err,error,*999)
    ENDIF
    
    EXITS("Computation_Initialise")        
    RETURN
999 CALL Computation_Finalise(dummyErr,dummyError,*998)
998 ERRORSEXITS("Computation_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE Computation_Initialise

  !
  !================================================================================================================================
  !

  !>Sets the world communicator to the given on. Note: This routine should be called straight after the main OpenCMISS initialise
  !>routine. If it is called after objects have started to be setup then good luck!
  SUBROUTINE ComputationEnvironment_WorldCommunicatorSet(worldCommunicator,err,error,*)
    
    !Argument Variables
    INTEGER(INTG), INTENT(IN) :: worldCommunicator !<The world communicator to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: compareResult,mpiIError
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("ComputationEnvironment_WorldCommunicatorSet",err,error,*999)
    
    !Perform a sanity check to see if the communicator is valid
    CALL MPI_COMM_COMPARE(worldCommunicator,computationEnvironment%mpiWorldCommunicator,compareResult,mpiIError)
    CALL MPI_ERROR_CHECK("MPI_COMM_COMPARE",mpiIError,err,error,*999)
    SELECT CASE(compareResult)
    CASE(MPI_IDENT,MPI_CONGRUENT,MPI_SIMILAR)
       !OK to use. 
       !!\TODO We should re-initialise the computation environment.   
       computationEnvironment%mpiWorldCommunicator=worldCommunicator
    CASE(MPI_UNEQUAL)
       !Vastly different to the current communicator. Stop.
       localError="The supplied world communicator of "//TRIM(NumberToVString(worldCommunicator,"*",err,error))// &      
            & " is unequal in structure to the current communicator of "// &
            & TRIM(NumberToVString(computationEnvironment%mpiWorldCommunicator,"*",err,error))//"."
       CALL FlagError(localError,err,error,*999)
    CASE DEFAULT
       localError="The MPI compare result of "//TRIM(NumberToVString(compareResult,"*",err,error))// &
            & " for world communicator "//TRIM(NumberToVString(worldCommunicator,"*",err,error))//" is invalid."
       CALL FlagError(localError,err,error,*999)
    END SELECT
 
    EXITS("ComputationEnvironment_WorldCommunicatorSet")
    RETURN
999 ERRORS("ComputationEnvironment_WorldCommunicatorSet",err,error)
    EXITS("ComputationEnvironment_WorldCommunicatorSet")
    RETURN 1
    
  END SUBROUTINE ComputationEnvironment_WorldCommunicatorSet

  !
  !================================================================================================================================
  !

  !>Gets the current world communicator.
  SUBROUTINE ComputationEnvironment_WorldCommunicatorGet(worldCommunicator,err,error,*)

    !Argument Variables
    INTEGER(INTG), INTENT(OUT) :: worldCommunicator !<On return, the current world communicator
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("ComputationEnvironment_WorldCommunicatorGet",err,error,*999)

    worldCommunicator=computationEnvironment%mpiWorldCommunicator
 
    EXITS("ComputationEnvironment_WorldCommunicatorGet")
    RETURN
999 ERRORS("ComputationEnvironment_WorldCommunicatorGet",err,error)
    EXITS("ComputationEnvironment_WorldCommunicatorGet")
    RETURN 1
    
  END SUBROUTINE ComputationEnvironment_WorldCommunicatorGet

  !
  !================================================================================================================================
  !

  !>Returns the number/rank of the computation nodes.
  FUNCTION ComputationEnvironment_NodeNumberGet(err,error)
      
    !Argument Variables
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    INTEGER(INTG) :: ComputationEnvironment_NodeNumberGet !<On exit the computation node number/rank of the current process.
    !Local Variables

    ENTERS("ComputationEnvironment_NodeNumberGet",err,error,*999)

    IF(ALLOCATED(computationEnvironment%computationNodes)) THEN
      ComputationEnvironment_NodeNumberGet=computationEnvironment%myWorldComputationNodeNumber
    ELSE
      CALL FlagError("Computation environment not initialised",err,error,*999)
    ENDIF
    
    EXITS("ComputationEnvironment_NodeNumberGet")
    RETURN
999 ERRORSEXITS("ComputationEnvironment_NodeNumberGet",err,error)
    RETURN
    
  END FUNCTION ComputationEnvironment_NodeNumberGet

  !
  !================================================================================================================================
  !
  
  !>Returns the number of computation nodes.
  FUNCTION ComputationEnvironment_NumberOfNodesGet(err,error)
     
    !Argument Variables
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    INTEGER(INTG) :: ComputationEnvironment_NumberOfNodesGet !<On exit, the number of computation nodes for the program.
    !Local Variables

    ENTERS("ComputationEnvironment_NumberOfNodesGet",err,error,*999)

    IF(ALLOCATED(computationEnvironment%computationNodes)) THEN
      ComputationEnvironment_NumberOfNodesGet=computationEnvironment%numberOfWorldComputationNodes
    ELSE
      CALL FlagError("Computation environment not initialised",err,error,*999)
    ENDIF
    
    EXITS("ComputationEnvironment_NumberOfNodesGet")
    RETURN
999 ERRORSEXITS("ComputationEnvironment_NumberOfNodesGet",err,error)
    RETURN
    
  END FUNCTION ComputationEnvironment_NumberOfNodesGet

  !
  !================================================================================================================================
  !

END MODULE ComputationRoutines
