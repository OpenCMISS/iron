!> \file
!> \author Chris Bradley
!> \brief This module contains all computational environment variables.
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

!> This module contains all computational environment variables.

MODULE ComputationEnvironment

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
  
  !>pointer type to ComputationalWorkGroupType
  TYPE :: ComputationalWorkGroupPtrType
     TYPE(ComputationalWorkGroupType), POINTER :: ptr
  END TYPE ComputationalWorkGroupPtrType
  
  !>Contains information on logical working groups (added by Robert on 01/04/2010)
  TYPE :: ComputationalWorkGroupType
    INTEGER(INTG) :: numberOfComputationalNodes !<size of the total compurational nodes belonging to this group
    INTEGER(INTG) :: numberOfSubWorkGroups !<size of sub working grous
    TYPE(ComputationalWorkGroupType), POINTER:: PARENT !<Parent of this working groups
    TYPE(ComputationalWorkGroupPtrType), ALLOCATABLE:: SUB_WORK_GROUPS(:) !<non-leaf node: The sub working groups
    
    TYPE(ComputationalEnvironmentType), POINTER :: computationalEnvironment !<pointer to the actual working environment
    LOGICAL :: computationalEnvironmentFinished !<Is .TURE. if the actual working environment has been generated, .FALSE. if not
  END TYPE ComputationalWorkGroupType

  !>Contains information on a cache heirarchy
  TYPE ComputationalCacheType
    INTEGER(INTG) :: NUMBER_LEVELS !<The number of levels in the cache hierarchy
    INTEGER(INTG),ALLOCATABLE :: SIZE(:) !<SIZE(level_idx). The size of the level_idx'th cache level.
  END TYPE ComputationalCacheType

  !>Contains information on a computational node containing a number of processors
  TYPE ComputationalNodeType
    INTEGER(INTG) :: NUMBER_PROCESSORS !<The number of processors for this computational node
    INTEGER(INTG) :: RANK !<The MPI rank of this computational node
   !TYPE(ComputationalCacheType) :: CACHE 
    INTEGER(INTG) :: NODE_NAME_LENGTH !<The length of the name of the computational node
    CHARACTER(LEN=MPI_MAX_PROCESSOR_NAME) :: NODE_NAME !<The name of the computational node
  END TYPE ComputationalNodeType

  !>Contains information on the MPI type to transfer information about a computational node
  TYPE MPIComputationalNodeType
    INTEGER(INTG) :: MPI_TYPE !<The MPI data type
    INTEGER(INTG) :: NUM_BLOCKS !<The number of blocks in the MPI data type. This will be equal to 4.
    INTEGER(INTG) :: BLOCK_LENGTHS(4) !<The length of each block.
    INTEGER(INTG) :: TYPES(4) !<The data types of each block.
    INTEGER(MPI_ADDRESS_KIND) :: DISPLACEMENTS(4) !<The address displacements to each block.
  END TYPE MPIComputationalNodeType

  !>Contains information on the computational environment the program is running in.
  TYPE ComputationalEnvironmentType
    LOGICAL :: cmissMPIInitialised !<Is .TRUE. if OpenCMISS has initialised MPI, .FALSE. if not.
    INTEGER(INTG) :: mpiCommunicator !<The MPI communicator for cmiss
    INTEGER(INTG) :: mpiCommWorld !<The MPI communicator for cmiss
    INTEGER(INTG) :: numberOfComputationalNodes !<The number of computational nodes
    INTEGER(INTG) :: myComputationalNodeNumber !<The index of the running process
    TYPE(ComputationalNodeType), ALLOCATABLE :: computationalNodes(:) !<computationalNodes(node_idx). Contains information on the node_idx'th computational node. 
  END TYPE ComputationalEnvironmentType

  !Module variables

  TYPE(ComputationalEnvironmentType), TARGET :: computationalEnvironment !<The computational environment the program is running in.
  TYPE(MPIComputationalNodeType) :: mpiComputationalNodeTypeData !<The MPI data on the computational nodes.

  !Interfaces

  INTERFACE ComputationalEnvironment_NodeNumberGet
    MODULE PROCEDURE COMPUTATIONAL_NODE_NUMBER_GET
  END INTERFACE ComputationalEnvironment_NodeNumberGet

  INTERFACE ComputationalEnvironment_NumberOfNodesGet
    MODULE PROCEDURE COMPUTATIONAL_NODES_NUMBER_GET
  END INTERFACE ComputationalEnvironment_NumberOfNodesGet

  INTERFACE ComputationalEnvironment_Initialise
    MODULE PROCEDURE COMPUTATIONAL_ENVIRONMENT_INITIALISE
  END INTERFACE ComputationalEnvironment_Initialise

  INTERFACE ComputationalEnvironment_Finalise
    MODULE PROCEDURE COMPUTATIONAL_ENVIRONMENT_FINALISE
  END INTERFACE ComputationalEnvironment_Finalise
  
  ! Access specifiers for subroutines and interfaces(if any)
  PUBLIC ComputationalEnvironmentType
  
  PUBLIC ComputationalNodeType

  PUBLIC ComputationalWorkGroupType
  
  PUBLIC ComputationalEnvironment_Initialise,ComputationalEnvironment_Finalise

  PUBLIC ComputationalEnvironment_WorldCommunicatorGet,ComputationalEnvironment_WorldCommunicatorSet
  
  PUBLIC ComputationalEnvironment_NodeNumberGet,ComputationalEnvironment_NumberOfNodesGet
  
  PUBLIC COMPUTATIONAL_WORK_GROUP_SUBGROUP_ADD, COMPUTATIONAL_WORK_GROUP_CREATE_START, COMPUTATIONAL_WORK_GROUP_CREATE_FINISH
  
  PUBLIC computationalEnvironment

CONTAINS
  !
  !================================================================================================================================
  !

  !>Add the work sub-group to the parent group based on the computational requirements (called by user)
  SUBROUTINE COMPUTATIONAL_WORK_GROUP_SUBGROUP_ADD(PARENT_WORK_GROUP, numberOfComputationalNodes, &
   & ADDED_WORK_GROUP,err,error,*)

    !Argument Variables
    TYPE(ComputationalWorkGroupType),POINTER, INTENT(INOUT) :: PARENT_WORK_GROUP
    TYPE(ComputationalWorkGroupType),POINTER, INTENT(INOUT) :: ADDED_WORK_GROUP
    INTEGER(INTG),INTENT(IN) :: numberOfComputationalNodes
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ComputationalWorkgroupPtrType) NEW_WORK_GROUP
    TYPE(ComputationalWorkGroupType),POINTER ::  TMP_PARENT_WORK_GROUP
    TYPE(ComputationalWorkgroupPtrType), ALLOCATABLE :: SUB_WORK_GROUPS(:)
    INTEGER(INTG):: I

    ENTERS("COMPUTATIONAL_WORK_GROUP_SUBGROUP_ADD",err,error,*999)

    ALLOCATE(NEW_WORK_GROUP%ptr)
    NEW_WORK_GROUP%ptr%numberOfComputationalNodes = numberOfComputationalNodes
    NEW_WORK_GROUP%ptr%numberOfSubWorkGroups = 0

    IF(ASSOCIATED(PARENT_WORK_GROUP)) THEN 
      ALLOCATE(SUB_WORK_GROUPS(PARENT_WORK_GROUP%numberOfSubWorkGroups+1))
      DO I=1,PARENT_WORK_GROUP%numberOfSubWorkGroups
        SUB_WORK_GROUPS(I)%ptr=>PARENT_WORK_GROUP%SUB_WORK_GROUPS(I)%ptr
      ENDDO
      !SUB_WORK_GROUPS(1:PARENT_WORK_GROUP%numberOfSubWorkGroups)=>PARENT_WORK_GROUP%SUB_WORK_GROUPS(:)
      
      IF(ALLOCATED(PARENT_WORK_GROUP%SUB_WORK_GROUPS)) THEN 
        DEALLOCATE(PARENT_WORK_GROUP%SUB_WORK_GROUPS)
      ENDIF
      SUB_WORK_GROUPS(1+PARENT_WORK_GROUP%numberOfSubWorkGroups)%ptr=>NEW_WORK_GROUP%ptr
      ALLOCATE(PARENT_WORK_GROUP%SUB_WORK_GROUPS(SIZE(SUB_WORK_GROUPS,1)))
      DO I=1,SIZE(SUB_WORK_GROUPS,1)
        PARENT_WORK_GROUP%SUB_WORK_GROUPS(I)%ptr => SUB_WORK_GROUPS(I)%ptr
      ENDDO
      !PARENT_WORK_GROUP%SUB_WORK_GROUPS(:) => SUB_WORK_GROUPS(:)
      
      DEALLOCATE(SUB_WORK_GROUPS)
      PARENT_WORK_GROUP%numberOfSubWorkGroups = 1+PARENT_WORK_GROUP%numberOfSubWorkGroups
      NEW_WORK_GROUP%ptr%PARENT => PARENT_WORK_GROUP
      TMP_PARENT_WORK_GROUP => PARENT_WORK_GROUP 
      DO WHILE(ASSOCIATED(TMP_PARENT_WORK_GROUP)) !Update the computational number of its ancestors
        TMP_PARENT_WORK_GROUP%numberOfComputationalNodes = TMP_PARENT_WORK_GROUP%numberOfComputationalNodes &
          & + NEW_WORK_GROUP%ptr%numberOfComputationalNodes
        TMP_PARENT_WORK_GROUP => TMP_PARENT_WORK_GROUP%PARENT
      ENDDO
    ELSE !Top level group
      CALL FlagError('PARENT_WORK_GROUP is not associated, call COMPUTATIONAL_WORK_GROUP_CREATE_START first',&
      & err,error,*999)
    ENDIF
    ADDED_WORK_GROUP => NEW_WORK_GROUP%ptr

    EXITS("COMPUTATIONAL_WORK_GROUP_SUBGROUP_ADD")
    RETURN
999 ERRORSEXITS("COMPUTATIONAL_WORK_GROUP_SUBGROUP_ADD",err,error)
    RETURN 1
  END SUBROUTINE COMPUTATIONAL_WORK_GROUP_SUBGROUP_ADD

  !
  !================================================================================================================================
  !

  !>Create the highest level work group (Default: GROUP_WORLD)
  SUBROUTINE COMPUTATIONAL_WORK_GROUP_CREATE_START(WORLD_WORK_GROUP,numberOfComputationalNodes,err,error,*)

    !Argument Variables
    TYPE(ComputationalWorkGroupType),POINTER, INTENT(INOUT) :: WORLD_WORK_GROUP
    INTEGER(INTG),INTENT(IN) :: numberOfComputationalNodes
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ComputationalWorkgroupPtrType) NEW_WORK_GROUP

    ENTERS("COMPUTATIONAL_WORK_GROUP_CREATE_START",err,error,*999)

    IF(ASSOCIATED(WORLD_WORK_GROUP)) THEN 
      CALL FlagError('WORLD_WORK_GROUP is already associated', ERR, ERROR, *999)
    ELSE
      ALLOCATE(NEW_WORK_GROUP%ptr)
      NEW_WORK_GROUP%ptr%numberOfComputationalNodes = numberOfComputationalNodes
      NEW_WORK_GROUP%ptr%numberOfSubWorkGroups = 0   
      NULLIFY(NEW_WORK_GROUP%ptr%PARENT) !It is the highest level work group already 
      NULLIFY(NEW_WORK_GROUP%ptr%computationalEnvironment) !Generate this later in COMPUTATIONAL_WORK_GROUP_CREATE_FINISH
      WORLD_WORK_GROUP=>NEW_WORK_GROUP%ptr
    ENDIF
    
    EXITS("COMPUTATIONAL_WORK_GROUP_CREATE_START")
    RETURN
999 ERRORSEXITS("COMPUTATIONAL_WORK_GROUP_CREATE_START",err,error)
    RETURN 1
  END SUBROUTINE COMPUTATIONAL_WORK_GROUP_CREATE_START

  !
  !================================================================================================================================
  !

  !>Generate computational environment for current level work group tree and all it's subgroups recursively 
  RECURSIVE SUBROUTINE Computational_WorkGroupGenerateCompEnviron(WORK_GROUP,AVAILABLE_RANK_LIST,err,error,*)

    !Argument Variables
    TYPE(ComputationalWorkGroupType),POINTER, INTENT(INOUT) :: WORK_GROUP
    INTEGER(INTG), ALLOCATABLE, INTENT(INOUT) :: AVAILABLE_RANK_LIST(:)
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: I,MPI_IERROR,RANK,ORIGINAL_GROUP,NEW_GROUP
    INTEGER(INTG), ALLOCATABLE :: NEW_AVAILABLE_RANK_LIST(:)
    
    ENTERS("Computational_WorkGroupGenerateCompEnviron",err,error,*999)
    
    ALLOCATE(WORK_GROUP%computationalEnvironment)

    !Set size of computational nodes in this communicator
    WORK_GROUP%computationalEnvironment%numberOfComputationalNodes = WORK_GROUP%numberOfComputationalNodes
    
    !Determine my processes rank
    CALL MPI_COMM_RANK(computationalEnvironment%mpiCommunicator,RANK,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_COMM_RANK",MPI_IERROR,err,error,*999)
    WORK_GROUP%computationalEnvironment%myComputationalNodeNumber=RANK
    
    !Fill in the information for every computational node in this group
    ALLOCATE(WORK_GROUP%computationalEnvironment%computationalNodes(WORK_GROUP%computationalEnvironment%numberOfComputationalNodes))
    I=SIZE(AVAILABLE_RANK_LIST,1)
    IF(SIZE(AVAILABLE_RANK_LIST,1)-WORK_GROUP%computationalEnvironment%numberOfComputationalNodes < 0) THEN
      CALL FlagError("NOT ENOUGH RANKS", ERR, ERROR, *999)
      GOTO 999
    ENDIF
    DO I=1,WORK_GROUP%computationalEnvironment%numberOfComputationalNodes, 1
      WORK_GROUP%computationalEnvironment%computationalNodes(I) = &
        & computationalEnvironment%computationalNodes(AVAILABLE_RANK_LIST(I))
    ENDDO

    !Create a communicator
    !CALL MPI_COMM_DUP(MPI_COMM_WORLD,WORK_GROUP%computationalEnvironment%mpiCommunicator,MPI_IERROR)
    CALL MPI_COMM_GROUP(MPI_COMM_WORLD,ORIGINAL_GROUP,MPI_IERROR);
    CALL MPI_ERROR_CHECK("MPI_COMM_RANK",MPI_IERROR,err,error,*999)    
    CALL MPI_GROUP_INCL(ORIGINAL_GROUP,I-1,AVAILABLE_RANK_LIST(1:I-1),NEW_GROUP,MPI_IERROR)  !Choose the first I-1 ranks
    CALL MPI_ERROR_CHECK("MPI_COMM_RANK",MPI_IERROR,err,error,*999)    
    CALL MPI_COMM_CREATE(MPI_COMM_WORLD,NEW_GROUP,WORK_GROUP%computationalEnvironment%mpiCommunicator,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_COMM_RANK",MPI_IERROR,err,error,*999)    
    CALL MPI_GROUP_FREE(ORIGINAL_GROUP,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_COMM_RANK",MPI_IERROR,err,error,*999)    
    CALL MPI_GROUP_FREE(NEW_GROUP,MPI_IERROR) 
    CALL MPI_ERROR_CHECK("MPI_COMM_RANK",MPI_IERROR,err,error,*999)    

    !Shrink the AVAILABLE_RANK_LIST
    ALLOCATE(NEW_AVAILABLE_RANK_LIST(SIZE(AVAILABLE_RANK_LIST,1)-WORK_GROUP%computationalEnvironment%numberOfComputationalNodes))  
    NEW_AVAILABLE_RANK_LIST(1:SIZE(NEW_AVAILABLE_RANK_LIST)) = AVAILABLE_RANK_LIST(I:SIZE(AVAILABLE_RANK_LIST,1))
    DEALLOCATE(AVAILABLE_RANK_LIST)
    ALLOCATE(AVAILABLE_RANK_LIST(SIZE(NEW_AVAILABLE_RANK_LIST,1)))
    AVAILABLE_RANK_LIST(:) = NEW_AVAILABLE_RANK_LIST(:)

    WORK_GROUP%computationalEnvironmentFinished = .TRUE.

    !Recursively do this to all its subgroups
    DO I=1,WORK_GROUP%numberOfSubWorkGroups,1
      CALL Computational_WorkGroupGenerateCompEnviron(WORK_GROUP%SUB_WORK_GROUPS(I)%ptr,&
        & AVAILABLE_RANK_LIST,err,error,*999)      
    ENDDO

    EXITS("Computational_WorkGroupGenerateCompEnviron")
    RETURN
999 ERRORSEXITS("Computational_WorkGroupGenerateCompEnviron",err,error)
    RETURN 1
    
  END SUBROUTINE Computational_WorkGroupGenerateCompEnviron

  !
  !================================================================================================================================
  !

  !>Generate the hierarchy computational environment based on work group tree
  SUBROUTINE COMPUTATIONAL_WORK_GROUP_CREATE_FINISH(WORLD_WORK_GROUP,err,error,*)

    !Argument Variables
    TYPE(ComputationalWorkGroupType),POINTER,INTENT(INOUT) :: WORLD_WORK_GROUP
    INTEGER(INTG),INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING),INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG),ALLOCATABLE:: AVAILABLE_RANK_LIST(:)
    INTEGER(INTG) :: I

    ENTERS("COMPUTATIONAL_WORK_GROUP_CREATE_FINISH",err,error,*999)

    !Set the computational environment of the world work group to be the global computational environment
    !(the default communicator in OpenCMISS)
    WORLD_WORK_GROUP%computationalEnvironment => computationalEnvironment 
    WORLD_WORK_GROUP%computationalEnvironmentFinished = .TRUE.

    !generate the communicators for subgroups if any
    ALLOCATE(AVAILABLE_RANK_LIST(WORLD_WORK_GROUP%computationalEnvironment%numberOfComputationalNodes))
    DO I=0,SIZE(AVAILABLE_RANK_LIST,1)-1
      AVAILABLE_RANK_LIST(I+1) = I
    END DO
    DO I=1,WORLD_WORK_GROUP%numberOfSubWorkGroups,1
      CALL Computational_WorkGroupGenerateCompEnviron(WORLD_WORK_GROUP%SUB_WORK_GROUPS(I)%ptr,AVAILABLE_RANK_LIST,err,error,*999)
    END DO

    EXITS("COMPUTATIONAL_WORK_GROUP_CREATE_FINISH")
    RETURN
999 ERRORSEXITS("COMPUTATIONAL_WORK_GROUP_CREATE_FINISH",err,error)
    RETURN 1
  END SUBROUTINE COMPUTATIONAL_WORK_GROUP_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Finalises the computational node data structures and deallocates all memory.
  SUBROUTINE COMPUTATIONAL_NODE_FINALISE(COMPUTATIONAL_NODE,err,error,*)
  
    !Argument Variables
    TYPE(ComputationalNodeType),INTENT(INOUT) :: COMPUTATIONAL_NODE !<The computational node to finalise
    INTEGER(INTG),INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING),INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("COMPUTATIONAL_NODE_FINALISE",err,error,*999)

    COMPUTATIONAL_NODE%NUMBER_PROCESSORS=0
    COMPUTATIONAL_NODE%RANK=-1
    COMPUTATIONAL_NODE%NODE_NAME_LENGTH=0
    COMPUTATIONAL_NODE%NODE_NAME=""    

    EXITS("COMPUTATIONAL_NODE_FINALISE")
    RETURN
999 ERRORSEXITS("COMPUTATIONAL_NODE_FINALISE",err,error)
    RETURN 1
  END SUBROUTINE COMPUTATIONAL_NODE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the computational node data structures.
  SUBROUTINE COMPUTATIONAL_NODE_INITIALISE(COMPUTATIONAL_NODE,RANK,err,error,*)
  
    !Argument Variables
    TYPE(ComputationalNodeType), INTENT(OUT) :: COMPUTATIONAL_NODE !<The computational node to initialise
    INTEGER(INTG), INTENT(IN) :: RANK !<The MPI rank of the computational node
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: MPI_IERROR

    ENTERS("COMPUTATIONAL_NODE_INITIALISE",err,error,*999)

!    COMPUTATIONAL_NODE%NUMBER_PROCESSORS=COMP_DETECT_NUMBER_PROCESSORS(ERR)
!    IF(ERR/=0) GOTO 999
    COMPUTATIONAL_NODE%NUMBER_PROCESSORS=1
    COMPUTATIONAL_NODE%RANK=RANK
    CALL MPI_GET_PROCESSOR_NAME(COMPUTATIONAL_NODE%NODE_NAME,COMPUTATIONAL_NODE%NODE_NAME_LENGTH,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_GET_PROCESSOR_NAME",MPI_IERROR,err,error,*999)
    
    EXITS("COMPUTATIONAL_NODE_INITIALISE")
    RETURN
999 ERRORSEXITS("COMPUTATIONAL_NODE_INITIALISE",err,error)
    RETURN 1
  END SUBROUTINE COMPUTATIONAL_NODE_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the data structure containing the MPI type information for the ComputationalNodeType.
  SUBROUTINE COMPUTATIONAL_NODE_MPI_TYPE_FINALISE(err,error,*)
  
    !Argument Variables
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: i,MPI_IERROR

    ENTERS("COMPUTATIONAL_NODE_MPI_TYPE_FINALISE",err,error,*999)

    DO i=1,mpiComputationalNodeTypeData%NUM_BLOCKS
      mpiComputationalNodeTypeData%TYPES(i)=0
      mpiComputationalNodeTypeData%BLOCK_LENGTHS(i)=0
      mpiComputationalNodeTypeData%DISPLACEMENTS(i)=0
    ENDDO !i
    mpiComputationalNodeTypeData%NUM_BLOCKS=0

    IF(mpiComputationalNodeTypeData%MPI_TYPE/=MPI_DATATYPE_NULL) THEN
      CALL MPI_TYPE_FREE(mpiComputationalNodeTypeData%MPI_TYPE,MPI_IERROR)
      CALL MPI_ERROR_CHECK("MPI_TYPE_FREE",MPI_IERROR,err,error,*999)
    ENDIF

    EXITS("COMPUTATIONAL_NODE_MPI_TYPE_FINALISE")
    RETURN
999 ERRORSEXITS("COMPUTATIONAL_NODE_MPI_TYPE_FINALISE",err,error)
    RETURN 1
  END SUBROUTINE COMPUTATIONAL_NODE_MPI_TYPE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the data structure containing the MPI type information for the ComputationalNodeType.
  SUBROUTINE COMPUTATIONAL_NODE_MPI_TYPE_INITIALISE(COMPUTATIONAL_NODE,err,error,*)
  
    !Argument Variables
    TYPE(ComputationalNodeType), INTENT(IN) :: COMPUTATIONAL_NODE !<The computational node containing the MPI type to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: I,MPI_IERROR

    ENTERS("COMPUTATIONAL_NODE_MPI_TYPE_INITIALISE",err,error,*999)

    mpiComputationalNodeTypeData%MPI_TYPE=MPI_DATATYPE_NULL
    
    mpiComputationalNodeTypeData%NUM_BLOCKS=4
    mpiComputationalNodeTypeData%TYPES=[MPI_INTEGER,MPI_INTEGER,MPI_INTEGER,MPI_CHARACTER]
    mpiComputationalNodeTypeData%BLOCK_LENGTHS=[1,1,1,MPI_MAX_PROCESSOR_NAME]
		
    CALL MPI_GET_ADDRESS(COMPUTATIONAL_NODE%NUMBER_PROCESSORS,mpiComputationalNodeTypeData%DISPLACEMENTS(1),MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_GET_ADDRESS",MPI_IERROR,err,error,*999)
    CALL MPI_GET_ADDRESS(COMPUTATIONAL_NODE%RANK,mpiComputationalNodeTypeData%DISPLACEMENTS(2),MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_GET_ADDRESS",MPI_IERROR,err,error,*999)
    CALL MPI_GET_ADDRESS(COMPUTATIONAL_NODE%NODE_NAME_LENGTH,mpiComputationalNodeTypeData%DISPLACEMENTS(3),MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_GET_ADDRESS",MPI_IERROR,err,error,*999)
    !CPB 19/02/07 AIX compiler complains about the type of the first parameter i.e., the previous 3 have been integers
    !and this one is not so cast the type.
    CALL MPI_GET_ADDRESS(COMPUTATIONAL_NODE%NODE_NAME,mpiComputationalNodeTypeData%DISPLACEMENTS(4),MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_GET_ADDRESS",MPI_IERROR,err,error,*999)

    DO i=4,1,-1
      mpiComputationalNodeTypeData%DISPLACEMENTS(I)=mpiComputationalNodeTypeData%DISPLACEMENTS(I)- &
        & mpiComputationalNodeTypeData%DISPLACEMENTS(1)
    ENDDO !i

    CALL MPI_TYPE_CREATE_STRUCT(mpiComputationalNodeTypeData%NUM_BLOCKS,mpiComputationalNodeTypeData%BLOCK_LENGTHS, &
      & mpiComputationalNodeTypeData%DISPLACEMENTS,mpiComputationalNodeTypeData%TYPES, &
      & mpiComputationalNodeTypeData%MPI_TYPE,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_TYPE_CREATE_STRUCT",MPI_IERROR,err,error,*999)

    CALL MPI_TYPE_COMMIT(mpiComputationalNodeTypeData%MPI_TYPE, MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_TYPE_COMMIT",MPI_IERROR,err,error,*999)
    
    IF(DIAGNOSTICS3) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"MPI Computational Node Type Data:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  MPI type = ",mpiComputationalNodeTypeData%MPI_TYPE,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number blocks  = ",mpiComputationalNodeTypeData%NUM_BLOCKS,err,error,*999)
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,mpiComputationalNodeTypeData%NUM_BLOCKS,4,4, &
        & mpiComputationalNodeTypeData%TYPES,'("  Block types =",4(X,I15))','(15X,4(X,I15))',err,error,*999)
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,mpiComputationalNodeTypeData%NUM_BLOCKS,8,8, &
        & mpiComputationalNodeTypeData%BLOCK_LENGTHS,'("  Block lengths =",8(X,I5))','(17X,8(X,I5))',err,error,*999)
      CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,mpiComputationalNodeTypeData%NUM_BLOCKS,8,8, &
        & mpiComputationalNodeTypeData%DISPLACEMENTS,'("  Displacements =",8(X,I5))','(17X,8(X,I5))',err,error,*999)
    ENDIF

    EXITS("COMPUTATIONAL_NODE_MPI_TYPE_INITIALISE")
    RETURN
999 CALL COMPUTATIONAL_NODE_MPI_TYPE_FINALISE(err,error,*998)
998 ERRORSEXITS("COMPUTATIONAL_NODE_MPI_TYPE_INITIALISE",err,error)
    RETURN 1
  END SUBROUTINE COMPUTATIONAL_NODE_MPI_TYPE_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the computational environment data structures and deallocates all memory.
  SUBROUTINE COMPUTATIONAL_ENVIRONMENT_FINALISE(err,error,*)

    !Argument Variables
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: COMPUTATIONAL_NODE,MPI_IERROR
    LOGICAL :: mpiFinalised

    ENTERS("COMPUTATIONAL_ENVIRONMENT_FINALISE",err,error,*999)

    IF(ALLOCATED(computationalEnvironment%computationalNodes)) THEN
      DO COMPUTATIONAL_NODE=0,computationalEnvironment%numberOfComputationalNodes-1
        CALL COMPUTATIONAL_NODE_FINALISE(computationalEnvironment%computationalNodes(COMPUTATIONAL_NODE),err,error,*999)
      ENDDO
      DEALLOCATE(computationalEnvironment%computationalNodes)
    ENDIF
    computationalEnvironment%numberOfComputationalNodes=0

    CALL COMPUTATIONAL_NODE_MPI_TYPE_FINALISE(err,error,*999)

    CALL MPI_COMM_FREE(computationalEnvironment%mpiCommWorld,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_COMM_FREE",MPI_IERROR,err,error,*999)

    !Finalise PetSc
    !Call this after MPI_COMM_FREE as PETSc routines are called when some
    !MPI comm attributes are freed.
    !CALL Petsc_LogView(PETSC_COMM_WORLD,"OpenCMISSTest.petsc",err,error,*999)
    CALL Petsc_Finalise(err,error,*999)

    IF(computationalEnvironment%cmissMPIInitialised) THEN
      !Check if MPI has been finalised
      CALL MPI_FINALIZED(mpiFinalised,MPI_IERROR)
      CALL MPI_ERROR_CHECK("MPI_FINALIZED",MPI_IERROR,err,error,*999)
      IF(.NOT.mpiFinalised) THEN
        CALL MPI_FINALIZE(MPI_IERROR)
        CALL MPI_ERROR_CHECK("MPI_FINALIZE",MPI_IERROR,err,error,*999)
      ENDIF
    ENDIF

    EXITS("COMPUTATIONAL_ENVIRONMENT_FINALISE")
    RETURN
999 ERRORSEXITS("COMPUTATIONAL_ENVIRONMENT_FINALISE",err,error)
    RETURN 1
    
  END SUBROUTINE COMPUTATIONAL_ENVIRONMENT_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the computational environment data structures.
  SUBROUTINE COMPUTATIONAL_ENVIRONMENT_INITIALISE(err,error,*)

    !Argument Variables
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: i,DUMMY_ERR,MPI_IERROR,RANK
    LOGICAL :: mpiInitialised
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    ENTERS("COMPUTATIONAL_ENVIRONMENT_INITIALISE",err,error,*999)

    computationalEnvironment%cmissMPIInitialised=.FALSE.
    
    !Check if MPI has been initialised
    mpiInitialised=.FALSE.
    CALL MPI_INITIALIZED(mpiInitialised,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_INITIALIZED",MPI_IERROR,err,error,*999)
    IF(.NOT.mpiInitialised) THEN
      !Initialise the MPI environment
      CALL MPI_INIT(MPI_IERROR)
      CALL MPI_ERROR_CHECK("MPI_INIT",MPI_IERROR,err,error,*999)
      computationalEnvironment%cmissMPIInitialised=.TRUE.
    ENDIF
      
    !Create a (private) communicator for OpenCMISS as a duplicate MPI_COMM_WORLD
    CALL MPI_COMM_DUP(MPI_COMM_WORLD,computationalEnvironment%mpiCommWorld,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_COMM_DUP",MPI_IERROR,err,error,*999)
    !Set the default MPI communicator to be the duplicate of MPI_COMM_WORLD
    computationalEnvironment%mpiCommunicator=computationalEnvironment%mpiCommWorld
    
    !Determine the number of ranks/computational nodes we have in our computational environment
    CALL MPI_COMM_SIZE(computationalEnvironment%mpiCommunicator,computationalEnvironment%numberOfComputationalNodes,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_COMM_SIZE",MPI_IERROR,err,error,*999)

    !Allocate the computational node data structures
    ALLOCATE(computationalEnvironment%computationalNodes(0:computationalEnvironment%numberOfComputationalNodes-1),STAT=ERR)
    IF(ERR /=0) CALL FlagError("Could not allocate computational nodes",err,error,*999)

    !Determine my processes rank
    CALL MPI_COMM_RANK(computationalEnvironment%mpiCommunicator,RANK,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_COMM_RANK",MPI_IERROR,err,error,*999)
    computationalEnvironment%myComputationalNodeNumber=RANK
    
#ifdef TAUPROF
    CALL TAU_PROFILE_SET_NODE(rank)
#endif

    !Create the MPI type information for the ComputationalNodeType
    CALL COMPUTATIONAL_NODE_MPI_TYPE_INITIALISE(computationalEnvironment%computationalNodes(RANK),err,error,*999)
    !Fill in all the computational node data structures for this rank at the root position (will be changed later with an
    !allgather call)
    CALL COMPUTATIONAL_NODE_INITIALISE(computationalEnvironment%computationalNodes(0),RANK,err,error,*999)

!     !Now transfer all the computational node information to the other computational nodes so that each rank has all the
!     !information.
! !!    CALL MPI_ALLGATHER(computationalEnvironment%computationalNodes(0),1,mpiComputationalNodeTypeData%MPI_TYPE, &
! !!      & computationalEnvironment%computationalNodes(0),1,mpiComputationalNodeTypeData%MPI_TYPE, &
! !!      & computationalEnvironment%mpiCommunicator,MPI_IERROR)
!     CALL MPI_ALLGATHER(MPI_IN_PLACE,1,mpiComputationalNodeTypeData%MPI_TYPE, &
!       & computationalEnvironment%computationalNodes(0),1,mpiComputationalNodeTypeData%MPI_TYPE, &
!       & computationalEnvironment%mpiCommunicator,MPI_IERROR)
!     CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"    Calling MPI_ERROR_CHECK...",err,error,*999)
!     CALL MPI_ERROR_CHECK("MPI_ALLGATHER",MPI_IERROR,err,error,*999)

    !Initialise node numbers in base routines.
    CALL ComputationalNodeNumbersSet(computationalEnvironment%myComputationalNodeNumber,computationalEnvironment% &
      & numberOfComputationalNodes,err,error,*999)
    
    !Initialise PETSc
    CALL Petsc_Initialise(PETSC_NULL_CHARACTER,err,error,*999)
    
    IF(DIAGNOSTICS1) THEN
      !Just let the master node write out this information
      IF(RANK==0) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Computational environment:",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  OpenCMISS MPI initialised = ", &
          & computationalEnvironment%cmissMPIInitialised,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of computational nodes = ", &
          & computationalEnvironment%numberOfComputationalNodes,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  My computational node number = ", &
          & computationalEnvironment%myComputationalNodeNumber,err,error,*999)
        IF(DIAGNOSTICS2) THEN
          DO i=0,computationalEnvironment%numberOfComputationalNodes-1
            CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Computational Node:",err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of Processors = ", &
              & computationalEnvironment%computationalNodes(i)%NUMBER_PROCESSORS,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    MPI rank = ", &
             & computationalEnvironment%computationalNodes(i)%RANK,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Node Name = ", &
              & computationalEnvironment%computationalNodes(i)%NODE_NAME,err,error,*999)
          ENDDO !i
        ENDIF
      ENDIF
    ENDIF

    EXITS("COMPUTATIONAL_ENVIRONMENT_INITIALISE")
    RETURN
999 CALL COMPUTATIONAL_ENVIRONMENT_FINALISE(DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("COMPUTATIONAL_ENVIRONMENT_INITIALISE",err,error)
    RETURN 1
    
  END SUBROUTINE COMPUTATIONAL_ENVIRONMENT_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets the world communicator to the given on. Note: This routine should be called straight after the main OpenCMISS initialise
  !>routine. If it is called after objects have started to be setup then good luck!
  SUBROUTINE ComputationalEnvironment_WorldCommunicatorSet(worldCommunicator,err,error,*)
    
    !Argument Variables
    INTEGER(INTG), INTENT(IN) :: worldCommunicator !<The world communicator to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: compareResult,mpiIError
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("ComputationalEnvironment_WorldCommunicatorSet",err,error,*999)
    
    !Perform a sanity check to see if the communicator is valid
    CALL MPI_COMM_COMPARE(worldCommunicator,computationalEnvironment%mpiCommunicator,compareResult,mpiIError)
    CALL MPI_ERROR_CHECK("MPI_COMM_COMPARE",mpiIError,err,error,*999)
    SELECT CASE(compareResult)
    CASE(MPI_IDENT,MPI_CONGRUENT,MPI_SIMILAR)
       !OK to use. 
       !!\TODO We should re-initialise the computational environment.   
       computationalEnvironment%mpiCommunicator=worldCommunicator
    CASE(MPI_UNEQUAL)
       !Vastly different to the current communicator. Stop.
       localError="The supplied world communicator of "//TRIM(NumberToVString(worldCommunicator,"*",err,error))// &      
            & " is unequal in structure to the current communicator of "// &
            & TRIM(NumberToVString(computationalEnvironment%mpiCommunicator,"*",err,error))//"."
       CALL FlagError(localError,err,error,*999)
    CASE DEFAULT
       localError="The MPI compare result of "//TRIM(NumberToVString(compareResult,"*",err,error))// &
            & " for world communicator "//TRIM(NumberToVString(worldCommunicator,"*",err,error))//" is invalid."
       CALL FlagError(localError,err,error,*999)
    END SELECT
 
    EXITS("ComputationalEnvironment_WorldCommunicatorSet")
    RETURN
999 ERRORS("ComputationalEnvironment_WorldCommunicatorSet",err,error)
    EXITS("ComputationalEnvironment_WorldCommunicatorSet")
    RETURN 1
    
  END SUBROUTINE ComputationalEnvironment_WorldCommunicatorSet

  !
  !================================================================================================================================
  !

  !>Gets the current world communicator.
  SUBROUTINE ComputationalEnvironment_WorldCommunicatorGet(worldCommunicator,err,error,*)

    !Argument Variables
    INTEGER(INTG), INTENT(OUT) :: worldCommunicator !<On return, the current world communicator
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("ComputationalEnvironment_WorldCommunicatorGet",err,error,*999)

    worldCommunicator=computationalEnvironment%mpiCommunicator
 
    EXITS("ComputationalEnvironment_WorldCommunicatorGet")
    RETURN
999 ERRORS("ComputationalEnvironment_WorldCommunicatorGet",err,error)
    EXITS("ComputationalEnvironment_WorldCommunicatorGet")
    RETURN 1
    
  END SUBROUTINE ComputationalEnvironment_WorldCommunicatorGet

  !
  !================================================================================================================================
  !

  !>Returns the number/rank of the computational nodes.
  FUNCTION COMPUTATIONAL_NODE_NUMBER_GET(err,error)
      
    !Argument Variables
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    INTEGER(INTG) :: COMPUTATIONAL_NODE_NUMBER_GET !<On exit the computational node number/rank of the current process.
    !Local Variables

    ENTERS("COMPUTATIONAL_NODE_NUMBER_GET",err,error,*999)

    IF(ALLOCATED(computationalEnvironment%computationalNodes)) THEN
      COMPUTATIONAL_NODE_NUMBER_GET=computationalEnvironment%myComputationalNodeNumber
    ELSE
      CALL FlagError("Computational environment not initialised",err,error,*999)
    ENDIF
    
    EXITS("COMPUTATIONAL_NODE_NUMBER_GET")
    RETURN
999 ERRORSEXITS("COMPUTATIONAL_NODE_NUMBER_GET",err,error)
    RETURN 
  END FUNCTION COMPUTATIONAL_NODE_NUMBER_GET

  !
  !================================================================================================================================
  !
  
  !>Returns the number of computational nodes.
  FUNCTION COMPUTATIONAL_NODES_NUMBER_GET(err,error)
     
    !Argument Variables
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    INTEGER(INTG) :: COMPUTATIONAL_NODES_NUMBER_GET !<On exit, the number of computational nodes for the program.
    !Local Variables

    ENTERS("COMPUTATIONAL_NODES_NUMBER_GET",err,error,*999)

    IF(ALLOCATED(computationalEnvironment%computationalNodes)) THEN
      COMPUTATIONAL_NODES_NUMBER_GET=computationalEnvironment%numberOfComputationalNodes
    ELSE
      CALL FlagError("Computational environment not initialised",err,error,*999)
    ENDIF
    
    EXITS("COMPUTATIONAL_NODES_NUMBER_GET")
    RETURN
999 ERRORSEXITS("COMPUTATIONAL_NODES_NUMBER_GET",err,error)
    RETURN 
  END FUNCTION COMPUTATIONAL_NODES_NUMBER_GET

  !
  !================================================================================================================================
  !

END MODULE ComputationEnvironment
