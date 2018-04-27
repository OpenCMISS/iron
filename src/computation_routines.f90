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
  USE ComputationAccessRoutines
  USE Constants
  USE Kinds
#ifndef NOMPIMOD
  USE MPI
#endif
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE Sorting
  USE Strings
  USE Types
 
#include "macros.h"

  IMPLICIT NONE

  PRIVATE
  
#ifdef NOMPIMOD
#include "mpif.h"
#endif

  !Module parameters  

  !Module types
  
  !Module variables

  LOGICAL, SAVE :: cmissMPIInitialised !<Is .TRUE. if OpenCMISS has initialised MPI

  !Interfaces

  INTERFACE WorkGroup_LabelSet
    MODULE PROCEDURE WorkGroup_LabelSetC
    MODULE PROCEDURE WorkGroup_LabelSetVS
  END INTERFACE WorkGroup_LabelSet
  
  PUBLIC Computation_Initialise,Computation_Finalise
  
  PUBLIC WorkGroup_CreateFinish,WorkGroup_CreateStart

  PUBLIC WorkGroup_Destroy

  PUBLIC WorkGroup_LabelSet

  PUBLIC WorkGroup_NumberOfGroupNodesSet
  
CONTAINS

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
    CALL MPI_ErrorCheck("MPI_GET_PROCESSOR_NAME",mpiIError,err,error,*999)
    
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
      CALL MPI_ErrorCheck("MPI_TYPE_FREE",mpiIError,err,error,*999)
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
    TYPE(ComputationEnvironmentType), POINTER, INTENT(INOUT) :: computationEnvironment !<A pointer to the computation environment to initialise the MPI computation node type for.
    INTEGER(INTG), INTENT(IN) :: rank !<The rank for the initailisation
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: blockIdx,dummyErr,mpiIError
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("Computation_MPIComputationNodeInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(computationEnvironment)) CALL FlagError("Computation environment is not associated.",err,error,*998)
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
    CALL MPI_ErrorCheck("MPI_GET_ADDRESS",mpiIError,err,error,*999)
    CALL MPI_GET_ADDRESS(computationEnvironment%computationNodes(rank)%rank, &
      & computationEnvironment%mpiComputationNode%displacements(2),mpiIError)
    CALL MPI_ErrorCheck("MPI_GET_ADDRESS",mpiIError,err,error,*999)
    CALL MPI_GET_ADDRESS(computationEnvironment%computationNodes(rank)%nodeNameLength, &
      & computationEnvironment%mpiComputationNode%displacements(3),mpiIError)
    CALL MPI_ErrorCheck("MPI_GET_ADDRESS",mpiIError,err,error,*999)
    !CPB 19/02/07 AIX compiler complains about the type of the first parameter i.e., the previous 3 have been integers
    !and this one is not so cast the type.
    CALL MPI_GET_ADDRESS(computationEnvironment%computationNodes(rank)%nodeName, &
      & computationEnvironment%mpiComputationNode%displacements(4),mpiIError)
    CALL MPI_ErrorCheck("MPI_GET_ADDRESS",mpiIError,err,error,*999)

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
    CALL MPI_ErrorCheck("MPI_TYPE_CREATE_STRUCT",mpiIError,err,error,*999)

    CALL MPI_TYPE_COMMIT(computationEnvironment%mpiComputationNode%mpiType,mpiIError)
    CALL MPI_ErrorCheck("MPI_TYPE_COMMIT",mpiIError,err,error,*999)
    
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

  !>Finalises the computation data structures and deallocates all memory.
  SUBROUTINE Computation_Finalise(computationEnvironment,err,error,*)

    !Argument Variables
    TYPE(ComputationEnvironmentType), POINTER :: computationEnvironment !<A pointer to the computation environment to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: mpiIError
    LOGICAL :: mpiFinalised

    ENTERS("Computation_Finalise",err,error,*999)

    CALL ComputationEnvironment_Finalise(computationEnvironment,err,error,*999)
    
    !Finalise PETSc
    !Call this after MPI_COMM_FREE as PETSc routines are called when some MPI comm attributes are freed.
    !CALL Petsc_LogView(PETSC_COMM_WORLD,"OpenCMISSTest.petsc",err,error,*999)
    CALL Petsc_Finalise(err,error,*999)

    IF(cmissMPIInitialised) THEN
      !Check if MPI has been finalised
      CALL MPI_FINALIZED(mpiFinalised,mpiIError)
      CALL MPI_ErrorCheck("MPI_FINALIZED",mpiIError,err,error,*999)
      IF(.NOT.mpiFinalised) THEN
        CALL MPI_FINALIZE(mpiIError)
        CALL MPI_ErrorCheck("MPI_FINALIZE",mpiIError,err,error,*999)
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
  SUBROUTINE Computation_Initialise(context,err,error,*)

    !Argument Variables
    TYPE(ContextType), POINTER :: context !<A pointer to the context to initialise the computation for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,mpiIError
    LOGICAL :: mpiInitialised
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("Computation_Initialise",err,error,*999)

    IF(.NOT.ASSOCIATED(context)) CALL FlagError("Context is not associated.",err,error,*999)
    
    !Check if MPI has been initialised
    cmissMPIInitialised=.FALSE.
    CALL MPI_INITIALIZED(mpiInitialised,mpiIError)
    CALL MPI_ErrorCheck("MPI_INITIALIZED",mpiIError,err,error,*999)
    IF(.NOT.mpiInitialised) THEN
      !Initialise the MPI environment
      CALL MPI_INIT(mpiIError)
      CALL MPI_ErrorCheck("MPI_INIT",mpiIError,err,error,*999)
      cmissMPIInitialised=.TRUE.
    ENDIF

    CALL ComputationEnvironment_Initialise(context,err,error,*999)

    !Initialise node numbers in base routines.
    CALL ComputationNodeNumbersSet(context%computationEnvironment%myWorldComputationNodeNumber, &
      & context%computationEnvironment%numberOfWorldComputationNodes,err,error,*999)
    
    !Initialise PETSc
    CALL Petsc_Initialise(PETSC_NULL_CHARACTER,err,error,*999)
    
    IF(diagnostics1) THEN
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"OpenCMISS MPI initialised = ",cmissMPIInitialised,err,error,*999)
    ENDIF
    
    EXITS("Computation_Initialise")        
    RETURN
999 CALL Computation_Finalise(context%computationEnvironment,dummyErr,dummyError,*998)
998 ERRORSEXITS("Computation_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE Computation_Initialise

  !
  !=================================================================================================================================
  !

  !>Finalises the computation environment data structures and deallocates all memory.
  SUBROUTINE ComputationEnvironment_Finalise(computationEnvironment,err,error,*)

    !Argument Variables
    TYPE(ComputationEnvironmentType), POINTER :: computationEnvironment !<A pointer to the computation environment to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: computationNodeIdx,mpiIError

    ENTERS("Computation_EnvironmentFinalise",err,error,*999)

    IF(ASSOCIATED(computationEnvironment)) THEN
      IF(ALLOCATED(computationEnvironment%computationNodes)) THEN      
        DO computationNodeIdx=0,computationEnvironment%numberOfWorldComputationNodes-1
          CALL Computation_ComputationNodeFinalise(computationEnvironment%computationNodes(computationNodeIdx), &
            & err,error,*999)
        ENDDO !computationNodeIdx
        DEALLOCATE(computationEnvironment%computationNodes)
      ENDIF
      computationEnvironment%numberOfWorldComputationNodes=0
      
      CALL Computation_MPIComputationNodeFinalise(computationEnvironment%mpiComputationNode,err,error,*999)
      
      IF(computationEnvironment%mpiGroupWorld /= MPI_GROUP_NULL) THEN
        CALL MPI_GROUP_FREE(computationEnvironment%mpiGroupWorld,mpiIError)
        CALL MPI_ErrorCheck("MPI_GROUP_FREE",mpiIError,err,error,*999)
      ENDIF
      IF(computationEnvironment%mpiCommWorld /= MPI_COMM_NULL) THEN
        CALL MPI_COMM_FREE(computationEnvironment%mpiCommWorld,mpiIError)
        CALL MPI_ErrorCheck("MPI_COMM_FREE",mpiIError,err,error,*999)
      ENDIF
    
      DEALLOCATE(computationEnvironment)
    ENDIF

    EXITS("ComputationEnvironment_Finalise")
    RETURN
999 ERRORS("ComputationEnvironment_Finalise",err,error)
    EXITS("ComputationEnvironment_Finalise")
    RETURN 1
    
  END SUBROUTINE ComputationEnvironment_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the computation environment data structures for a context.
  SUBROUTINE ComputationEnvironment_Initialise(context,err,error,*)

    !Argument Variables
    TYPE(ContextType), POINTER :: context !<A pointer to the context to initialise the computation environment for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: computationNodeIdx,dummyErr,mpiIError,rank,rankIdx
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("ComputationEnvironment_Initialise",err,error,*999)

    IF(.NOT.ASSOCIATED(context)) CALL FlagError("Context is not associated.",err,error,*999)
    IF(ASSOCIATED(context%computationEnvironment)) &
      & CALL FlagError("The context computation environment is already associated.",err,error,*999)

    ALLOCATE(context%computationEnvironment,STAT=err)
    IF(err/=0) CALL FlagError("Computation environment could not be allocated.",err,error,*999)
    !Initialise
    context%computationEnvironment%context=>context
    context%computationEnvironment%mpiVersion=0
    context%computationEnvironment%mpiSubversion=0
    context%computationEnvironment%mpiCommWorld=MPI_COMM_NULL
    context%computationEnvironment%mpiGroupWorld=MPI_GROUP_NULL
    context%computationEnvironment%numberOfWorldComputationNodes=0
    context%computationEnvironment%myWorldComputationNodeNumber=-1
    NULLIFY(context%computationEnvironment%worldWorkGroup)

    !Get the version of MPI that we are running with
    CALL MPI_GET_VERSION(context%computationEnvironment%mpiVersion,context%computationEnvironment%mpiSubversion,mpiIError)
    CALL MPI_ErrorCheck("MPI_GET_VERSION",mpiIError,err,error,*999)
    
    !Create a (private) communicator for OpenCMISS as a duplicate MPI_COMM_WORLD
    CALL MPI_COMM_DUP(MPI_COMM_WORLD,context%computationEnvironment%mpiCommWorld,mpiIError)
    CALL MPI_ErrorCheck("MPI_COMM_DUP",mpiIError,err,error,*999)
    !Get the default MPI_COMM_GROUP for the world communicator
    CALL MPI_COMM_GROUP(context%computationEnvironment%mpiCommWorld,context%computationEnvironment%mpiGroupWorld,mpiIError)
    CALL MPI_ErrorCheck("MPI_COMM_GROUP",mpiIError,err,error,*999)
    !Set the default MPI world communicator to be the cloned communicator
    context%computationEnvironment%mpiCommWorld=context%computationEnvironment%mpiCommWorld
    
    !Determine the number of ranks/computation nodes we have in our world computation environment
    CALL MPI_COMM_SIZE(context%computationEnvironment%mpiCommWorld,context%computationEnvironment%numberOfWorldComputationNodes, &
      & mpiIError)
    CALL MPI_ErrorCheck("MPI_COMM_SIZE",mpiIError,err,error,*999)

    !Allocate the computation node data structures
    ALLOCATE(context%computationEnvironment%computationNodes(0:context%computationEnvironment%numberOfWorldComputationNodes-1), &
      & STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate computation environment computation nodes.",err,error,*999)

    !Determine my processes rank in the world communicator
    CALL MPI_COMM_RANK(context%computationEnvironment%mpiCommWorld,rank,mpiIError)
    CALL MPI_ErrorCheck("MPI_COMM_RANK",mpiIError,err,error,*999)
    context%computationEnvironment%myWorldComputationNodeNumber=rank
    
#ifdef TAUPROF
    CALL TAU_PROFILE_SET_NODE(rank)
#endif
    
    !Create the MPI type information for the ComputationNodeType
    CALL Computation_MPIComputationNodeInitialise(context%computationEnvironment,rank,err,error,*999)
    !Fill in all the computation node data structures for this rank at the root position (will be changed later with an
    !allgather call)
    CALL Computation_ComputationNodeInitialise(context%computationEnvironment%computationNodes(0),rank,err,error,*999)

    !Now transfer all the computation node information to the other computation nodes so that each rank has all the
    !information.
    CALL MPI_ALLGATHER(MPI_IN_PLACE,1,context%computationEnvironment%mpiComputationNode%mpiType, &
      & context%computationEnvironment%computationNodes(0),1,context%computationEnvironment%mpiComputationNode%mpiType, &
      & context%computationEnvironment%mpiCommWorld,mpiIError)
    CALL MPI_ErrorCheck("MPI_ALLGATHER",mpiIError,err,error,*999)
    
    !Setup the world work group.
    !Initialise
    CALL WorkGroup_Initialise(context%computationEnvironment%worldWorkGroup,err,error,*999)
    !Set the world work group to have all the ranks in the world communicator
    context%computationEnvironment%worldWorkGroup%userNumber=0
    context%computationEnvironment%worldWorkGroup%label="World Work Group"
    context%computationEnvironment%worldWorkGroup%parentWorkGroup=>context%computationEnvironment%worldWorkGroup
    context%computationEnvironment%worldWorkGroup%numberOfGroupComputationNodes= &
      & context%computationEnvironment%numberOfWorldComputationNodes
    context%computationEnvironment%worldWorkGroup%computationEnvironment=>context%computationEnvironment
    ALLOCATE(context%computationEnvironment%worldWorkGroup%worldRanks( &
      & context%computationEnvironment%numberOfWorldComputationNodes),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate world work group world ranks.",err,error,*999)
    ALLOCATE(context%computationEnvironment%worldWorkGroup%availableRanks( &
      & context%computationEnvironment%numberOfWorldComputationNodes),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate world work group available ranks.",err,error,*999)
    DO rankIdx=1,context%computationEnvironment%numberOfWorldComputationNodes
      context%computationEnvironment%worldWorkGroup%worldRanks(rankIdx)=rankIdx-1
      context%computationEnvironment%worldWorkGroup%availableRanks(rankIdx)=rankIdx-1
    ENDDO !rankIdx
    context%computationEnvironment%worldWorkGroup%numberOfAvailableRanks= &
      & context%computationEnvironment%numberOfWorldComputationNodes
    !Create a new MPI group
    CALL MPI_GROUP_INCL(context%computationEnvironment%mpiGroupWorld, &
      & context%computationEnvironment%worldWorkGroup%numberOfGroupComputationNodes, &
      & context%computationEnvironment%worldWorkGroup%worldRanks, &
      & context%computationEnvironment%worldWorkGroup%mpiGroup,mpiIError)
    CALL MPI_ErrorCheck("MPI_GROUP_INCL",mpiIError,err,error,*999)    
    CALL MPI_COMM_CREATE(context%computationEnvironment%mpiCommWorld, &
      & context%computationEnvironment%worldWorkGroup%mpiGroup, &
      & context%computationEnvironment%worldWorkGroup%mpiGroupCommunicator,mpiIError)
    CALL MPI_ErrorCheck("MPI_COMM_CREATE",mpiIError,err,error,*999)
    !Determine ranks
    CALL MPI_COMM_RANK(context%computationEnvironment%mpiCommWorld,rank,mpiIError)
    CALL MPI_ErrorCheck("MPI_COMM_RANK",mpiIError,err,error,*999)
    context%computationEnvironment%worldWorkGroup%myWorldComputationNodeNumber=rank
    CALL MPI_COMM_RANK(context%computationEnvironment%worldWorkGroup%mpiGroupCommunicator,rank,mpiIError)
    CALL MPI_ErrorCheck("MPI_COMM_RANK",mpiIError,err,error,*999)
    context%computationEnvironment%worldWorkGroup%myGroupComputationNodeNumber=rank
    
    context%computationEnvironment%worldWorkGroup%workGroupFinished=.TRUE.
    
    IF(diagnostics1) THEN
      !Just let the master node write out this information
      IF(rank==0) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Computation environment:",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  MPI version = ", &
          & context%computationEnvironment%mpiVersion,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  MPI subversion = ", &
          & context%computationEnvironment%mpiSubversion,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  MPI world communicator = ", &
          & context%computationEnvironment%mpiCommWorld,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  MPI world group = ", &
          & context%computationEnvironment%mpiGroupWorld,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of world computation nodes = ", &
          & context%computationEnvironment%numberOfWorldComputationNodes,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  My world computation node number = ", &
          & context%computationEnvironment%myWorldComputationNodeNumber,err,error,*999)
        IF(diagnostics2) THEN
          DO computationNodeIdx=0,context%computationEnvironment%numberOfWorldComputationNodes-1
            CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Computation Node:",err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of Processors = ", &
              & context%computationEnvironment%computationNodes(computationNodeIdx)%numberOfProcessors, &
              & err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    MPI rank = ", &
              & context%computationEnvironment%computationNodes(computationNodeIdx)%rank,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Node Name = ", &
              & context%computationEnvironment%computationNodes(computationNodeIdx)%nodeName,err,error,*999)
          ENDDO !computationNodeIdx
        ENDIF
      ENDIF
    ENDIF
    
    EXITS("ComputationEnvironment_Initialise")    
    RETURN
999 CALL ComputationEnvironment_Finalise(context%computationEnvironment,dummyErr,dummyError,*998)
998 ERRORS("ComputationEnvironment_Initialise",err,error)
    EXITS("ComputationEnvironment_Initialise")
    RETURN 1
    
  END SUBROUTINE ComputationEnvironment_Initialise

  !
  !=================================================================================================================================
  !

  !>Start the creation of a work group \see OpenCMISS::Iron::cmfe_WorkGroup_CreateStart
  SUBROUTINE WorkGroup_CreateStart(userNumber,parentWorkGroup,workGroup,err,error,*)

    !Argument Variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the work group.
    TYPE(WorkGroupType), POINTER, INTENT(INOUT) :: parentWorkGroup !<The parent work group to create the work group as a sub group for.
    TYPE(WorkGroupType), POINTER, INTENT(OUT) :: workGroup !<On exit, the created work group. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: subGroupIdx
    TYPE(ComputationEnvironmentType), POINTER :: computationEnvironment
    TYPE(WorkGroupPtrType), ALLOCATABLE :: newSubGroups(:)
    TYPE(VARYING_STRING) :: localError

    ENTERS("WorkGroup_CreateStart",err,error,*999)

    IF(.NOT.ASSOCIATED(parentWorkGroup)) CALL FlagError('Parent work group is not associated.',err,error,*999)
    IF(ASSOCIATED(workGroup)) CALL FlagError("Work group is already associated.",err,error,*999)
    NULLIFY(computationEnvironment)
    CALL WorkGroup_ComputationEnvironmentGet(parentWorkGroup,computationEnvironment,err,error,*999)
    CALL WorkGroup_UserNumberFind(computationEnvironment,userNumber,workGroup,err,error,*999)
    IF(ASSOCIATED(workGroup)) THEN
      localError="The user number of "//TRIM(NumberToVString(userNumber,"*",err,error))//" is already is use."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    CALL WorkGroup_Initialise(workGroup,err,error,*999)
    !Set the defaults
    workGroup%userNumber=userNumber
    workGroup%label="Work Group "//TRIM(NumberToVString(userNumber,"*",err,error))
    workGroup%numberOfGroupComputationNodes=1    
    workGroup%computationEnvironment=>parentWorkGroup%computationEnvironment
    !Add the work group to the list of parent sub groups
    ALLOCATE(newSubGroups(parentWorkGroup%numberOfSubGroups+1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new sub groups.",err,error,*999)
    workGroup%parentWorkGroup=>parentWorkGroup
    DO subGroupIdx=1,parentWorkGroup%numberOfSubGroups
      newSubGroups(subGroupIdx)%ptr=>parentWorkGroup%subGroups(subGroupIdx)%ptr
    ENDDO !subGroupIdx
    newSubGroups(parentWorkGroup%numberOfSubGroups+1)%ptr=>workGroup
    CALL MOVE_ALLOC(newSubGroups,parentWorkGroup%subGroups)
    parentWorkGroup%numberOfSubGroups=parentWorkGroup%numberOfSubGroups+1
    
    EXITS("WorkGroup_CreateStart")
    RETURN
999 ERRORSEXITS("WorkGroup_CreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE WorkGroup_CreateStart

  !
  !================================================================================================================================
  !

  !>Finish the creation of a work group \see OpenCMISS::Iron::cmfe_WorkGroup_CreateFinish
  SUBROUTINE WorkGroup_CreateFinish(workGroup,err,error,*)

    !Argument Variables
    TYPE(WorkGroupType), POINTER, INTENT(INOUT) :: workGroup !<A pointer to the work group to finish the creation of
    INTEGER(INTG),INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING),INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG),ALLOCATABLE:: newAvailableRanks(:)
    INTEGER(INTG) :: groupRank,mpiIError,newNumberOfAvailableRanks,rankIdx,worldRank
    TYPE(ComputationEnvironmentType), POINTER :: computationEnvironment
    TYPE(WorkGroupType), POINTER :: parentWorkGroup
    TYPE(VARYING_STRING) :: localError

    ENTERS("WorkGroup_CreateFinish",err,error,*999)

    IF(.NOT.ASSOCIATED(workGroup)) CALL FlagError("Work group is not associated.",err,error,*999)
    IF(workGroup%workGroupFinished) CALL FlagError("Work group has already been finished.",err,error,*999)
    NULLIFY(parentWorkGroup)
    CALL WorkGroup_ParentWorkGroupGet(workGroup,parentWorkGroup,err,error,*999)
    IF(.NOT.ALLOCATED(parentWorkGroup%availableRanks)) &
      & CALL FlagError("Parent work group available ranks is not allocated.",err,error,*999)
    IF(workGroup%numberOfGroupComputationNodes>SIZE(parentWorkGroup%availableRanks,1)) THEN
      localError="There are insufficient parent work group available ranks. There are "// &
        & TRIM(NumberToVString(SIZE(parentWorkGroup%availableRanks,1),"*",err,error))// &
        & " parent ranks available and "// &
        & TRIM(NumberToVString(workGroup%numberOfGroupComputationNodes,"*",err,error))//" ranks are required."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    !Get the ranks from the list of available ranks of the parent.
    ALLOCATE(workGroup%worldRanks(workGroup%numberOfGroupComputationNodes),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate work group world ranks.",err,error,*999)
    newNumberOfAvailableRanks=parentWorkGroup%numberOfAvailableRanks-workGroup%numberOfGroupComputationNodes
    ALLOCATE(newAvailableRanks(newNumberOfAvailableRanks),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate work group parent new available ranks.",err,error,*999)
    DO rankIdx=1,workGroup%numberOfGroupComputationNodes
      workGroup%worldRanks(rankIdx)=parentWorkGroup%availableRanks(rankIdx)
    ENDDO !rankIdx
    DO rankIdx=1,newNumberOfAvailableRanks
      newAvailableRanks(rankIdx)=parentWorkGroup%availableRanks(rankIdx+workGroup%numberOfGroupComputationNodes)
    ENDDO !rankIdx
    CALL Sorting_HeapSort(newAvailableRanks,err,error,*999)
    CALL MOVE_ALLOC(newAvailableRanks,parentWorkGroup%availableRanks)
    parentWorkGroup%numberOfAvailableRanks=newNumberOfAvailableRanks
    !Set up the available ranks for this work group
    ALLOCATE(workGroup%availableRanks(workGroup%numberOfGroupComputationNodes),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate work group available ranks.",err,error,*999)
    DO rankIdx=1,workGroup%numberOfGroupComputationNodes
      workGroup%availableRanks(rankIdx)=workGroup%worldRanks(rankIdx)
    ENDDO !rankIdx
    workGroup%numberOfAvailableRanks=workGroup%numberOfGroupComputationNodes

    NULLIFY(computationEnvironment)
    CALL WorkGroup_ComputationEnvironmentGet(workGroup,computationEnvironment,err,error,*999)
    !Create a new MPI group
    CALL MPI_GROUP_INCL(workGroup%computationEnvironment%mpiGroupWorld,workGroup%numberOfGroupComputationNodes, &
      & workGroup%worldRanks,workGroup%mpiGroup,mpiIError)
    CALL MPI_ErrorCheck("MPI_GROUP_INCL",mpiIError,err,error,*999)    
    CALL MPI_COMM_CREATE(computationEnvironment%mpiCommWorld,workGroup%mpiGroup,workGroup%mpiGroupCommunicator,mpiIError)
    CALL MPI_ErrorCheck("MPI_COMM_CREATE",mpiIError,err,error,*999)
    
    !Determine my processes rank in the group communicator
    IF(workGroup%mpiGroupCommunicator /= MPI_COMM_NULL) THEN
      CALL MPI_COMM_RANK(workGroup%mpiGroupCommunicator,groupRank,mpiIError)
      CALL MPI_ErrorCheck("MPI_COMM_RANK",mpiIError,err,error,*999)
      workGroup%myGroupComputationNodeNumber=groupRank
    ELSE
      workGroup%myGroupComputationNodeNumber=-1
    ENDIF
    !Determine my process rank in the world communicator
    CALL MPI_COMM_RANK(computationEnvironment%mpiCommWorld,worldRank,mpiIError)
    CALL MPI_ErrorCheck("MPI_COMM_RANK",mpiIError,err,error,*999)
    workGroup%myWorldComputationNodeNumber=worldRank


    workGroup%workGroupFinished=.TRUE.
    
    EXITS("WorkGroup_CreateFinish")
    RETURN
999 IF(ALLOCATED(newAvailableRanks)) DEALLOCATE(newAvailableRanks)
    ERRORSEXITS("WorkGroup_CreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE WorkGroup_CreateFinish

  !
  !=================================================================================================================================
  !

  !>Destroy a work group \see OpenCMISS::Iron::cmfe_WorkGroup_Destroy
  SUBROUTINE WorkGroup_Destroy(workGroup,err,error,*)

    !Argument Variables
    TYPE(WorkGroupType), POINTER, INTENT(INOUT) :: workGroup !<The work group to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: userNumber
    TYPE(ComputationEnvironmentType), POINTER :: computationEnvironment
 
    ENTERS("WorkGroup_Destroy",err,error,*999)

    IF(.NOT.ASSOCIATED(workGroup)) CALL FlagError("Work group is not associated.",err,error,*999)
    IF(workGroup%userNumber==0) CALL FlagError("Cannot destroy the world work group.",err,error,*999)

    NULLIFY(computationEnvironment)
    CALL WorkGroup_ComputationEnvironmentGet(workGroup,computationEnvironment,err,error,*999)
    userNumber=workGroup%userNumber
    CALL WorkGroup_DestroyNumber(computationEnvironment,userNumber,err,error,*999)
    
    EXITS("WorkGroup_Destroy")
    RETURN
999 ERRORSEXITS("WorkGroup_Destroy",err,error)
    RETURN 1
    
  END SUBROUTINE WorkGroup_Destroy

  !
  !=================================================================================================================================
  !

  !>Destroy a work group given by a user number and all sub groups under it
  RECURSIVE SUBROUTINE WorkGroup_DestroyNumber(computationEnvironment,workGroupUserNumber,err,error,*)

    !Argument Variables
    TYPE(ComputationEnvironmentType), POINTER :: computationEnvironment !<The computation environment in which to destroy the work group
    INTEGER(INTG), INTENT(IN) :: workGroupUserNumber !<The work group user number to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: count,subGroupIdx,rankIdx
    INTEGER(INTG), ALLOCATABLE :: newAvailableRanks(:)
    TYPE(WorkGroupType), POINTER :: parentWorkGroup,subGroup,workGroup,workGroup2
    TYPE(WorkGroupPtrType), ALLOCATABLE :: newSubGroups(:)
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("WorkGroup_DestroyNumber",err,error,*999)

    IF(.NOT.ASSOCIATED(computationEnvironment)) CALL FlagError("Computation environment is not associated.",err,error,*999)

    NULLIFY(workGroup)
    CALL WorkGroup_UserNumberFind(computationEnvironment,workGroupUserNumber,workGroup,err,error,*999)
    IF(.NOT.ASSOCIATED(workGroup)) THEN
      localError="A work group with a user number of "//TRIM(NumberToVString(workGroupUserNumber,"*",err,error))// &
        & " does not exist."
      CALL FlagError(localError,err,error,*999)
    ENDIF

!!NOTE: We have to find a pointer to the work group to destroy within this routine rather than passing in a pointer to a
!!WorkGroup_DestroyPtr routine because we need to change workGroup%subGroups of the parent work group and this would violate section
!!12.4.1.6 of the Fortran standard if the dummy workGroup pointer argument was associated with the subGroups(x)%ptr actual
!!argument.
    
    IF(workGroup%numberOfSubGroups==0) THEN
      !No work sub groups so delete this instance
      NULLIFY(parentWorkGroup)
      CALL WorkGroup_ParentWorkGroupGet(workGroup,parentWorkGroup,err,error,*999)
      IF(parentWorkGroup%numberOfSubGroups>1) THEN
        !If the parentWorkGroup has more than one sub groups then remove this work group from the list of sub groups.
        ALLOCATE(newSubGroups(parentWorkGroup%numberOfSubGroups-1),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate new sub groups.",err,error,*999)
        count=0
        DO subGroupIdx=1,parentWorkGroup%numberOfSubGroups
          NULLIFY(subGroup)
          CALL WorkGroup_WorkSubGroupGet(parentWorkGroup,subGroupIdx,subGroup,err,error,*999)
          IF(subGroup%userNumber/=workGroup%userNumber) THEN
            count=count+1
            newSubGroups(count)%ptr=>parentWorkGroup%subGroups(subGroupIdx)%ptr
          ENDIF
        ENDDO !subGroupIdx
        CALL MOVE_ALLOC(newSubGroups,parentWorkGroup%subGroups)
        parentWorkGroup%numberOfSubGroups=parentWorkGroup%numberOfSubGroups-1
      ELSE
        IF(ALLOCATED(parentWorkGroup%subGroups)) DEALLOCATE(parentWorkGroup%subGroups)
        parentWorkGroup%numberOfSubGroups=0
      ENDIF
      !Put the work group ranks back into the parent work group available ranks
      ALLOCATE(newAvailableRanks(parentWorkGroup%numberOfAvailableRanks+workGroup%numberOfGroupComputationNodes),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new available ranks.",err,error,*999)
      DO rankIdx=1,parentWorkGroup%numberOfAvailableRanks
        newAvailableRanks(rankIdx)=parentWorkGroup%availableRanks(rankIdx)
      ENDDO !rankIdx
      DO rankIdx=1,workGroup%numberOfGroupComputationNodes
        newAvailableRanks(parentWorkGroup%numberOfAvailableRanks+rankIdx)=workGroup%worldRanks(rankIdx)
      ENDDO !rankIdx
      CALL Sorting_HeapSort(newAvailableRanks,err,error,*999)
      CALL MOVE_ALLOC(newAvailableRanks,parentWorkGroup%availableRanks)
      parentWorkGroup%numberOfAvailableRanks=parentWorkGroup%numberOfAvailableRanks+workGroup%numberOfGroupComputationNodes
      !Finalise the work group
      CALL WorkGroup_Finalise(workGroup,err,error,*999)
    ELSE
      !Recursively delete the sub groups first
      DO WHILE(workGroup%numberOfSubGroups>0)
        NULLIFY(workGroup2)
        CALL WorkGroup_WorkSubGroupGet(workGroup,1,workGroup2,err,error,*999)
        CALL WorkGroup_DestroyNumber(computationEnvironment,workGroup2%userNumber,err,error,*999)
      ENDDO
      !Now delete this instance
      CALL WorkGroup_DestroyNumber(computationEnvironment,workGroup%userNumber,err,error,*999)     
    ENDIF

    EXITS("WorkGroup_DestroyNumber")
    RETURN
999 IF(ALLOCATED(newSubGroups)) DEALLOCATE(newSubGroups)
    IF(ALLOCATED(newAvailableRanks)) DEALLOCATE(newAvailableRanks)
    ERRORSEXITS("WorkGroup_DestroyNumber",err,error)
    RETURN 1
    
  END SUBROUTINE WorkGroup_DestroyNumber

  !
  !=================================================================================================================================
  !

  !>Finalise a work group and deallocate all memory
  RECURSIVE SUBROUTINE WorkGroup_Finalise(workGroup,err,error,*)

    !Argument Variables
    TYPE(WorkGroupType),POINTER :: workGroup !<A pointer to the work group to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: mpiIError,subGroupIdx

    ENTERS("WorkGroup_Finalise",err,error,*999)

    IF(ASSOCIATED(workGroup)) THEN
      workGroup%label=""
      IF(ALLOCATED(workGroup%availableRanks)) DEALLOCATE(workGroup%availableRanks)
      IF(ALLOCATED(workGroup%subGroups)) THEN
        DO subGroupIdx=1,SIZE(workGroup%subGroups,1)
          CALL WorkGroup_Finalise(workGroup%subGroups(subGroupIdx)%ptr,err,error,*999)
        ENDDO !subGroupIdx
        DEALLOCATE(workGroup%subGroups)
      ENDIF
      IF(workGroup%mpiGroupCommunicator/=MPI_COMM_NULL) THEN
        CALL MPI_COMM_FREE(workGroup%mpiGroupCommunicator,mpiIError) 
        CALL MPI_ErrorCheck("MPI_COMM_FREE",mpiIError,err,error,*999)
      ENDIF
      IF(workGroup%mpiGroup/=MPI_GROUP_NULL) THEN
        CALL MPI_GROUP_FREE(workGroup%mpiGroup,mpiIError) 
        CALL MPI_ErrorCheck("MPI_GROUP_FREE",mpiIError,err,error,*999)
      ENDIF
      DEALLOCATE(workGroup)
    ENDIF
    
    EXITS("WorkGroup_Finalise")
    RETURN
999 ERRORSEXITS("WorkGroup_Finalise",err,error)    
    RETURN 1
    
  END SUBROUTINE WorkGroup_Finalise
  
  !
  !=================================================================================================================================
  !

  !>Add the work sub-group to the parent group based on the computation requirements (called by user)
  SUBROUTINE WorkGroup_Initialise(workGroup,err,error,*)

    !Argument Variables
    TYPE(WorkGroupType),POINTER, INTENT(OUT) :: workGroup !<A pointer to the work group to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG):: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("WorkGroup_Initialise",err,error,*998)

    IF(ASSOCIATED(workGroup)) CALL FlagError("Work group is already associated.",err,error,*998)
    
    ALLOCATE(workGroup,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate work group.",err,error,*999)
    workGroup%userNumber=0
    workGroup%workGroupFinished=.FALSE.    
    NULLIFY(workGroup%parentWorkGroup)
    workGroup%label=""
    workGroup%numberOfGroupComputationNodes=0
    workGroup%numberOfAvailableRanks=0
    workGroup%numberOfSubGroups=0
    NULLIFY(workGroup%computationEnvironment)
    workGroup%mpiGroupCommunicator=MPI_COMM_NULL
    workGroup%mpiGroup=MPI_GROUP_NULL
    workGroup%myGroupComputationNodeNumber=0
    workGroup%myWorldComputationNodeNumber=0
    
    EXITS("WorkGroup_Initialise")
    RETURN
999 CALL WorkGroup_Finalise(workGroup,dummyErr,dummyError,*998)
998 ERRORSEXITS("WorkGroup_Initialise",err,error)    
    RETURN 1
    
  END SUBROUTINE WorkGroup_Initialise
  
  !
  !=================================================================================================================================
  !

  !>Set the character label of a work group \see OpenCMISS::Iron::cmfe_WorkGroup_LabelSet
  SUBROUTINE WorkGroup_LabelSetC(workGroup,label,err,error,*)

    !Argument Variables
    TYPE(WorkGroupType), POINTER, INTENT(IN) :: workGroup !<The work group to set the label for
    CHARACTER(LEN=*), INTENT(IN) :: label !<The label to set for the work group
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("WorkGroup_LabelSetC",err,error,*999)

    IF(.NOT.ASSOCIATED(workGroup)) CALL FlagError("Work group is not associated.",err,error,*999)

    workGroup%label=label
    
    EXITS("WorkGroup_LabelSetC")
    RETURN
999 ERRORSEXITS("WorkGroup_LabelSetC",err,error)
    RETURN 1
    
  END SUBROUTINE WorkGroup_LabelSetC

  !
  !=================================================================================================================================
  !

  !>Set the varying string label of a work group \see OpenCMISS::Iron::cmfe_WorkGroup_LabelSet
  SUBROUTINE WorkGroup_LabelSetVS(workGroup,label,err,error,*)

    !Argument Variables
    TYPE(WorkGroupType), POINTER, INTENT(IN) :: workGroup !<The work group to set the label for
    TYPE(VARYING_STRING), INTENT(IN) :: label !<The label to set for the work group
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("WorkGroup_LabelSetVS",err,error,*999)

    IF(.NOT.ASSOCIATED(workGroup)) CALL FlagError("Work group is not associated.",err,error,*999)

    workGroup%label=label
    
    EXITS("WorkGroup_LabelSetVS")
    RETURN
999 ERRORSEXITS("WorkGroup_LabelSetVS",err,error)
    RETURN 1
    
  END SUBROUTINE WorkGroup_LabelSetVS

  !
  !=================================================================================================================================
  !

  !>Set the number of group nodes in a work group \see OpenCMISS::Iron::cmfe_WorkGroup_NumberOfGroupNodesSet
  SUBROUTINE WorkGroup_NumberOfGroupNodesSet(workGroup,numberOfGroupComputationNodes,err,error,*)

    !Argument Variables
    TYPE(WorkGroupType), POINTER, INTENT(OUT) :: workGroup !<The work group to set the number of group nodes for
    INTEGER(INTG), INTENT(IN) :: numberOfGroupComputationNodes !<The number of group nodes to set for the work group
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(WorkGroupType), POINTER :: parentWorkGroup
    TYPE(VARYING_STRING) :: localError

    ENTERS("WorkGroup_NumberOfGroupNodesSet",err,error,*999)

    IF(.NOT.ASSOCIATED(workGroup)) CALL FlagError("Work group is not associated.",err,error,*999)
    IF(workGroup%workGroupFinished) CALL FlagError("Work group has already been finished.",err,error,*999)
    NULLIFY(parentWorkGroup)
    CALL WorkGroup_ParentWorkGroupGet(workGroup,parentWorkGroup,err,error,*999)
    IF(numberOfGroupComputationNodes<1.OR.numberOfGroupComputationNodes>parentWorkGroup%numberOfAvailableRanks) THEN
      localError="The number of group nodes of "//TRIM(NumberToVString(numberOfGroupComputationNodes,"*",err,error))// &
        & " is invalid. The number of group nodes should be > 0 and <= "// &
        & TRIM(NumberToVString(parentWorkGroup%numberOfAvailableRanks,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    workGroup%numberOfGroupComputationNodes=numberOfGroupComputationNodes
    
    EXITS("WorkGroup_NumberOfGroupNodesSet")
    RETURN
999 ERRORSEXITS("WorkGroup_NumberOfGroupNodesSet",err,error)
    RETURN 1
    
  END SUBROUTINE WorkGroup_NumberOfGroupNodesSet

  !
  !=================================================================================================================================
  !

END MODULE ComputationRoutines
