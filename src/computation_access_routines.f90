!> \file
!> \author Chris Bradley
!> \brief This module contains all computation access method routines.
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

!> This module contains all computation access method routines.
MODULE ComputationAccessRoutines
  
  USE BaseRoutines
  USE ISO_VARYING_STRING
  USE Kinds
#ifndef NOMPIMOD
  USE MPI
#endif
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
 
  !Interfaces

  INTERFACE WorkGroup_LabelGet
    MODULE PROCEDURE WorkGroup_LabelGetC
    MODULE PROCEDURE WorkGroup_LabelGetVS
  END INTERFACE WorkGroup_LabelGet

  PUBLIC ComputationEnvironment_NumberOfWorldNodesGet

  PUBLIC ComputationEnvironment_WorldCommunicatorGet

  PUBLIC ComputationEnvironment_WorldNodeNumberGet

  PUBLIC ComputationEnvironment_WorldWorkGroupGet

  PUBLIC WorkGroup_ComputationEnvironmentGet

  PUBLIC WorkGroup_ContextGet

  PUBLIC WorkGroup_Get

  PUBLIC WorkGroup_GroupCommunicatorGet

  PUBLIC WorkGroup_GroupNodeNumberGet

  PUBLIC WorkGroup_LabelGet
  
  PUBLIC WorkGroup_NumberOfGroupNodesGet

  PUBLIC WorkGroup_ParentWorkGroupGet

  PUBLIC WorkGroup_UserNumberFind

  PUBLIC WorkGroup_UserNumberGet

  PUBLIC WorkGroup_WorkSubGroupGet

CONTAINS

  !
  !================================================================================================================================
  !

  !>Gets the current world communicator.
  SUBROUTINE ComputationEnvironment_WorldCommunicatorGet(computationEnvironment,worldCommunicator,err,error,*)

    !Argument Variables
    TYPE(ComputationEnvironmentType), POINTER, INTENT(IN) :: computationEnvironment !<The computational environment to get the world node number for.
    INTEGER(INTG), INTENT(OUT) :: worldCommunicator !<On return, the current world communicator
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("ComputationEnvironment_WorldCommunicatorGet",err,error,*999)

    IF(.NOT.ASSOCIATED(computationEnvironment)) CALL FlagError("Computation environment is not associated.",err,error,*999)
    
    worldCommunicator=computationEnvironment%mpiCommWorld
 
    EXITS("ComputationEnvironment_WorldCommunicatorGet")
    RETURN
999 ERRORS("ComputationEnvironment_WorldCommunicatorGet",err,error)
    EXITS("ComputationEnvironment_WorldCommunicatorGet")
    RETURN 1
    
  END SUBROUTINE ComputationEnvironment_WorldCommunicatorGet

  !
  !================================================================================================================================
  !
  
  !>Returns the number/rank of the computation node in the world communicator  
  SUBROUTINE ComputationEnvironment_WorldNodeNumberGet(computationEnvironment,worldNodeNumber,err,error,*)
      
    !Argument Variables
    TYPE(ComputationEnvironmentType), POINTER, INTENT(IN) :: computationEnvironment !<The computational environment to get the world node number for.
    INTEGER(INTG), INTENT(OUT) :: worldNodeNumber !<On return, the node number in the world communicator.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("ComputationEnvironment_WorldNodeNumberGet",err,error,*999)

    IF(.NOT.ASSOCIATED(computationEnvironment)) CALL FlagError("Computation environment is not associated.",err,error,*999)
    
    worldNodeNumber=computationEnvironment%myWorldComputationNodeNumber
        
    EXITS("ComputationEnvironment_WorldNodeNumberGet")
    RETURN
999 ERRORS("ComputationEnvironment_WorldNodeNumberGet",err,error)
    EXITS("ComputationEnvironment_WorldNodeNumberGet")
    RETURN
    
  END SUBROUTINE ComputationEnvironment_WorldNodeNumberGet

  !
  !================================================================================================================================
  !
  
  !>Gets the number of computation nodes in the world communicator.
  SUBROUTINE ComputationEnvironment_NumberOfWorldNodesGet(computationEnvironment,numberOfWorldNodes,err,error,*)
     
    !Argument Variables
    TYPE(ComputationEnvironmentType), POINTER, INTENT(IN) :: computationEnvironment !<The computational environment to get the world number of nodes for.
    INTEGER(INTG), INTENT(OUT) :: numberOfWorldNodes !<On return, the number of nodes in the world communicator.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("ComputationEnvironment_NumberOfWorldNodesGet",err,error,*999)

    IF(.NOT.ASSOCIATED(computationEnvironment)) CALL FlagError("Computation environment is not associated.",err,error,*999)
    
    numberOfWorldNodes=computationEnvironment%numberOfWorldComputationNodes    
    
    EXITS("ComputationEnvironment_NumberOfWorldNodesGet")
    RETURN
999 ERRORS("ComputationEnvironment_NumberOfWorldNodesGet",err,error)
    EXITS("ComputationEnvironment_NumberOfWorldNodesGet")
    RETURN
    
  END SUBROUTINE ComputationEnvironment_NumberOfWorldNodesGet

  !
  !================================================================================================================================
  !

  !>Gets the world work group from a computational environment.
  SUBROUTINE ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup,err,error,*)

    !Argument variables
    TYPE(ComputationEnvironmentType), POINTER, INTENT(IN) :: computationEnvironment !<The computational environment to get the world work group for.
    TYPE(WorkGroupType), POINTER, INTENT(OUT) :: worldWorkGroup !<On exit, a pointer to the world work group for the computation environment. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("ComputationEnvironment_WorldWorkGroupGet",err,error,*999)

    !Check input arguments
    IF(ASSOCIATED(worldWorkGroup)) CALL FlagError("World work group is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(computationEnvironment)) CALL FlagError("Computation environment is not associated.",err,error,*999)

    !Get the world work group
    worldWorkGroup=>computationEnvironment%worldWorkGroup
    !Check world work group is associated.
    IF(.NOT.ASSOCIATED(worldWorkGroup)) &
      & CALL FlagError("World work group is not associated for the computation environment.",err,error,*999)
    
    EXITS("ComputationEnvironment_WorldWorkGroupGet")
    RETURN
999 NULLIFY(worldWorkGroup)
998 ERRORSEXITS("ComputationEnvironment_WorldWorkGroupGet",err,error)
    RETURN 1
    
  END SUBROUTINE ComputationEnvironment_WorldWorkGroupGet

  !
  !================================================================================================================================
  !

  !>Returns the computation environment associated with a work group.
  SUBROUTINE WorkGroup_ComputationEnvironmentGet(workGroup,computationEnvironment,err,error,*)

    !Argument variables
    TYPE(WorkGroupType), POINTER, INTENT(IN) :: workGroup !<A pointer to the work group to get the computation environment for.
    TYPE(ComputationEnvironmentType), POINTER, INTENT(OUT) :: computationEnvironment !<On return, the computational environment for the work group. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("WorkGroup_ComputationEnvironmentGet",err,error,*998)

    IF(ASSOCIATED(computationEnvironment)) CALL FlagError("Computation environment is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(workGroup)) CALL FlagError("Work group is not associated.",err,error,*999)

    computationEnvironment=>workGroup%computationEnvironment
    IF(.NOT.ASSOCIATED(computationEnvironment)) THEN      
      localError="The computation environment is not associated for work group number "// &
        & TRIM(NumberToVString(workGroup%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("WorkGroup_ComputationEnvironmentGet")
    RETURN
999 NULLIFY(computationEnvironment)
998 ERRORSEXITS("WorkGroup_ComputationEnvironmentGet",err,error)
    RETURN 1
    
  END SUBROUTINE WorkGroup_ComputationEnvironmentGet

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the context for a work group.
  SUBROUTINE WorkGroup_ContextGet(workGroup,context,err,error,*)

    !Argument variables
    TYPE(WorkGroupType), POINTER :: workGroup !<A pointer to the work group to get the context for
    TYPE(ContextType), POINTER :: context !<On exit, a pointer to the context for the work group. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("WorkGroup_ContextGet",err,error,*998)

    IF(ASSOCIATED(context)) CALL FlagError("Context is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(workGroup)) CALL FlagError("Work group is not associated.",err,error,*999)

    IF(.NOT.ASSOCIATED(workGroup%computationEnvironment)) THEN
      localError="Computation environment is not associated for work group number "// &
        & TRIM(NumberToVString(workGroup%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    context=>workGroup%computationEnvironment%context
    IF(.NOT.ASSOCIATED(context)) THEN
      localError="The context is not associated for the computation environment for work group number "// &
        & TRIM(NumberToVString(workGroup%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("WorkGroup_ContextGet")
    RETURN
999 NULLIFY(context)
998 ERRORSEXITS("WorkGroup_ContextGet",err,error)
    RETURN 1
    
  END SUBROUTINE WorkGroup_ContextGet

  !
  !================================================================================================================================
  !

  !>Finds and returns a pointer to the work group with the given user number.
  SUBROUTINE WorkGroup_Get(computationEnvironment,workGroupUserNumber,workGroup,err,error,*)

    !Argument variables
    TYPE(ComputationEnvironmentType), POINTER, INTENT(IN) :: computationEnvironment !<The computational environment
    INTEGER(INTG), INTENT(IN) :: workGroupUserNumber !<The user number of the work group to get
    TYPE(WorkGroupType), POINTER, INTENT(OUT) :: workGroup !<On return, the work group corresponding to the user number. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("WorkGroup_Get",err,error,*999)

    CALL WorkGroup_UserNumberFind(computationEnvironment,workGroupUserNumber,workGroup,err,error,*999)
    IF(.NOT.ASSOCIATED(workGroup)) THEN
      localError="A work group with an user number of "//TRIM(NumberToVString(workGroupUserNumber,"*",err,error))// &
        & " does not exist."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("WorkGroup_Get")
    RETURN
999 ERRORSEXITS("WorkGroup_Get",err,error)
    RETURN 1
    
  END SUBROUTINE WorkGroup_Get

  !
  !================================================================================================================================
  !

  !>Gets the group communicator from a work group.
  SUBROUTINE WorkGroup_GroupCommunicatorGet(workGroup,groupCommunicator,err,error,*)

    !Argument variables
    TYPE(WorkGroupType), POINTER, INTENT(IN) :: workGroup !<The work group to get the group communicator for.
    INTEGER(INTG), INTENT(OUT) :: groupCommunicator !<On exit, the group communicator for the work group.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("WorkGroup_GroupCommunicatorGet",err,error,*999)

    !Check input arguments
    IF(.NOT.ASSOCIATED(workGroup)) CALL FlagError("Work group is not associated.",err,error,*999)

    groupCommunicator=workGroup%mpiGroupCommunicator
    
    EXITS("WorkGroup_GroupCommunicatorGet")
    RETURN
999 ERRORSEXITS("WorkGroup_GroupCommunicatorGet",err,error)
    RETURN 1
    
  END SUBROUTINE WorkGroup_GroupCommunicatorGet

  !
  !================================================================================================================================
  !

  !>Gets the number/rank of the computation node in the work group communicator
  SUBROUTINE WorkGroup_GroupNodeNumberGet(workGroup,groupNodeNumber,err,error,*)
      
    !Argument Variables
    TYPE(WorkGroupType), POINTER, INTENT(IN) :: workGroup !<The work group to get the group node number for.
    INTEGER(INTG), INTENT(OUT) :: groupNodeNumber !<On return, the node number in the group communicator.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("WorkGroup_GroupNodeNumberGet",err,error,*999)

    IF(.NOT.ASSOCIATED(workGroup)) CALL FlagError("Work group is not associated.",err,error,*999)
    
    groupNodeNumber=workGroup%myGroupComputationNodeNumber
        
    EXITS("WorkGroup_GroupNodeNumberGet")
    RETURN
999 ERRORSEXITS("WorkGroup_GroupNodeNumberGet",err,error)
    RETURN
    
  END SUBROUTINE WorkGroup_GroupNodeNumberGet

  !
  !=================================================================================================================================
  !

  !>Returns the character label of a work group \see OpenCMISS::Iron::cmfe_WorkGroup_LabelGet
  SUBROUTINE WorkGroup_LabelGetC(workGroup,label,err,error,*)

    !Argument Variables
    TYPE(WorkGroupType), POINTER, INTENT(IN) :: workGroup !<The work group to get the label for
    CHARACTER(LEN=*), INTENT(OUT) :: label !<On exit, the label for the work group
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("WorkGroup_LabelGetC",err,error,*999)

    IF(.NOT.ASSOCIATED(workGroup)) CALL FlagError("Work group is not associated.",err,error,*999)

    label=CHAR(workGroup%label)
    
    EXITS("WorkGroup_LabelGetC")
    RETURN
999 ERRORSEXITS("WorkGroup_LabelGetC",err,error)
    RETURN 1
    
  END SUBROUTINE WorkGroup_LabelGetC

  !
  !=================================================================================================================================
  !

  !>Returns the varying string label of a work group \see OpenCMISS::Iron::cmfe_WorkGroup_LabelGet
  SUBROUTINE WorkGroup_LabelGetVS(workGroup,label,err,error,*)

    !Argument Variables
    TYPE(WorkGroupType), POINTER, INTENT(IN) :: workGroup !<The work group to get the label for
    TYPE(VARYING_STRING), INTENT(OUT) :: label !<On exit, the label for the work group
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("WorkGroup_LabelGetVS",err,error,*999)

    IF(.NOT.ASSOCIATED(workGroup)) CALL FlagError("Work group is not associated.",err,error,*999)

    label=workGroup%label
    
    EXITS("WorkGroup_LabelGetVS")
    RETURN
999 ERRORSEXITS("WorkGroup_LabelGetVS",err,error)
    RETURN 1
    
  END SUBROUTINE WorkGroup_LabelGetVS

  !
  !================================================================================================================================
  !
  
  !>Gets the number of computation nodes in the work group communicator.
  SUBROUTINE WorkGroup_NumberOfGroupNodesGet(workGroup,numberOfGroupNodes,err,error,*)
     
    !Argument Variables
    TYPE(WorkGroupType), POINTER, INTENT(IN) :: workGroup !<The work group to get the group number of nodes for.
    INTEGER(INTG), INTENT(OUT) :: numberOfGroupNodes !<On return, the number of nodes in the group communicator.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("WorkGroup_NumberOfGroupNodesGet",err,error,*999)

    numberOfGroupNodes=workGroup%numberOfGroupComputationNodes
    
    EXITS("WorkGroup_NumberOfGroupNodesGet")
    RETURN
999 ERRORSEXITS("WorkGroup_NumberOfGroupNodesGet",err,error)
    RETURN
    
  END SUBROUTINE WorkGroup_NumberOfGroupNodesGet

  !
  !================================================================================================================================
  !

  !>Gets the parent work group from a work group.
  SUBROUTINE WorkGroup_ParentWorkGroupGet(workGroup,parentWorkGroup,err,error,*)

    !Argument variables
    TYPE(WorkGroupType), POINTER, INTENT(IN) :: workGroup !<The work group to get the parent work group for.
    TYPE(WorkGroupType), POINTER, INTENT(OUT) :: parentWorkGroup !<On exit, a pointer to the parent work group for the work group. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("WorkGroup_ParentWorkGroupGet",err,error,*999)

    !Check input arguments
    IF(ASSOCIATED(parentWorkGroup)) CALL FlagError("Parent work group is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(workGroup)) CALL FlagError("Work group is not associated.",err,error,*999)

    !Get the parent work group
    parentWorkGroup=>workGroup%parentWorkGroup
    !Check parent work group is associated.
    IF(.NOT.ASSOCIATED(parentWorkGroup)) THEN
      localError="Parent work group is not associated for work group "// &
        & TRIM(NumberToVString(workGroup%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("WorkGroup_ParentWorkGroupGet")
    RETURN
999 NULLIFY(parentWorkGroup)
998 ERRORSEXITS("WorkGroup_ParentWorkGroupGet",err,error)
    RETURN 1
    
  END SUBROUTINE WorkGroup_ParentWorkGroupGet

  !
  !================================================================================================================================
  !

  !>Finds and returns a pointer to the work group with the given user number. If no work group with that number exists work group is left nullified.
  SUBROUTINE WorkGroup_UserNumberFind(computationEnvironment,userNumber,workGroup,err,error,*)

    !Argument variables
    TYPE(ComputationEnvironmentType), POINTER :: computationEnvironment !<The computation environment containing the work group.
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the work Group to find
    TYPE(WorkGroupType), POINTER, INTENT(OUT) :: workGroup !<On exit, a pointer to the work group with the specified user number if it exists. If no work group exists with the specified user number a NULL pointer is returned. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(WorkGroupType), POINTER :: worldWorkGroup
    
    ENTERS("WorkGroup_UserNumberFind",err,error,*999)

    IF(ASSOCIATED(workGroup)) CALL FlagError("Work group is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(computationEnvironment)) CALL FlagError("Computation environment is not associated.",err,error,*999)
    
    worldWorkGroup=>computationEnvironment%worldWorkGroup
    IF(.NOT.ASSOCIATED(worldWorkGroup)) CALL FlagError("World work group is not associated.",err,error,*999)
    
    NULLIFY(workGroup)
    IF(userNumber==0) THEN
      workGroup=>worldWorkGroup
    ELSE
      CALL WorkGroup_UserNumberFindPtr(userNumber,worldWorkGroup,workGroup,err,error,*999)
    ENDIF
  
    EXITS("WorkGroup_UserNumberFind")
    RETURN
999 ERRORSEXITS("WorkGroup_UserNumberFind",err,error)
    RETURN 1
    
  END SUBROUTINE WorkGroup_UserNumberFind

  !
  !================================================================================================================================
  !

  !>Finds and returns a pointer to the work group with the given user number starting from the given start work group and searching all sub-groups under the start work group. If no work group with that number exists work group is left nullified.
  RECURSIVE SUBROUTINE WorkGroup_UserNumberFindPtr(userNumber,startWorkGroup,workGroup,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number to find
    TYPE(WorkGroupType), POINTER, INTENT(IN) :: startWorkGroup !<A pointer to the work group to start the search from
    TYPE(WorkGroupType), POINTER, INTENT(OUT) :: workGroup !<On exit, a pointer to the work group with the specified user number if it exists. If no work group exists with the specified user number a NULL pointer is returned.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: workGroupIdx

    ENTERS("WorkGroup_UserNumberFindPtr",err,error,*999)

    IF(.NOT.ASSOCIATED(startWorkGroup)) CALL FlagError("Start work group is not associated",err,error,*999)
    IF(ASSOCIATED(workGroup)) CALL FlagError("Work group is already associated.",err,error,*999)

    NULLIFY(workGroup)
    IF(startWorkGroup%userNumber==userNumber) THEN
      workGroup=>startWorkGroup
    ELSE
      IF(ALLOCATED(startWorkGroup%subGroups)) THEN
        DO workGroupIdx=1,startWorkGroup%numberOfSubGroups
          CALL WorkGroup_UserNumberFindPtr(userNumber,startWorkGroup%subGroups(workGroupIdx)%ptr,workGroup,err,error,*999)
          IF(ASSOCIATED(workGroup)) EXIT
        ENDDO !workGroupIdx
      ENDIF
    ENDIF
    
    EXITS("WorkGroup_UserNumberFindPtr")
    RETURN
999 ERRORSEXITS("WorkGroup_UserNumberFindPtr",err,error)
    RETURN 1
    
  END SUBROUTINE WorkGroup_UserNumberFindPtr

  !
  !================================================================================================================================
  !

  !>Returns the user number for a work group.
  SUBROUTINE WorkGroup_UserNumberGet(workGroup,userNumber,err,error,*)

    !Argument variables
    TYPE(WorkGroupType), POINTER, INTENT(IN) :: workGroup !<A pointer to the work group to get the user number for
    INTEGER(INTG), INTENT(OUT) :: userNumber !<On exit, the user number of the work group.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("WorkGroup_UserNumberGet",err,error,*999)

    IF(.NOT.ASSOCIATED(workGroup)) CALL FlagError("Work group is not associated.",err,error,*999)

    userNumber=workGroup%userNumber
  
    EXITS("WorkGroup_UserNumberGet")
    RETURN
999 ERRORSEXITS("WorkGroup_UserNumberGet",err,error)
    RETURN 1
    
  END SUBROUTINE WorkGroup_UserNumberGet

  !
  !================================================================================================================================
  !

  !>Gets a sub work group from a work group.
  SUBROUTINE WorkGroup_WorkSubGroupGet(workGroup,subGroupIdx,subWorkGroup,err,error,*)

    !Argument variables
    TYPE(WorkGroupType), POINTER, INTENT(IN) :: workGroup !<The work group to get the sub group for.
    INTEGER(INTG), INTENT(IN) :: subGroupIdx !<The sub group index to get
    TYPE(WorkGroupType), POINTER, INTENT(OUT) :: subWorkGroup !<On exit, a pointer to the sub work group of the work group for the subGroupIdx'th sub group. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("WorkGroup_WorkSubGroupGet",err,error,*999)

    !Check input arguments
    IF(ASSOCIATED(subWorkGroup)) CALL FlagError("Sub work group is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(workGroup)) CALL FlagError("Work group is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(workGroup%subGroups)) CALL FlagError("Work group sub groups is not allocated.",err,error,*999)
    IF(subGroupIdx<1.OR.subGroupIdx>workGroup%numberOfSubGroups) THEN
      localError="The specified sub group index of "//TRIM(NumberToVString(subGroupIdx,"*",err,error))// &
        & " is invalid. The sub group index must be >=1 and <= "// &
        & TRIM(NumberToVString(workGroup%numberOfSubGroups,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    !Get the parent work group
    subWorkGroup=>workGroup%subGroups(subGroupIdx)%ptr
    !Check sub work group is associated.
    IF(.NOT.ASSOCIATED(subWorkGroup)) THEN
      localError="The sub work group is not associated for sub group index "// &
        & TRIM(NumberToVString(subGroupIdx,"*",err,error))//" of work group "// &
        & TRIM(NumberToVString(workGroup%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("WorkGroup_WorkSubGroupGet")
    RETURN
999 NULLIFY(subWorkGroup)
998 ERRORSEXITS("WorkGroup_WorkSubGroupGet",err,error)
    RETURN 1
    
  END SUBROUTINE WorkGroup_WorkSubGroupGet

  !
  !================================================================================================================================
  !

END MODULE ComputationAccessRoutines
