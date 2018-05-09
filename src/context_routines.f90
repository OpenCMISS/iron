!> \file
!> \author Chris Bradley
!> \brief This module handles all context routines
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

!> This module handles all contex routines.
MODULE ContextRoutines

  USE BaseRoutines
  USE BasisRoutines
  USE ComputationRoutines
  USE ContextAccessRoutines
  USE COORDINATE_ROUTINES
  USE ISO_VARYING_STRING
  USE Kinds
  USE PROBLEM_ROUTINES
  USE REGION_ROUTINES
  USE Strings
  USE Types

#include "macros.h"

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC Context_Create

  PUBLIC Context_Destroy

  PUBLIC Context_RandomSeedsSet

  PUBLIC Contexts_Finalise,Contexts_Initialise
 
CONTAINS

  !
  !=================================================================================================================================
  !

  !>Create a context
  SUBROUTINE Context_Create(contexts,context,err,error,*)

    !Argument Variables
    TYPE(ContextsType), TARGET :: contexts !<the contexts to create the context for
    TYPE(ContextType), POINTER :: context !<On return, a pointer to the created context. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: contextIdx,dummyErr,i,randomSeedsSize,time(8)
    TYPE(ContextType), POINTER :: newContext
    TYPE(ContextPtrType), ALLOCATABLE :: newContexts(:)
    TYPE(VARYING_STRING) :: dummyError

    NULLIFY(newContext)

    ENTERS("Context_Create",err,error,*999)

    IF(ASSOCIATED(context)) CALL FlagError("Context is already associated.",err,error,*999)

    CALL Context_Initialise(newContext,err,error,*999)
    !Setup the random seeds based on the time
    CALL RANDOM_SEED(SIZE=randomSeedsSize)
    ALLOCATE(newContext%cmissRandomSeeds(randomSeedsSize),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate random seeds.",err,error,*999)
    newContext%cmissRandomSeeds(1:randomSeedsSize)=[(i,i=1,randomSeedsSize)]
    CALL DATE_AND_TIME(VALUES=time)
    newContext%cmissRandomSeeds(1)=3600000*time(5)+60000*time(6)+1000*time(7)+time(8)
    CALL RANDOM_SEED(PUT=newContext%cmissRandomSeeds)
    !Intialise the computation
    CALL Computation_Initialise(newContext,err,error,*999)
    !Intialise the basis functions
    CALL BasisFunctions_Initialise(newContext,err,error,*999)
    !Initialise the coordinate systems
    CALL CoordinateSystems_Initialise(newContext,err,error,*999)
    !Initialise the regions 
    CALL Regions_Initialise(newContext,err,error,*999)
    !Initialise the problems
    CALL Problems_Initialise(newContext,err,error,*999)
    !Add the new context to the list of contexts    
    contexts%lastContextUserNumber=contexts%lastContextUserNumber+1
    newContext%userNumber=contexts%lastContextUserNumber    
    newContext%contexts=>contexts    
    ALLOCATE(newContexts(contexts%numberOfContexts+1),STAT=err)
    DO contextIdx=1,contexts%numberOfContexts
      newContexts(contextIdx)%ptr=>contexts%contexts(contextIdx)%ptr
    ENDDO !contextIdx
    newContexts(contexts%numberOfContexts+1)%ptr=>newContext
    CALL MOVE_ALLOC(newContexts,contexts%contexts)
    contexts%numberOfContexts=contexts%numberOfContexts+1
    context=>newContext
    
    EXITS("Context_Create")
    RETURN
999 CALL Context_Finalise(newContext,dummyErr,dummyError,*998)
998 IF(ALLOCATED(newContexts)) DEALLOCATE(newContexts)
    ERRORSEXITS("Context_Create",err,error)    
    RETURN 1
    
  END SUBROUTINE Context_Create
  
  !
  !=================================================================================================================================
  !

  !>Destroy a context
  SUBROUTINE Context_Destroy(context,err,error,*)

    !Argument Variables
    TYPE(ContextType), POINTER :: context !<A pointer to the context to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: contextIdx,position
    TYPE(ContextPtrType), ALLOCATABLE :: newContexts(:)
    TYPE(ContextsType), POINTER :: contexts
    TYPE(VARYING_STRING) :: localError

    ENTERS("Context_Destroy",err,error,*999)

    IF(.NOT.ASSOCIATED(context)) CALL FlagError("Context is not associated.",err,error,*999)
    
    NULLIFY(contexts)
    CALL Context_ContextsGet(context,contexts,err,error,*999)
    position=0
    DO contextIdx=1,contexts%numberOfContexts
      IF(context%userNumber==contexts%contexts(contextIdx)%ptr%userNumber) THEN
        position=contextIdx
        EXIT
      ENDIF
    ENDDO !contextIdx
    IF(position==0) THEN
      localError="The specified context with a user number of "// &
        & TRIM(NumberToVString(context%userNumber,"*",err,error))// &
        & " could not be found in the list of contexts."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    ALLOCATE(newContexts(contexts%numberOfContexts-1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new contexts.",err,error,*999)
    DO contextIdx=1,contexts%numberOfContexts
      IF(contextIdx<position) THEN
        newContexts(contextIdx)%ptr=>contexts%contexts(contextIdx)%ptr
      ELSE IF(contextIdx>position) THEN
        newContexts(contextIdx-1)%ptr=>contexts%contexts(contextIdx)%ptr
      ENDIF
    ENDDO !contextIdx
    CALL MOVE_ALLOC(newContexts,contexts%contexts)
    contexts%numberOfContexts=contexts%numberOfContexts-1
    CALL Context_Finalise(context,err,error,*999)
    
    EXITS("Context_Destroy")
    RETURN
999 IF(ALLOCATED(newContexts)) DEALLOCATE(newContexts)
    ERRORSEXITS("Context_Destroy",err,error)    
    RETURN 1
    
  END SUBROUTINE Context_Destroy
  
  !
  !=================================================================================================================================
  !

  !>Finalise a context and deallocate all memory
  SUBROUTINE Context_Finalise(context,err,error,*)

    !Argument Variables
    TYPE(ContextType),POINTER :: context !<A pointer to the context to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Context_Finalise",err,error,*999)

    IF(ASSOCIATED(context)) THEN
      IF(ALLOCATED(context%cmissRandomSeeds)) DEALLOCATE(context%cmissRandomSeeds)
      !Finalise the problems
      CALL Problems_Finalise(context%problems,err,error,*999)
      !Finalise the regions
      CALL Regions_Finalise(context%regions,err,error,*999)
      !Finalise the coordinate systems
      CALL CoordinateSystems_Finalise(context%coordinateSystems,err,error,*999)
      !Finalise basis functions
      CALL BasisFunctions_Finalise(context%basisFunctions,err,error,*999)
      !Finalise computation
      CALL Computation_Finalise(context%computationEnvironment,err,error,*999)
      !Deallocate context
      DEALLOCATE(context)
    ENDIF
    
    EXITS("Context_Finalise")
    RETURN
999 ERRORSEXITS("Context_Finalise",err,error)    
    RETURN 1
    
  END SUBROUTINE Context_Finalise
  
  !
  !=================================================================================================================================
  !

  !>Initialise a context
  SUBROUTINE Context_Initialise(context,err,error,*)

    !Argument Variables
    TYPE(ContextType), POINTER :: context !<A pointer to the context to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("Context_Initialise",err,error,*998)

    IF(ASSOCIATED(context)) CALL FlagError("Context is already associated.",err,error,*998)
      
    ALLOCATE(context,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate context.",err,error,*999)
    context%userNumber=0
    NULLIFY(context%contexts)
    NULLIFY(context%coordinateSystems)
    NULLIFY(context%basisFunctions)
    NULLIFY(context%computationEnvironment)
    NULLIFY(context%problems)
    NULLIFY(context%regions)
   
    EXITS("Context_Initialise")
    RETURN
999 CALL Context_Finalise(context,dummyErr,dummyError,*998)
998 ERRORSEXITS("Context_Initialise",err,error)    
    RETURN 1
    
  END SUBROUTINE Context_Initialise
  
  !
  !================================================================================================================================
  !

  !>Sets the random seeds for a context \see OpenCMISS::Iron::cmfe_Context_RandomSeedsSet
  SUBROUTINE Context_RandomSeedsSet(context,randomSeeds,err,error,*)
  
    !Argument variables 
    TYPE(ContextType), POINTER :: context !<A pointer to the context to set the random seeds for
    INTEGER(INTG), INTENT(IN) :: randomSeeds(:) !<The random seeds to set. 
    INTEGER(INTG), INTENT(INOUT) :: err !<The error string
    TYPE(VARYING_STRING), INTENT(INOUT) :: error !<The error code
    !Local Variables
    
    ENTERS("Context_RandomSeedsSet",err,error,*999)

    IF(.NOT.ASSOCIATED(context)) CALL FlagError("Context is not associated.",err,error,*999)
    
    IF(SIZE(randomSeeds,1)>SIZE(context%cmissRandomSeeds,1)) THEN
      context%cmissRandomSeeds(1:SIZE(context%cmissRandomSeeds,1))=randomSeeds(1:SIZE(context%cmissRandomSeeds,1))
    ELSE
      context%cmissRandomSeeds(1:SIZE(randomSeeds,1))=randomSeeds(1:SIZE(randomSeeds,1))
    ENDIF

    EXITS("Context_RandomSeedsSet")
    RETURN
999 ERRORSEXITS("Context_RandomSeedsSet",err,error)
    RETURN 1
    
  END SUBROUTINE Context_RandomSeedsSet

  !
  !=================================================================================================================================
  !

  !>Finalise the contexts
  SUBROUTINE Contexts_Finalise(contexts,err,error,*)

    !Argument Variables
    TYPE(ContextsType) :: contexts !<The contexts to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Contexts_Finalise",err,error,*999)

    contexts%numberOfContexts=0
    
    EXITS("Contexts_Finalise")
    RETURN
999 ERRORSEXITS("Contexts_Finalise",err,error)    
    RETURN 1
    
  END SUBROUTINE Contexts_Finalise
  
  !
  !=================================================================================================================================
  !

  !>Initialise the contexts
  SUBROUTINE Contexts_Initialise(contexts,err,error,*)

    !Argument Variables
    TYPE(ContextsType) :: contexts !<The contexts to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Contexts_Initialise",err,error,*999)

    contexts%numberOfContexts=0
    contexts%lastContextUserNumber=0
    
    EXITS("Contexts_Initialise")
    RETURN
999 ERRORSEXITS("Contexts_Initialise",err,error)    
    RETURN 1
    
  END SUBROUTINE Contexts_Initialise
  
  !
  !================================================================================================================================
  !

END MODULE ContextRoutines
