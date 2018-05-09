!> \file
!> \author Chris Bradley
!> \brief This module handles all context access routines
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

!> This module handles all contex access routines.
MODULE ContextAccessRoutines

  USE BaseRoutines
  USE ISO_VARYING_STRING
  USE Kinds
  USE Strings
  USE Types

#include "macros.h"

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC Context_BasisFunctionsGet

  PUBLIC Context_ComputationEnvironmentGet

  PUBLIC Context_ContextsGet

  PUBLIC Context_CoordinateSystemsGet
  
  PUBLIC Context_ProblemsGet

  PUBLIC Context_RandomSeedsGet

  PUBLIC Context_RandomSeedsSizeGet

  PUBLIC Context_RegionsGet
  
  PUBLIC Context_UserNumberFind

  PUBLIC Context_UserNumberGet

  PUBLIC Context_Get

CONTAINS

  !
  !================================================================================================================================
  !
  
  !>Gets the basis functions for a context.
  SUBROUTINE Context_BasisFunctionsGet(context,basisFunctions,err,error,*)

    !Argument Variables
    TYPE(ContextType), POINTER, INTENT(IN) :: context !<The context to get the basis functions for.
    TYPE(BasisFunctionsType), POINTER :: basisFunctions !<On return, the basis functions for the context. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Context_BasisFunctionsGet",err,error,*998)

    IF(ASSOCIATED(basisFunctions)) CALL FlagError("Basis functions is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(context)) CALL FlagError("Context is not associated.",err,error,*999)
    
    basisFunctions=>context%basisFunctions
    IF(.NOT.ASSOCIATED(basisFunctions)) THEN
      localError="The basis functions for context number "// &
        & TRIM(NumberToVString(context%userNumber,"*",err,error))//" is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
 
    EXITS("Context_BasisFunctionsGet")
    RETURN
999 NULLIFY(basisFunctions)
998 ERRORS("Context_BasisFunctionsGet",err,error)
    EXITS("Context_BasisFunctionsGet")
    RETURN 1
    
  END SUBROUTINE Context_BasisFunctionsGet

  !
  !================================================================================================================================
  !
  
  !>Gets the computation environment for a context.
  SUBROUTINE Context_ComputationEnvironmentGet(context,computationEnvironment,err,error,*)

    !Argument Variables
    TYPE(ContextType), POINTER, INTENT(IN) :: context !<The context to get the computation environment for.
    TYPE(ComputationEnvironmentType), POINTER :: computationEnvironment !<On return, the computation environment for the context. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Context_ComptationEnvironmentGet",err,error,*998)

    IF(ASSOCIATED(computationEnvironment)) CALL FlagError("Computation environment is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(context)) CALL FlagError("Context is not associated.",err,error,*999)
    
    computationEnvironment=>context%computationEnvironment
    IF(.NOT.ASSOCIATED(computationEnvironment)) THEN
      localError="The computation environment for context number "// &
        & TRIM(NumberToVString(context%userNumber,"*",err,error))//" is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
 
    EXITS("Context_ComputationEnvironmentGet")
    RETURN
999 NULLIFY(computationEnvironment)
998 ERRORS("Context_ComputationEnvironmentGet",err,error)
    EXITS("Context_ComputationEnvironmentGet")
    RETURN 1
    
  END SUBROUTINE Context_ComputationEnvironmentGet

  !
  !================================================================================================================================
  !
  
  !>Gets the contexts for a context.
  SUBROUTINE Context_ContextsGet(context,contexts,err,error,*)

    !Argument Variables
    TYPE(ContextType), POINTER, INTENT(IN) :: context !<The context to get the contexts for.
    TYPE(ContextsType), POINTER :: contexts !<On return, the contexts for the context. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Context_ContextsGet",err,error,*998)

    IF(ASSOCIATED(contexts)) CALL FlagError("Contexts is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(context)) CALL FlagError("Context is not associated.",err,error,*999)
    
    contexts=>context%contexts
    IF(.NOT.ASSOCIATED(contexts)) THEN
      localError="The contexts for context number "// &
        & TRIM(NumberToVString(context%userNumber,"*",err,error))//" is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
 
    EXITS("Context_ContextsGet")
    RETURN
999 NULLIFY(contexts)
998 ERRORS("Context_ContextsGet",err,error)
    EXITS("Context_ContextsGet")
    RETURN 1
    
  END SUBROUTINE Context_ContextsGet

  !
  !================================================================================================================================
  !
  
  !>Gets the coordinate systems for a context.
  SUBROUTINE Context_CoordinateSystemsGet(context,coordinateSystems,err,error,*)

    !Argument Variables
    TYPE(ContextType), POINTER, INTENT(IN) :: context !<The context to get the coordinate systems for.
    TYPE(CoordinateSystemsType), POINTER :: coordinateSystems !<On return, the coordinate systems for the context. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Context_CoordinateSystemsGet",err,error,*998)

    IF(ASSOCIATED(coordinateSystems)) CALL FlagError("Coordinate systems is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(context)) CALL FlagError("Context is not associated.",err,error,*999)
    
    coordinateSystems=>context%coordinateSystems
    IF(.NOT.ASSOCIATED(coordinateSystems)) THEN
      localError="The coordinate systems for context number "// &
        & TRIM(NumberToVString(context%userNumber,"*",err,error))//" is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
 
    EXITS("Context_CoordinateSystemsGet")
    RETURN
999 NULLIFY(coordinateSystems)
998 ERRORS("Context_CoordinateSystemsGet",err,error)
    EXITS("Context_CoordinateSystemsGet")
    RETURN 1
    
  END SUBROUTINE Context_CoordinateSystemsGet

  !
  !================================================================================================================================
  !
  
  !>Gets the problems for a context.
  SUBROUTINE Context_ProblemsGet(context,problems,err,error,*)

    !Argument Variables
    TYPE(ContextType), POINTER, INTENT(IN) :: context !<The context to get the problems for.
    TYPE(ProblemsType), POINTER :: problems !<On return, the problems for the context. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Context_ProblemsGet",err,error,*998)

    IF(ASSOCIATED(problems)) CALL FlagError("Problems is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(context)) CALL FlagError("Context is not associated.",err,error,*999)
    
    problems=>context%problems
    IF(.NOT.ASSOCIATED(problems)) THEN
      localError="The problems for context number "// &
        & TRIM(NumberToVString(context%userNumber,"*",err,error))//" is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
 
    EXITS("Context_ProblemsGet")
    RETURN
999 NULLIFY(problems)
998 ERRORS("Context_ProblemsGet",err,error)
    EXITS("Context_ProblemsGet")
    RETURN 1
    
  END SUBROUTINE Context_ProblemsGet

  !
  !================================================================================================================================
  !
  
  !>Gets the random seeds for a context. \see OpenCMISS::Iron::cmfe_Context_RandomSeedsGet
  SUBROUTINE Context_RandomSeedsGet(context,randomSeeds,err,error,*)

    !Argument Variables
    TYPE(ContextType), POINTER, INTENT(IN) :: context !<The context to get the random seeds for.
    INTEGER(INTG), INTENT(OUT) :: randomSeeds(:) !<On return, the random seeds for the context. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Context_RandomSeedsGet",err,error,*999)

    IF(.NOT.ASSOCIATED(context)) CALL FlagError("Context is not associated.",err,error,*999)
    
    IF(SIZE(randomSeeds,1)<SIZE(context%cmissRandomSeeds,1)) THEN
      localError="The size of the supplied random seeds array of "// &
        & TRIM(NumberToVString(SIZE(randomSeeds,1),"*",err,error))//" is too small. The size must be >= "// &
        & TRIM(NumberToVString(SIZE(context%cmissRandomSeeds,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    randomSeeds(1:SIZE(context%cmissRandomSeeds,1))=context%cmissRandomSeeds(1:SIZE(context%cmissRandomSeeds,1))
  
    EXITS("Context_RandomSeedsGet")
    RETURN
999 ERRORS("Context_RandomSeedsGet",err,error)
    EXITS("Context_RandomSeedsGet")
    RETURN 1
    
  END SUBROUTINE Context_RandomSeedsGet

  !
  !================================================================================================================================
  !

  !>Returns the size of the random seeds array for a context \see OpenCMISS::Iron::cmfe_Context_RandomSeedsSizeGet
  SUBROUTINE Context_RandomSeedsSizeGet(context,randomSeedsSize,err,error,*)
  
    !Argument variables
    TYPE(ContextType), POINTER, INTENT(IN) :: context !<The context to get the size of the random seeds array for.
    INTEGER(INTG), INTENT(OUT) :: randomSeedsSize !<On return, the size of the random seeds array.
    INTEGER(INTG), INTENT(INOUT) :: err !<The error string
    TYPE(VARYING_STRING), INTENT(INOUT) :: error !<The error code
    !Local Variables

    ENTERS("Context_RandomSeedsSizeGet",err,error,*999)

    IF(.NOT.ASSOCIATED(context)) CALL FlagError("Context is not associated.",err,error,*999)
    
    randomSeedsSize=SIZE(context%cmissRandomSeeds,1)
    
    EXITS("Context_RandomSeedsSizeGet")
    RETURN
999 ERRORSEXITS("Context_RandomSeedsSizeGet",err,error)
    RETURN 1
    
  END SUBROUTINE Context_RandomSeedsSizeGet

  !
  !================================================================================================================================
  !
  
  !>Gets the regions for a context.
  SUBROUTINE Context_RegionsGet(context,regions,err,error,*)

    !Argument Variables
    TYPE(ContextType), POINTER, INTENT(IN) :: context !<The context to get the regions for.
    TYPE(RegionsType), POINTER :: regions !<On return, the regions for the context. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Context_RegionsGet",err,error,*998)

    IF(ASSOCIATED(regions)) CALL FlagError("Regions is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(context)) CALL FlagError("Context is not associated.",err,error,*999)
    
    regions=>context%regions
    IF(.NOT.ASSOCIATED(regions)) THEN
      localError="The regions for context number "// &
        & TRIM(NumberToVString(context%userNumber,"*",err,error))//" is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
 
    EXITS("Context_RegionsGet")
    RETURN
999 NULLIFY(regions)
998 ERRORS("Context_RegionsGet",err,error)
    EXITS("Context_RegionsGet")
    RETURN 1
    
  END SUBROUTINE Context_RegionsGet

  !
  !================================================================================================================================
  !

  !>Finds and returns a pointer to the context identified by a user number. If no context with that user number exists context is left nullified.
  SUBROUTINE Context_UserNumberFind(contexts,userNumber,context,err,error,*)

    !Argument variables
    TYPE(ContextsType), INTENT(IN) :: contexts !<The contexts to find the user number in. 
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number to find.
    TYPE(ContextType), POINTER :: context !<On return, a pointer to the context with the given user number. If no context with that user number exists then the pointer is returned as NULL. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: contextIdx
    TYPE(VARYING_STRING) :: localError

    ENTERS("Context_UserNumberFind",err,error,*999)

    IF(ASSOCIATED(context)) CALL FlagError("Context is already associated.",err,error,*999)
    IF(contexts%numberOfContexts>0) THEN
      IF(.NOT.ALLOCATED(contexts%contexts)) CALL FlagError("Contexts contexts is not allocated.",err,error,*999)   
      NULLIFY(context)
      DO contextIdx=1,contexts%numberOfContexts
        IF(.NOT.ASSOCIATED(contexts%contexts(contextIdx)%ptr)) THEN
          localError="The context pointer in contexts is not associated for context index "// &
            & TRIM(NumberToVString(contextIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        IF(contexts%contexts(contextIdx)%ptr%userNumber==userNumber) THEN
          context=>contexts%contexts(contextIdx)%ptr
          EXIT
        ENDIF
      ENDDO !contextIdx
    ENDIF
    
    EXITS("Context_UserNumberFind")
    RETURN
999 ERRORSEXITS("Context_UserNumberFind",err,error)
    RETURN 1
    
  END SUBROUTINE Context_UserNumberFind

  !
  !================================================================================================================================
  !

  !>Returns the user number for a context.
  SUBROUTINE Context_UserNumberGet(context,userNumber,err,error,*)

    !Argument variables
    TYPE(ContextType), POINTER :: context !<A pointer to the context to get the user number for
    INTEGER(INTG), INTENT(OUT) :: userNumber !<On exit, the user number of the context.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Context_UserNumberGet",err,error,*999)

    IF(.NOT.ASSOCIATED(context)) CALL FlagError("Context is not associated.",err,error,*999)

    userNumber=context%userNumber
  
    EXITS("Context_UserNumberGet")
    RETURN
999 ERRORSEXITS("Context_UserNumberGet",err,error)
    RETURN 1
    
  END SUBROUTINE Context_UserNumberGet

  !
  !================================================================================================================================
  !

  !>Finds and returns a pointer to the context with the given user number. 
  SUBROUTINE Context_Get(contexts,userNumber,context,err,error,*)

    !Argument variables
    TYPE(ContextsType), INTENT(IN) :: contexts !<The contexts to get the user number for.
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the context to find
    TYPE(ContextType), POINTER :: context !<On exit, a pointer to the context with the specified user number if it exists. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Context_Get",err,error,*999)

    CALL Context_UserNumberFind(contexts,userNumber,context,err,error,*999)
    IF(.NOT.ASSOCIATED(context)) THEN
      localError="A context with an user number of "//TRIM(NumberToVString(userNumber,"*",err,error))//" does not exist."
      CALL FlagError(localError,err,error,*999)
    ENDIF
  
    EXITS("Context_Get")
    RETURN
999 ERRORSEXITS("Context_Get",err,error)
    RETURN 1
    
  END SUBROUTINE Context_Get

  !
  !================================================================================================================================
  !    

END MODULE ContextAccessRoutines
