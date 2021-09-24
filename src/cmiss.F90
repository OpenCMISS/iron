!> \file
!> \author Chris Bradley
!> \brief The top level OpenCMISS Iron module.
!>
!> \mainpage OpenCMISS Documentation
!>
!> An open source interactive computer program for Continuum Mechanics, Image analysis, Signal processing and System
!> Identification. Target usage: Bioengineering application of finite element analysis, boundary element and collocation
!> techniques.
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
!>

!> \defgroup OpenCMISS_CMISS OpenCMISS::Iron::CMISS
!> The top level cmiss module.
MODULE Cmiss

  USE ISO_C_BINDING
  
  USE BaseRoutines
  USE Constants
  USE ContextRoutines
  USE ComputationAccessRoutines
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

  !> \addtogroup CMISS_Constants OpenCMISS::Iron::CMISS::Constants
  !> \brief OpenCMISS CMISS constants
  !>@{
  
  !> \addtogroup CMISS_Versions OpenCMISS::Iron::CMISS::Constants::Versions
  !> \brief CMISS version parameters
  !>@{
  INTEGER(INTG), PARAMETER :: CMFE_MAJOR_VERSION = 0
  INTEGER(INTG), PARAMETER :: CMFE_MINOR_VERSION = 4
  INTEGER(INTG), PARAMETER :: CMFE_REVISION_VERSION = 0

  CHARACTER(LEN=MAXSTRLEN), PARAMETER :: CMFE_BUILD_VERSION = "$Rev"

  !> \addtogroup CMFE_ErrorHandlingModes OpenCMISS::Iron::CMISS::Constants::ErrorHandlingModes
  !> \brief Error handling mode parameters
  !> \see OpenCMISS
  !>@{
  INTEGER(INTG), PARAMETER :: CMFE_RETURN_ERROR_CODE = 0 !<Just return the error code \see cmfe_ErrorHandlingModes,OpenCMISS
  INTEGER(INTG), PARAMETER :: CMFE_OUTPUT_ERROR = 1 !<Output the error traceback and return the error code \see cmfe_ErrorHandlingModes,OpenCMISS
  INTEGER(INTG), PARAMETER :: CMFE_TRAP_ERROR = 2 !<Trap the error by outputing the error traceback and stopping the program \see cmfe_ErrorHandlingModes,OpenCMISS
  !>@}
  !>@}
  
  !Module types

  !Module variables

  LOGICAL, SAVE :: cmfeFirstInit = .FALSE. !<cmfeFirstInit will be .TRUE. if cmfe has ever been initialised, .FALSE. if not.

  INTEGER(INTG), SAVE :: cmfe_ErrorHandlingMode !<The current error handling mode for OpenCMISS \see cmfe_ErrorHandlingModes

  TYPE(ContextsType), SAVE, TARGET :: contexts !<The contexts for OpenCMISS
 
  !Interfaces

  INTERFACE

    SUBROUTINE cmfe_InitFatalHandler() BIND(C,NAME="cmfe_InitFatalHandler")
    END SUBROUTINE cmfe_InitFatalHandler

    SUBROUTINE cmfe_ResetFatalHandler() BIND(C,NAME="cmfe_ResetFatalHandler")
    END SUBROUTINE cmfe_ResetFatalHandler
    
    SUBROUTINE cmfe_SetFatalHandler() BIND(C,NAME="cmfe_SetFatalHandler")
    END SUBROUTINE cmfe_SetFatalHandler

  END INTERFACE

  PUBLIC CMFE_MAJOR_VERSION,CMFE_MINOR_VERSION,CMFE_REVISION_VERSION,CMFE_BUILD_VERSION

  PUBLIC CMFE_RETURN_ERROR_CODE,CMFE_OUTPUT_ERROR,CMFE_TRAP_ERROR

  PUBLIC contexts

  PUBLIC cmfe_ErrorHandlingModeGet_,cmfe_ErrorHandlingModeSet_
  
  PUBLIC cmfe_HandleError
  
  PUBLIC cmfe_Finalise_,cmfe_Initialise_

CONTAINS

  !
  !================================================================================================================================
  !

!!TODO Underscore to avoid name clash. Can be removed upon prefix rename.

  !>Returns the error handling mode for OpenCMISS \see OpenCMISS::Iron::cmfe_ErrorHandlingModeGet
  SUBROUTINE cmfe_ErrorHandlingModeGet_(errorHandlingMode,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: errorHandlingMode !<On return, the error handling mode. \see cmfe_ErrorHandlingModes,OpenCMISS
    INTEGER(INTG), INTENT(INOUT) :: err !<The error string
    TYPE(VARYING_STRING), INTENT(INOUT) :: error !<The error code
    !Local Variables

    ENTERS("cmfe_ErrorHandlingModeGet_",err,error,*999)

    errorHandlingMode=cmfe_ErrorHandlingMode
    
    EXITS("cmfe_ErrorHandlingModeGet_")
    RETURN
999 ERRORSEXITS("cmfe_ErrorHandlingModeGet_",err,error)
    RETURN 1
    
  END SUBROUTINE cmfe_ErrorHandlingModeGet_

  !
  !================================================================================================================================
  !

!!TODO Underscore to avoid name clash. Can be removed upon prefix rename.

  !>Sets the error handling mode for cmiss \see OpenCMISS::Iron::cmfe_ErrorHandlingModeSet
  SUBROUTINE cmfe_ErrorHandlingModeSet_(errorHandlingMode,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: errorHandlingMode !<The error handling mode to set. \see cmfe_ErrorHandlingModes,OpenCMISS
    INTEGER(INTG), INTENT(INOUT) :: err !<The error string
    TYPE(VARYING_STRING), INTENT(INOUT) :: error !<The error code
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("cmfe_ErrorHandlingModeSet_",err,error,*999)

    SELECT CASE(errorHandlingMode)
    CASE(CMFE_RETURN_ERROR_CODE)
      cmfe_ErrorHandlingMode=CMFE_RETURN_ERROR_CODE
    CASE(CMFE_OUTPUT_ERROR)
      cmfe_ErrorHandlingMode=CMFE_OUTPUT_ERROR
    CASE(CMFE_TRAP_ERROR)
      cmfe_ErrorHandlingMode=CMFE_TRAP_ERROR
    CASE DEFAULT
      localError="The supplied error handling mode of "//TRIM(NumberToVString(errorHandlingMode,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("cmfe_ErrorHandlingModeSet_")
    RETURN
999 ERRORSEXITS("cmfe_ErrorHandlingModeSet_",err,error)
    RETURN 1
    
  END SUBROUTINE cmfe_ErrorHandlingModeSet_

  !
  !================================================================================================================================
  !

!!TODO Underscore to avoid name clash. Can be removed upon prefix rename.
  
  !>Finalises OpenCMISS. \see OpenCMISS::Iron::cmfe_Finalise
  SUBROUTINE cmfe_Finalise_(context,err,error,*)
  
    !Argument variables
    TYPE(ContextType), POINTER :: context !<The context to finalise
    INTEGER(INTG), INTENT(INOUT) :: err !<The error string
    TYPE(VARYING_STRING), INTENT(INOUT) :: error !<The error code
    !Local Variables

    !Destroy context
    CALL Context_Destroy(context,err,error,*999)

    IF(contexts%numberOfContexts==0) THEN
      !Reset the signal handler
      CALL cmfe_ResetFatalHandler()
      !Finalise contexts
      CALL Contexts_Finalise(contexts,err,error,*999)
      !Finalise the base routines
      CALL BaseRoutines_Finalise(err,error,*999)
      !Reset first init
      cmfeFirstInit=.FALSE.
    ENDIF
     
    RETURN
999 RETURN 1
    
  END SUBROUTINE cmfe_Finalise_

  !
  !================================================================================================================================
  !

!!TODO Underscore to avoid name clash. Can be removed upon prefix rename.

  !>Initialises OpenCMISS. \see OpenCMISS::Iron::cmfe_Initialise
  SUBROUTINE cmfe_Initialise_(newContext,err,error,*)
  
    !Argument variables
    TYPE(ContextType), POINTER :: newContext !<On return, a pointer to the new context. Must not be associated on entry.
    INTEGER(INTG), INTENT(INOUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(INOUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: myWorldComputationNodeNumber
    TYPE(VARYING_STRING) :: versionString

    IF(.NOT.cmfeFirstInit) THEN
      !Initialise error mode
      cmfe_ErrorHandlingMode = CMFE_OUTPUT_ERROR !Default for now, maybe make CMFE_RETURN_ERROR_CODE the default
      !Initialise the base routines
      CALL BaseRoutines_Initialise(err,error,*999)
      !Initialise contexts
      CALL Contexts_Initialise(contexts,err,error,*999)
      !Setup signal handling
      CALL cmfe_InitFatalHandler()
      CALL cmfe_SetFatalHandler()
    ENDIF
    
    !Create new context
    CALL Context_Create(contexts,newContext,err,error,*999)

    IF(.NOT.cmfeFirstInit) THEN
      !Write out the CMISS version
      CALL ComputationEnvironment_WorldNodeNumberGet(newContext%computationEnvironment,myWorldComputationNodeNumber, &
        & err,error,*999)
      IF(myWorldComputationNodeNumber==0) THEN
        versionString="OpenCMISS(Iron) version "//TRIM(NumberToVString(CMFE_MAJOR_VERSION,"*",err,error))
        versionString=versionString//"."
        versionString=versionString//TRIM(NumberToVString(CMFE_MINOR_VERSION,"*",err,error))
        versionString=versionString//"."
        versionString=versionString//TRIM(NumberToVString(CMFE_REVISION_VERSION,"*",err,error))
        !versionString=versionString//" ("
        !versionString=versionString//TRIM(CMFE_BUILD_VERSION(6:))
        !versionString=versionString//" )"        
        !WRITE(*,'(A)') CHAR(versionString)
        versionString=""
      ENDIF
      !Set first initalised
      cmfeFirstInit = .TRUE.
    ENDIF
    
    RETURN
999 RETURN 1
    
  END SUBROUTINE cmfe_Initialise_

  !
  !================================================================================================================================
  !

  !>Handle an error condition
  SUBROUTINE cmfe_HandleError(err,error)
  
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(INOUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: mpiError
    
    SELECT CASE(cmfe_ErrorHandlingMode)
    CASE(CMFE_RETURN_ERROR_CODE)
      !Do nothing
    CASE(CMFE_OUTPUT_ERROR)
      CALL WriteError(err,error,*999)
    CASE(CMFE_TRAP_ERROR)
      CALL WriteError(err,error,*999)
      CALL MPI_ABORT(MPI_COMM_WORLD,err,mpiError)
      STOP
    CASE DEFAULT
      !Do nothing
    END SELECT

    RETURN
999 RETURN

  END SUBROUTINE cmfe_HandleError
  
  !
  !================================================================================================================================
  !

END MODULE Cmiss
