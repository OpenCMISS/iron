!> \file
!> \author Chris Bradley
!> \brief This module contains OpenCMISS MPI routines.
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

!> This module contains OpenCMISS MPI routines.
MODULE CmissMPI
  
  USE BaseRoutines
  USE Constants
  USE Kinds
#ifndef NOMPIMOD
  USE MPI
#endif
  USE ISO_VARYING_STRING
  USE Strings

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

  PUBLIC MPI_ErrorCheck

CONTAINS

  !
  !================================================================================================================================
  !

  !>Checks to see if an MPI error has occured during an MPI call and flags a OpenCMISS error it if it has.
  SUBROUTINE MPI_ErrorCheck(routine,MPIErrCode,err,error,*)
  
    !Argument Variables
    CHARACTER(LEN=*) :: routine !<The name of the MPI routine that has just been called.
    INTEGER(INTG), INTENT(IN) :: MPIErrCode !<The MPI error code returned from the MPI routine.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: MPIIerror, MPIErrStrLength
    CHARACTER(LEN=MAXSTRLEN) :: MPIErrStr
    TYPE(VARYING_STRING) :: localError

    ENTERS("MPI_ErrorCheck",err,error,*999)

    IF(MPIErrCode/=MPI_SUCCESS) THEN
      CALL MPI_ERROR_STRING(MPIErrCode,MPIErrStr,MPIErrStrLength,MPIIerror)
      localError="MPI error "//TRIM(NumberToVString(MPIErrCode,"*",err,error))//" ("// &
        & MPIErrStr(1:MPIErrStrLength)//") in "//routine(1:LEN_TRIM(routine))
      CALL FlagError(localError,err,error,*999)
    ENDIF

    EXITS("MPI_ErrorCheck")
    RETURN
999 ERRORSEXITS("MPI_ErrorCheck",err,error)
    RETURN 1
    
  END SUBROUTINE MPI_ErrorCheck

  !
  !================================================================================================================================
  !
    
END MODULE CmissMPI
