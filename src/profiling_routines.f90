!> \file
!> \author Chris Bradley
!> \brief This module contains all profiling routines.
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

!> This module contains all profiling routines.
MODULE ProfilingRoutines
  
  USE BaseRoutines
  USE INPUT_OUTPUT
  USE Kinds
  USE Strings
  USE Types
  USE ISO_VARYING_STRING
  
#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables
  
  !Interfaces
  
  PUBLIC Profiling_TimingsOutput

CONTAINS

  !
  !================================================================================================================================
  !

  !>Output timing information .
  SUBROUTINE Profiling_TimingsOutput(index,string,userTime,systemTime,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: index !<The output line index. 0 will output a header, 1 will output a line format with total time, 2 will output a line format without totaltime
    CHARACTER(LEN=*), INTENT(IN) :: string !<The string to output
    REAL(SP), INTENT(IN) :: userTime !<The user time
    REAL(SP), INTENT(IN) :: systemTime !<The system time
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: stringLength
    CHARACTER(LEN=45) :: localString
    CHARACTER(LEN=MAXSTRLEN) :: outputString
    TYPE(VARYING_STRING) :: localError

    ENTERS("Profiling_TimingsOutput",err,error,*999) 

    SELECT CASE(index)
    CASE(0)
      WRITE(outputString,'(49X,"    User time   System time    Total time")')
      CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
    CASE(1)
      stringLength=LEN_TRIM(string)
      IF(stringLength<=45) THEN
        localString=ADJUSTL(string(1:LEN_TRIM(string)))
      ELSE
        localString=ADJUSTL(string(1:45))
      ENDIF
      WRITE(outputString,'(A,4X,E13.6,X,E13.6,X,E13.6)') localString,userTime,systemTime,userTime+systemTime
    CASE(2)
      stringLength=LEN_TRIM(string)
      IF(stringLength<=45) THEN
        localString=ADJUSTL(string(1:LEN_TRIM(string)))
      ELSE
        localString=ADJUSTL(string(1:45))
      ENDIF
      WRITE(outputString,'(A,4X,E13.6,X,E13.6)') localString,userTime,systemTime
    CASE DEFAULT
      localError="The specified index of "//TRIM(NumberToVString(index,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    CALL WriteString(GENERAL_OUTPUT_TYPE,outputString,err,error,*999)      
      
       
    EXITS("Profiling_TimingsOutput")
    RETURN
999 ERRORSEXITS("Profiling_TimingsOutput",err,error)
    RETURN 1
    
  END SUBROUTINE Profiling_TimingsOutput
  
  !
  !================================================================================================================================
  !

END MODULE ProfilingRoutines
