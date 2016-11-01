!> \file
!> \author Chris Bradley
!> \brief This module contains all data projection access method routines.
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

!> This module contains all data projection access method routines.
MODULE DataProjectionAccessRoutines
  
  USE BASE_ROUTINES
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

  PUBLIC DataProjection_ResultMaximumErrorGet

  PUBLIC DataProjection_ResultMinimumErrorGet

  PUBLIC DataProjection_ResultRMSErrorGet

CONTAINS

  !
  !================================================================================================================================
  !

  !>Gets the projection maximum error for a data projection.
  SUBROUTINE DataProjection_ResultMaximumErrorGet(dataProjection,maximumDataPoint,maximumError,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to get the maximum error for
    INTEGER(INTG), INTENT(OUT) :: maximumDataPoint !<On exit, the data point number of the maximum error.
    REAL(DP), INTENT(OUT) :: maximumError !<On exit, the maximum error for the data projection.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
   
    ENTERS("DataProjection_ResultMaximumErrorGet",err,error,*999)

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionFinished) CALL FlagError("Data projection has not been finished.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionProjected) CALL FlagError("Data projection has not been projected.",err,error,*999)
    
    maximumDataPoint=dataProjection%maximumErrorDataPoint
    maximumError=dataProjection%maximumError
 
    EXITS("DataProjection_ResultMaximumErrorGet")
    RETURN
999 ERRORSEXITS("DataProjection_ResultMaximumErrorGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ResultMaximumErrorGet

  !
  !================================================================================================================================
  !

  !>Gets the projection minimum error for a data projection.
  SUBROUTINE DataProjection_ResultMinimumErrorGet(dataProjection,minimumDataPoint,minimumError,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to get the minimum error for
    INTEGER(INTG), INTENT(OUT) :: minimumDataPoint !<On exit, the data point number of the minimum error.
    REAL(DP), INTENT(OUT) :: minimumError !<On exit, the maximum error for the data projection.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
   
    ENTERS("DataProjection_ResultMinimumErrorGet",err,error,*999)

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionFinished) CALL FlagError("Data projection has not been finished.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionProjected) CALL FlagError("Data projection has not been projected.",err,error,*999)
    
    minimumDataPoint=dataProjection%minimumErrorDataPoint
    minimumError=dataProjection%minimumError
 
    EXITS("DataProjection_ResultMinimumErrorGet")
    RETURN
999 ERRORSEXITS("DataProjection_ResultMinimumErrorGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ResultMinimumErrorGet

  !
  !================================================================================================================================
  !

  !>Gets the projection Root Mean Squared (RMS) error for a data projection.
  SUBROUTINE DataProjection_ResultRMSErrorGet(dataProjection,rmsError,err,error,*)

    !Argument variables
    TYPE(DataProjectionType), POINTER :: dataProjection !<A pointer to the data projection to get the RMS error for
    REAL(DP), INTENT(OUT) :: rmsError !<On exit, the RMS error for the data projection.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("DataProjection_ResultRMSErrorGet",err,error,*999)

    IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionFinished) CALL FlagError("Data projection has not been finished.",err,error,*999)
    IF(.NOT.dataProjection%dataProjectionProjected) CALL FlagError("Data projection has not been projected.",err,error,*999)
    
    rmsError=dataProjection%rmsError
 
    EXITS("DataProjection_ResultRMSErrorGet")
    RETURN
999 ERRORSEXITS("DataProjection_ResultRMSErrorGet",err,error)    
    RETURN 1

  END SUBROUTINE DataProjection_ResultRMSErrorGet

  !
  !================================================================================================================================
  !

END MODULE DataProjectionAccessRoutines
