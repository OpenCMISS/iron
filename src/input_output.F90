!> \file
!> \author Chris Bradley
!> \brief This module handles all formating and input and output.
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

!> This module handles all formating and input and output.
MODULE InputOutput

  USE BaseRoutines
  USE Constants
  USE Kinds
  USE ISO_VARYING_STRING
  USE Strings

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE
  
  !Module parameters

  !> \addtogroup InputOutput_MatrixNameIndexFormat InputOutput::MatrixNameIndexFormat
  !> \brief Output type parameter
  !> \see InputOutput
  !>@{  
  INTEGER(INTG), PARAMETER :: WRITE_STRING_MATRIX_NAME_ONLY=1 !<Write the matrix name with out any indices \see InputOutput_MatrixNameIndexFormat,InputOutput::MatrixNameIndexFormat
  INTEGER(INTG), PARAMETER :: WRITE_STRING_MATRIX_NAME_AND_INDICES=2 !<Write the matrix name together with the matrix indices \see InputOutput_MatrixNameIndexFormat,InputOutput::MatrixNameIndexFormat
  !>@}

  !Module types

  !Interfaces

  !>Write a string to a given output stream
  INTERFACE WRITE_STRING
    MODULE PROCEDURE WriteStringC
    MODULE PROCEDURE WriteStringVS
  END INTERFACE WRITE_STRING

  !>Write a string to a given output stream
  INTERFACE WriteString
    MODULE PROCEDURE WriteStringC
    MODULE PROCEDURE WriteStringVS
  END INTERFACE WriteString

  !>Write a string followed by a value to a given output stream
  INTERFACE WRITE_STRING_VALUE
    MODULE PROCEDURE WriteStringValueC
    MODULE PROCEDURE WriteStringValueDP
    MODULE PROCEDURE WriteStringValueIntg
    MODULE PROCEDURE WriteStringValueLIntg
    MODULE PROCEDURE WriteStringValueL
    MODULE PROCEDURE WriteStringValueSP
    MODULE PROCEDURE WriteStringValueVS
  END INTERFACE WRITE_STRING_VALUE

  !>Write a string followed by a value to a given output stream
  INTERFACE WriteStringValue
    MODULE PROCEDURE WriteStringValueC
    MODULE PROCEDURE WriteStringValueDP
    MODULE PROCEDURE WriteStringValueIntg
    MODULE PROCEDURE WriteStringValueLIntg
    MODULE PROCEDURE WriteStringValueL
    MODULE PROCEDURE WriteStringValueSP
    MODULE PROCEDURE WriteStringValueVS
  END INTERFACE WriteStringValue

  !>Write a string, value, string then a value to a given output stream
  INTERFACE WRITE_STRING_TWO_VALUE
    MODULE PROCEDURE WriteStringTwoValueCC
    MODULE PROCEDURE WriteStringTwoValueCDP
    MODULE PROCEDURE WriteStringTwoValueCIntg
    MODULE PROCEDURE WriteStringTwoValueCL
    MODULE PROCEDURE WriteStringTwoValueCSP
    MODULE PROCEDURE WriteStringTwoValueCVS
    MODULE PROCEDURE WriteStringTwoValueDPC
    MODULE PROCEDURE WriteStringTwoValueDPDP
    MODULE PROCEDURE WriteStringTwoValueDPIntg
    MODULE PROCEDURE WriteStringTwoValueDPL
    MODULE PROCEDURE WriteStringTwoValueDPSP
    MODULE PROCEDURE WriteStringTwoValueDPVS
    MODULE PROCEDURE WriteStringTwoValueIntgC
    MODULE PROCEDURE WriteStringTwoValueIntgDP
    MODULE PROCEDURE WriteStringTwoValueIntgIntg
    MODULE PROCEDURE WriteStringTwoValueIntgL
    MODULE PROCEDURE WriteStringTwoValueIntgSP
    MODULE PROCEDURE WriteStringTwoValueIntgVS
    MODULE PROCEDURE WriteStringTwoValueLC
    MODULE PROCEDURE WriteStringTwoValueLDP
    MODULE PROCEDURE WriteStringTwoValueLIntg
    MODULE PROCEDURE WriteStringTwoValueLL
    MODULE PROCEDURE WriteStringTwoValueLSP
    MODULE PROCEDURE WriteStringTwoValueLVS
    MODULE PROCEDURE WriteStringTwoValueSPC
    MODULE PROCEDURE WriteStringTwoValueSPDP
    MODULE PROCEDURE WriteStringTwoValueSPIntg
    MODULE PROCEDURE WriteStringTwoValueSPL
    MODULE PROCEDURE WriteStringTwoValueSPSP
    MODULE PROCEDURE WriteStringTwoValueSPVS
    MODULE PROCEDURE WriteStringTwoValueVSC
    MODULE PROCEDURE WriteStringTwoValueVSDP
    MODULE PROCEDURE WriteStringTwoValueVSIntg
    MODULE PROCEDURE WriteStringTwoValueVSL
    MODULE PROCEDURE WriteStringTwoValueVSSP
    MODULE PROCEDURE WriteStringTwoValueVSVS
  END INTERFACE WRITE_STRING_TWO_VALUE

  !>Write a string, value, string then a value to a given output stream
  INTERFACE WriteStringTwoValue
    MODULE PROCEDURE WriteStringTwoValueCC
    MODULE PROCEDURE WriteStringTwoValueCDP
    MODULE PROCEDURE WriteStringTwoValueCIntg
    MODULE PROCEDURE WriteStringTwoValueCL
    MODULE PROCEDURE WriteStringTwoValueCSP
    MODULE PROCEDURE WriteStringTwoValueCVS
    MODULE PROCEDURE WriteStringTwoValueDPC
    MODULE PROCEDURE WriteStringTwoValueDPDP
    MODULE PROCEDURE WriteStringTwoValueDPIntg
    MODULE PROCEDURE WriteStringTwoValueDPL
    MODULE PROCEDURE WriteStringTwoValueDPSP
    MODULE PROCEDURE WriteStringTwoValueDPVS
    MODULE PROCEDURE WriteStringTwoValueIntgC
    MODULE PROCEDURE WriteStringTwoValueIntgDP
    MODULE PROCEDURE WriteStringTwoValueIntgIntg
    MODULE PROCEDURE WriteStringTwoValueIntgL
    MODULE PROCEDURE WriteStringTwoValueIntgSP
    MODULE PROCEDURE WriteStringTwoValueIntgVS
    MODULE PROCEDURE WriteStringTwoValueLC
    MODULE PROCEDURE WriteStringTwoValueLDP
    MODULE PROCEDURE WriteStringTwoValueLIntg
    MODULE PROCEDURE WriteStringTwoValueLL
    MODULE PROCEDURE WriteStringTwoValueLSP
    MODULE PROCEDURE WriteStringTwoValueLVS
    MODULE PROCEDURE WriteStringTwoValueSPC
    MODULE PROCEDURE WriteStringTwoValueSPDP
    MODULE PROCEDURE WriteStringTwoValueSPIntg
    MODULE PROCEDURE WriteStringTwoValueSPL
    MODULE PROCEDURE WriteStringTwoValueSPSP
    MODULE PROCEDURE WriteStringTwoValueSPVS
    MODULE PROCEDURE WriteStringTwoValueVSC
    MODULE PROCEDURE WriteStringTwoValueVSDP
    MODULE PROCEDURE WriteStringTwoValueVSIntg
    MODULE PROCEDURE WriteStringTwoValueVSL
    MODULE PROCEDURE WriteStringTwoValueVSSP
    MODULE PROCEDURE WriteStringTwoValueVSVS
  END INTERFACE WriteStringTwoValue

  !>Write a string followed by a value formatted in a particular way to a specified output stream
  INTERFACE WRITE_STRING_FMT_VALUE
    MODULE PROCEDURE WriteStringFmtValueC
    MODULE PROCEDURE WriteStringFmtValueDP
    MODULE PROCEDURE WriteStringFmtValueIntg
    MODULE PROCEDURE WriteStringFmtValueLIntg
    MODULE PROCEDURE WriteStringFmtValueL
    MODULE PROCEDURE WriteStringFmtValueSP
    MODULE PROCEDURE WriteStringFmtValueVS
  END INTERFACE WRITE_STRING_FMT_VALUE
  
  !>Write a string followed by a value formatted in a particular way to a specified output stream
  INTERFACE WriteStringFmtValue
    MODULE PROCEDURE WriteStringFmtValueC
    MODULE PROCEDURE WriteStringFmtValueDP
    MODULE PROCEDURE WriteStringFmtValueIntg
    MODULE PROCEDURE WriteStringFmtValueLIntg
    MODULE PROCEDURE WriteStringFmtValueL
    MODULE PROCEDURE WriteStringFmtValueSP
    MODULE PROCEDURE WriteStringFmtValueVS
  END INTERFACE WriteStringFmtValue
  
  !>Write a string, value, string then a value with the values formatted in a particular way to a given output stream
  INTERFACE WRITE_STRING_FMT_TWO_VALUE
    MODULE PROCEDURE WriteStringFmtTwoValueCC
    MODULE PROCEDURE WriteStringFmtTwoValueCDP
    MODULE PROCEDURE WriteStringFmtTwoValueCIntg
    MODULE PROCEDURE WriteStringFmtTwoValueCL
    MODULE PROCEDURE WriteStringFmtTwoValueCSP
    MODULE PROCEDURE WriteStringFmtTwoValueCVS
    MODULE PROCEDURE WriteStringFmtTwoValueDPC
    MODULE PROCEDURE WriteStringFmtTwoValueDPDP
    MODULE PROCEDURE WriteStringFmtTwoValueDPIntg
    MODULE PROCEDURE WriteStringFmtTwoValueDPL
    MODULE PROCEDURE WriteStringFmtTwoValueDPSP
    MODULE PROCEDURE WriteStringFmtTwoValueDPVS
    MODULE PROCEDURE WriteStringFmtTwoValueIntgC
    MODULE PROCEDURE WriteStringFmtTwoValueIntgDP
    MODULE PROCEDURE WriteStringFmtTwoValueIntgIntg
    MODULE PROCEDURE WriteStringFmtTwoValueIntgL
    MODULE PROCEDURE WriteStringFmtTwoValueIntgSP
    MODULE PROCEDURE WriteStringFmtTwoValueIntgVS
    MODULE PROCEDURE WriteStringFmtTwoValueLC
    MODULE PROCEDURE WriteStringFmtTwoValueLDP
    MODULE PROCEDURE WriteStringFmtTwoValueLIntg
    MODULE PROCEDURE WriteStringFmtTwoValueLL
    MODULE PROCEDURE WriteStringFmtTwoValueLSP
    MODULE PROCEDURE WriteStringFmtTwoValueLVS
    MODULE PROCEDURE WriteStringFmtTwoValueSPC
    MODULE PROCEDURE WriteStringFmtTwoValueSPDP
    MODULE PROCEDURE WriteStringFmtTwoValueSPIntg
    MODULE PROCEDURE WriteStringFmtTwoValueSPL
    MODULE PROCEDURE WriteStringFmtTwoValueSPSP
    MODULE PROCEDURE WriteStringFmtTwoValueSPVS
    MODULE PROCEDURE WriteStringFmtTwoValueVSC
    MODULE PROCEDURE WriteStringFmtTwoValueVSDP
    MODULE PROCEDURE WriteStringFmtTwoValueVSIntg
    MODULE PROCEDURE WriteStringFmtTwoValueVSL
    MODULE PROCEDURE WriteStringFmtTwoValueVSSP
    MODULE PROCEDURE WriteStringFmtTwoValueVSVS
  END INTERFACE WRITE_STRING_FMT_TWO_VALUE

  !>Write a string, value, string then a value with the values formatted in a particular way to a given output stream
  INTERFACE WriteStringFmtTwoValue
    MODULE PROCEDURE WriteStringFmtTwoValueCC
    MODULE PROCEDURE WriteStringFmtTwoValueCDP
    MODULE PROCEDURE WriteStringFmtTwoValueCIntg
    MODULE PROCEDURE WriteStringFmtTwoValueCL
    MODULE PROCEDURE WriteStringFmtTwoValueCSP
    MODULE PROCEDURE WriteStringFmtTwoValueCVS
    MODULE PROCEDURE WriteStringFmtTwoValueDPC
    MODULE PROCEDURE WriteStringFmtTwoValueDPDP
    MODULE PROCEDURE WriteStringFmtTwoValueDPIntg
    MODULE PROCEDURE WriteStringFmtTwoValueDPL
    MODULE PROCEDURE WriteStringFmtTwoValueDPSP
    MODULE PROCEDURE WriteStringFmtTwoValueDPVS
    MODULE PROCEDURE WriteStringFmtTwoValueIntgC
    MODULE PROCEDURE WriteStringFmtTwoValueIntgDP
    MODULE PROCEDURE WriteStringFmtTwoValueIntgIntg
    MODULE PROCEDURE WriteStringFmtTwoValueIntgL
    MODULE PROCEDURE WriteStringFmtTwoValueIntgSP
    MODULE PROCEDURE WriteStringFmtTwoValueIntgVS
    MODULE PROCEDURE WriteStringFmtTwoValueLC
    MODULE PROCEDURE WriteStringFmtTwoValueLDP
    MODULE PROCEDURE WriteStringFmtTwoValueLIntg
    MODULE PROCEDURE WriteStringFmtTwoValueLL
    MODULE PROCEDURE WriteStringFmtTwoValueLSP
    MODULE PROCEDURE WriteStringFmtTwoValueLVS
    MODULE PROCEDURE WriteStringFmtTwoValueSPC
    MODULE PROCEDURE WriteStringFmtTwoValueSPDP
    MODULE PROCEDURE WriteStringFmtTwoValueSPIntg
    MODULE PROCEDURE WriteStringFmtTwoValueSPL
    MODULE PROCEDURE WriteStringFmtTwoValueSPSP
    MODULE PROCEDURE WriteStringFmtTwoValueSPVS
    MODULE PROCEDURE WriteStringFmtTwoValueVSC
    MODULE PROCEDURE WriteStringFmtTwoValueVSDP
    MODULE PROCEDURE WriteStringFmtTwoValueVSIntg
    MODULE PROCEDURE WriteStringFmtTwoValueVSL
    MODULE PROCEDURE WriteStringFmtTwoValueVSSP
    MODULE PROCEDURE WriteStringFmtTwoValueVSVS
  END INTERFACE WriteStringFmtTwoValue

  !>Write a string followed by a vector to a specified output stream.
  INTERFACE WRITE_STRING_VECTOR
    MODULE PROCEDURE WriteStringVectorDP
    MODULE PROCEDURE WriteStringVectorIntg
    MODULE PROCEDURE WriteStringVectorLIntg
    MODULE PROCEDURE WriteStringVectorL
    MODULE PROCEDURE WriteStringVectorSP
  END INTERFACE WRITE_STRING_VECTOR
  
  !>Write a string followed by a vector to a specified output stream.
  INTERFACE WriteStringVector
    MODULE PROCEDURE WriteStringVectorDP
    MODULE PROCEDURE WriteStringVectorIntg
    MODULE PROCEDURE WriteStringVectorLIntg
    MODULE PROCEDURE WriteStringVectorL
    MODULE PROCEDURE WriteStringVectorSP
  END INTERFACE WriteStringVector
  
  !>Write a string followed by a indexed vector to a specified output stream.
  INTERFACE WRITE_STRING_IDX_VECTOR
    MODULE PROCEDURE WriteStringIdxVectorDP
    MODULE PROCEDURE WriteStringIdxVectorIntg
    MODULE PROCEDURE WriteStringIdxVectorLIntg
    MODULE PROCEDURE WriteStringIdxVectorL
    MODULE PROCEDURE WriteStringIdxVectorSP
  END INTERFACE WRITE_STRING_IDX_VECTOR

  !>Write a string followed by a indexed vector to a specified output stream.
  INTERFACE WriteStringIdxVector
    MODULE PROCEDURE WriteStringIdxVectorDP
    MODULE PROCEDURE WriteStringIdxVectorIntg
    MODULE PROCEDURE WriteStringIdxVectorLIntg
    MODULE PROCEDURE WriteStringIdxVectorL
    MODULE PROCEDURE WriteStringIdxVectorSP
  END INTERFACE WriteStringIdxVector

  !>Write a string followed by a matrix to a specified output stream
  INTERFACE WRITE_STRING_MATRIX
    MODULE PROCEDURE WriteStringMatrixDP
    MODULE PROCEDURE WriteStringMatrixIntg
    MODULE PROCEDURE WriteStringMatrixLIntg
    MODULE PROCEDURE WriteStringMatrixL
    MODULE PROCEDURE WriteStringMatrixSP
  END INTERFACE WRITE_STRING_MATRIX

  !>Write a string followed by a matrix to a specified output stream
  INTERFACE WriteStringMatrix
    MODULE PROCEDURE WriteStringMatrixDP
    MODULE PROCEDURE WriteStringMatrixIntg
    MODULE PROCEDURE WriteStringMatrixLIntg
    MODULE PROCEDURE WriteStringMatrixL
    MODULE PROCEDURE WriteStringMatrixSP
  END INTERFACE WriteStringMatrix

  PUBLIC WRITE_STRING_MATRIX_NAME_ONLY,WRITE_STRING_MATRIX_NAME_AND_INDICES
  
  PUBLIC WRITE_STRING,WRITE_STRING_VALUE,WRITE_STRING_TWO_VALUE,WRITE_STRING_FMT_VALUE,WRITE_STRING_FMT_TWO_VALUE, &
    & WRITE_STRING_VECTOR,WRITE_STRING_IDX_VECTOR,WRITE_STRING_MATRIX

  PUBLIC WriteString

  PUBLIC WriteStringValue

  PUBLIC WriteStringTwoValue

  PUBLIC WriteStringFmtValue

  PUBLIC WriteStringFmtTwoValue

  PUBLIC WriteStringVector
  
  PUBLIC WriteStringIdxVector

  PUBLIC WriteStringMatrix

  !Module variables

CONTAINS

  !!TODO: put back enters,exits etc.
 
  !
  !================================================================================================================================
  !

  !>Writes the character string to the given output stream specified by id.
  SUBROUTINE WriteStringC(id,string,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream to write to \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: string !<The string to write
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
        
    WRITE(outputString,'(A)') string
    CALL WriteStr(id,err,error,*999)
      
    RETURN
999 ERRORS("WriteStringC",err,error)
    RETURN 1
    
  END SUBROUTINE WriteStringC

  !
  !================================================================================================================================
  !

  !>Writes the varying string to the given output stream specified by id.
  SUBROUTINE WriteStringVS(id,string,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream to write to \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    TYPE(VARYING_STRING), INTENT(IN) :: string !<The string to write
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
        
    WRITE(outputString,'(A)') CHAR(string)
    CALL WriteStr(id,err,error,*999)
      
    RETURN
999 ERRORS("WriteStringVS",err,error)
    RETURN 1
    
  END SUBROUTINE WriteStringVS

  !
  !================================================================================================================================
  !

  !>Writes the first string followed by a formatted character value to the given output stream specified by id. Free format is used to format the value.
  SUBROUTINE WriteStringValueC(id,firstString,value,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    CHARACTER(LEN=*), INTENT(IN) :: value !<The value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString
        
    localString=firstString//value(1:LEN_TRIM(value))
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
    RETURN
999 ERRORS("WriteStringValueC",err,error)
    RETURN 1
    
  END SUBROUTINE WriteStringValueC

  !
  !================================================================================================================================
  !

  !>Writes the first string followed by a formatted double precision value to the given output stream specified by ID. Free format is used to format the value.
  SUBROUTINE WriteStringValueDP(id,firstString,value,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    REAL(DP), INTENT(IN) :: value !<The value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    CHARACTER(LEN=1) :: formatString = "*"
    TYPE(VARYING_STRING) :: localString
       
    localString=firstString//NumberToVString(value,formatString,err,error)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
    RETURN
999 ERRORS("WriteStringValueDP",err,error)
    RETURN 1
    
  END SUBROUTINE WriteStringValueDP

  !
  !================================================================================================================================
  !

  !>Writes the first string followed by a formatted integer value to the given output stream specified by ID. Free format is used to format the value.
  SUBROUTINE WriteStringValueIntg(id,firstString,value,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    INTEGER(INTG), INTENT(IN) :: value !<The value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString
        
    localString=firstString//NumberToVString(value,"*",err,error)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
    RETURN
999 ERRORS("WriteStringValueIntg",err,error)
    RETURN 1
    
  END SUBROUTINE WriteStringValueIntg

  !
  !================================================================================================================================
  !

  !>Writes the first string followed by a formatted long integer value to the given output stream specified by ID. Free format is used to format the value.
  SUBROUTINE WriteStringValueLIntg(id,firstString,VALUE,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    INTEGER(LINTG), INTENT(IN) :: value !<The value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

    localString=firstString//NumberToVString(value,"*",err,error)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
    RETURN
999 ERRORS("WriteStringValueLIntg",err,error)
    RETURN 1
    
  END SUBROUTINE WriteStringValueLIntg

  !
  !================================================================================================================================
  !

  !>Writes the first string followed by a formatted logical value to the given output stream specified by ID. Free format is used to format the value.
  SUBROUTINE WriteStringValueL(id,firstString,value,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    LOGICAL, INTENT(IN) :: value !<The value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

    localString=firstString//LogicalToVString(value,err,error)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
    RETURN
999 ERRORS("WriteStringValueL",err,error)
    RETURN 1
    
  END SUBROUTINE WriteStringValueL

  !
  !================================================================================================================================
  !

  !>Writes the first string followed by a formatted single precision value to the given output stream specified by ID. Free format is used to format the value.
  SUBROUTINE WriteStringValueSP(id,firstString,VALUE,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    REAL(SP), INTENT(IN) :: value !<The value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    CHARACTER(LEN=1) :: formatString = "*"
    TYPE(VARYING_STRING) :: localString
       
    localString=firstString//NumberToVString(value,formatString,err,error)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
    RETURN
999 ERRORS("WriteStringValueSP",err,error)
    RETURN 1
    
  END SUBROUTINE WriteStringValueSP

  !
  !================================================================================================================================
  !

  !>Writes the first string followed by a formatted varying string value to the given output stream specified by ID. Free format is used to format the value.
  SUBROUTINE WriteStringValueVS(id,firstString,value,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: value !<The value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

    IF(value==" ") THEN
      localString=firstString
    ELSE
      localString=firstString//value
    ENDIF
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)

    RETURN
999 ERRORS("WriteStringValueVS",err,error)
    RETURN 1
    
  END SUBROUTINE WriteStringValueVS

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted character firstValue and the the secondString followed by a formatted character secondValue to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WriteStringTwoValueCC(id,firstString,firstValue,secondString,secondValue,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondValue !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringTwoValueCC",err,error,*999)
        
    localString=firstString//firstValue//secondString//secondValue
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringTwoValueCC")
    RETURN
999 ERRORS("WriteStringTwoValueCC",err,error)
!    EXITS("WriteStringTwoValueCC")
    RETURN 1
    
  END SUBROUTINE WriteStringTwoValueCC

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted character firstValue and the the secondString followed by a formatted double precision secondValue to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WriteStringTwoValueCDP(id,firstString,firstValue,secondString,secondValue,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    REAL(DP), INTENT(IN) :: secondValue !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

 !   ENTERS("WriteStringTwoValueCDP",err,error,*999)
        
    localString=firstString//firstValue//secondString//NumberToVString(secondValue,"*",err,error)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringTwoValueCDP")
    RETURN
999 ERRORS("WriteStringTwoValueCDP",err,error)
!    EXITS("WriteStringTwoValueCDP")
    RETURN 1
    
  END SUBROUTINE WriteStringTwoValueCDP

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted character firstValue and the the secondString followed by a formatted integer secondValue to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WriteStringTwoValueCIntg(id,firstString,firstValue,secondString,secondValue,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    INTEGER(INTG), INTENT(IN) :: secondValue !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringTwoValueCIntg",err,error,*999)
        
    localString=firstString//firstValue//secondString//NumberToVString(secondValue,"*",err,error)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringTwoValueCIntg")
    RETURN
999 ERRORS("WriteStringTwoValueCIntg",err,error)
!    EXITS("WriteStringTwoValueCIntg")
    RETURN 1
    
  END SUBROUTINE WriteStringTwoValueCIntg

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted character firstValue and the the secondString followed by a formatted logical secondValue to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WriteStringTwoValueCL(id,firstString,firstValue,secondString,secondValue,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    LOGICAL, INTENT(IN) :: secondValue !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

    !ENTERS("WriteStringTwoValueCL",err,error,*999)
        
    localString=firstString//firstValue//secondString//LogicalToVString(secondValue,err,error)
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
    !EXITS("WriteStringTwoValueCL")
    RETURN
999 ERRORS("WriteStringTwoValueCL",err,error)
    !EXITS("WriteStringTwoValueCL")
    
    RETURN 1   
  END SUBROUTINE WriteStringTwoValueCL

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted character firstValue and the the secondString followed by a formatted single precision secondValue to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WriteStringTwoValueCSP(id,firstString,firstValue,secondString,secondValue,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    REAL(SP), INTENT(IN) :: secondValue !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

 !   ENTERS("WriteStringTwoValueCSP",err,error,*999)
        
    localString=firstString//firstValue//secondString//NumberToVString(secondValue,"*",err,error)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringTwoValueCSP")
    RETURN
999 ERRORS("WriteStringTwoValueCSP",err,error)
!    EXITS("WriteStringTwoValueCSP")
    RETURN 1
    
  END SUBROUTINE WriteStringTwoValueCSP

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted character firstValue and the the secondString followed by a formatted varying string secondValue to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WriteStringTwoValueCVS(id,firstString,firstValue,secondString,secondValue,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: secondValue !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringTwoValueCVS",err,error,*999)
        
    localString=firstString//firstValue//secondString//secondValue
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringTwoValueCVS")
    RETURN
999 ERRORS("WriteStringTwoValueCVS",err,error)
!    EXITS("WriteStringTwoValueCVS")
    RETURN 1
    
  END SUBROUTINE WriteStringTwoValueCVS

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted double precision firstValue and the the secondString followed by a formatted character secondValue to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WriteStringTwoValueDPC(id,firstString,firstValue,secondString,secondValue,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    REAL(DP), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondValue !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringTwoValueDPC",err,error,*999)
        
    localString=firstString//NumberToVString(firstValue,"*",err,error)//secondString//secondValue
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringTwoValueDPC")
    RETURN
999 ERRORS("WriteStringTwoValueDPC",err,error)
!    EXITS("WriteStringTwoValueDPC")
    RETURN 1
    
  END SUBROUTINE WriteStringTwoValueDPC

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted double precision firstValue and the the secondString followed by a formatted double precision secondValue to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WriteStringTwoValueDPDP(id,firstString,firstValue,secondString,secondValue,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    REAL(DP), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    REAL(DP), INTENT(IN) :: secondValue !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString,localString2

 !   ENTERS("WriteStringTwoValueDPDP",err,error,*999)
        
    localString=firstString//NumberToVString(firstValue,"*",err,error)
    IF(err/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    localString2=localString//secondString
    localString=localString2//NumberToVString(secondValue,"*",err,error)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringTwoValueDPDP")
    RETURN
999 ERRORS("WriteStringTwoValueDPDP",err,error)
!    EXITS("WriteStringTwoValueDPDP")
    RETURN 1
    
  END SUBROUTINE WriteStringTwoValueDPDP

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted double precision firstValue and the the secondString followed by a formatted integer secondValue to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WriteStringTwoValueDPIntg(id,firstString,firstValue,secondString,secondValue,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    REAL(DP), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    INTEGER(INTG), INTENT(IN) :: secondValue !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString,localString2

!    ENTERS("WriteStringTwoValueDPIntg",err,error,*999)
        
    localString=firstString//NumberToVString(firstValue,"*",err,error)
    IF(err/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    localString2=localString//secondString
    localString=localString2//NumberToVString(secondValue,"*",err,error)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringTwoValueDPIntg")
    RETURN
999 ERRORS("WriteStringTwoValueDPIntg",err,error)
!    EXITS("WriteStringTwoValueDPIntg")
    RETURN 1
    
  END SUBROUTINE WriteStringTwoValueDPIntg

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted double precision firstValue and the the secondString followed by a formatted logical secondValue to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WriteStringTwoValueDPL(id,firstString,firstValue,secondString,secondValue,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits7
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    REAL(DP), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    LOGICAL, INTENT(IN) :: secondValue !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString,localString2

!    ENTERS("WriteStringTwoValueDPL",err,error,*999)
        
    localString=firstString//NumberToVString(firstValue,"*",err,error)
    IF(err/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    localString2=localString//secondString
    localString=localString2//LogicalToVString(secondValue,err,error)
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringTwoValueDPL")
    RETURN
999 ERRORS("WriteStringTwoValueDPL",err,error)
!    EXITS("WriteStringTwoValueDPL")
    RETURN 1
    
  END SUBROUTINE WriteStringTwoValueDPL

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted double precision firstValue and the the secondString followed by a formatted single precision secondValue to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WriteStringTwoValueDPSP(id,firstString,firstValue,secondString,secondValue,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    REAL(DP), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    REAL(SP), INTENT(IN) :: secondValue !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString,localString2

 !   ENTERS("WriteStringTwoValueDPSP",err,error,*999)
        
    localString=firstString//NumberToVString(firstValue,"*",err,error)
    IF(err/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    localString2=localString//secondString
    localString=localString2//NumberToVString(secondValue,"*",err,error)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringTwoValueDPSP")
    RETURN
999 ERRORS("WriteStringTwoValueDPSP",err,error)
!    EXITS("WriteStringTwoValueDPSP")
    RETURN 1
    
  END SUBROUTINE WriteStringTwoValueDPSP

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted double precision firstValue and the the secondString followed by a formatted single precision secondValue to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WriteStringTwoValueDPVS(id,firstString,firstValue,secondString,secondValue,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    REAL(DP), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: secondValue !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringTwoValueDPVS",err,error,*999)
        
    localString=firstString//NumberToVString(firstValue,"*",err,error)//secondString//secondValue
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringTwoValueDPVS")
    RETURN
999 ERRORS("WriteStringTwoValueDPVS",err,error)
!    EXITS("WriteStringTwoValueDPVS")
    RETURN 1
    
  END SUBROUTINE WriteStringTwoValueDPVS

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted integer firstValue and the the secondString followed by a formatted character secondValue to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WriteStringTwoValueIntgC(id,firstString,firstValue,secondString,secondValue,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    INTEGER(INTG), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondValue !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringTwoValueIntgC",err,error,*999)
        
    localString=firstString//NumberToVString(firstValue,"*",err,error)//secondString//secondValue
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringTwoValueIntgC")
    RETURN
999 ERRORS("WriteStringTwoValueIntgC",err,error)
!    EXITS("WriteStringTwoValueIntgC")
    RETURN 1
    
  END SUBROUTINE WriteStringTwoValueIntgC

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted integer firstValue and the the secondString followed by a formatted double precision secondValue to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WriteStringTwoValueIntgDP(id,firstString,firstValue,secondString,secondValue,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    INTEGER(INTG), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    REAL(DP), INTENT(IN) :: secondValue !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString,localString2

 !   ENTERS("WriteStringTwoValueIntgDP",err,error,*999)
        
    localString=firstString//NumberToVString(firstValue,"*",err,error)
    IF(err/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    localString2=localString//secondString
    localString=localString2//NumberToVString(secondValue,"*",err,error)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringTwoValueIntgDP")
    RETURN
999 ERRORS("WriteStringTwoValueIntgDP",err,error)
!    EXITS("WriteStringTwoValueIntgDP")
    RETURN 1
    
  END SUBROUTINE WriteStringTwoValueIntgDP

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted integer firstValue and the the secondString followed by a formatted integer secondValue to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WriteStringTwoValueIntgIntg(id,firstString,firstValue,secondString,secondValue,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    INTEGER(INTG), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    INTEGER(INTG), INTENT(IN) :: secondValue !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString,localString2

!    ENTERS("WriteStringTwoValueIntgIntg",err,error,*999)
        
    localString=firstString//NumberToVString(firstValue,"*",err,error)
    IF(err/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    localString2=localString//secondString
    localString=localString2//NumberToVString(secondValue,"*",err,error)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringTwoValueIntgIntg")
    RETURN
999 ERRORS("WriteStringTwoValueIntgIntg",err,error)
!    EXITS("WriteStringTwoValueIntgIntg")
    RETURN 1
    
  END SUBROUTINE WriteStringTwoValueIntgIntg

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted integer firstValue and the the secondString followed by a formatted logical secondValue to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WriteStringTwoValueIntgL(id,firstString,firstValue,secondString,secondValue,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    INTEGER(INTG), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    LOGICAL, INTENT(IN) :: secondValue !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString,localString2

!    ENTERS("WriteStringTwoValueIntgL",err,error,*999)
        
    localString=firstString//NumberToVString(firstValue,"*",err,error)
    IF(err/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    localString2=localString//secondString
    localString=localString2//LogicalToVString(secondValue,err,error)
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringTwoValueIntgL")
    RETURN
999 ERRORS("WriteStringTwoValueIntgL",err,error)
!    EXITS("WriteStringTwoValueIntgL")
    RETURN 1
    
  END SUBROUTINE WriteStringTwoValueIntgL

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted integer firstValue and the the secondString followed by a formatted single precision secondValue to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WriteStringTwoValueIntgSP(id,firstString,firstValue,secondString,secondValue,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    INTEGER(INTG), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    REAL(SP), INTENT(IN) :: secondValue !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString,localString2

 !   ENTERS("WriteStringTwoValueIntgSP",err,error,*999)
        
    localString=firstString//NumberToVString(firstValue,"*",err,error)
    IF(err/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    localString2=localString//secondString
    localString=localString2//NumberToVString(secondValue,"*",err,error)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringTwoValueIntgSP")
    RETURN
999 ERRORS("WriteStringTwoValueIntgSP",err,error)
!    EXITS("WriteStringTwoValueIntgSP")
    RETURN 1
    
  END SUBROUTINE WriteStringTwoValueIntgSP

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted integer firstValue and the the secondString followed by a formatted single precision secondValue to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WriteStringTwoValueIntgVS(id,firstString,firstValue,secondString,secondValue,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    INTEGER(INTG), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: secondValue !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringTwoValueIntgVS",err,error,*999)
        
    localString=firstString//NumberToVString(firstValue,"*",err,error)//secondString//secondValue
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringTwoValueIntgVS")
    RETURN
999 ERRORS("WriteStringTwoValueIntgVS",err,error)
!    EXITS("WriteStringTwoValueIntgVS")
    RETURN 1
    
  END SUBROUTINE WriteStringTwoValueIntgVS

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted logical firstValue and the the secondString followed by a formatted character secondValue to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WriteStringTwoValueLC(id,firstString,firstValue,secondString,secondValue,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    LOGICAL, INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondValue !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringTwoValueLC",err,error,*999)
        
    localString=firstString//LogicalToVString(firstValue,err,error)//secondString//secondValue
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringTwoValueLC")
    RETURN
999 ERRORS("WriteStringTwoValueLC",err,error)
!    EXITS("WriteStringTwoValueLC")
    RETURN 1
    
  END SUBROUTINE WriteStringTwoValueLC

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted logical firstValue and the the secondString followed by a formatted double precision secondValue to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WriteStringTwoValueLDP(id,firstString,firstValue,secondString,secondValue,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    LOGICAL, INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    REAL(DP), INTENT(IN) :: secondValue !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString,localString2

 !   ENTERS("WriteStringTwoValueLDP",err,error,*999)
        
    localString=firstString//LogicalToVString(firstValue,err,error)
    IF(err/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    localString2=localString//secondString
    localString=localString2//NumberToVString(secondValue,"*",err,error)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringTwoValueLDP")
    RETURN
999 ERRORS("WriteStringTwoValueLDP",err,error)
!    EXITS("WriteStringTwoValueLDP")
    RETURN 1
    
  END SUBROUTINE WriteStringTwoValueLDP

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted logical firstValue and the the secondString followed by a formatted integer secondValue to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WriteStringTwoValueLIntg(id,firstString,firstValue,secondString,secondValue,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    LOGICAL, INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    INTEGER(INTG), INTENT(IN) :: secondValue !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringTwoValueLIntg",err,error,*999)
        
    localString=firstString//LogicalToVString(firstValue,err,error)
    IF(err/=0) GOTO 999
    localString=localString//secondString//NumberToVString(secondValue,"*",err,error)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringTwoValueLIntg")
    RETURN
999 ERRORS("WriteStringTwoValueLIntg",err,error)
!    EXITS("WriteStringTwoValueLIntg")
    RETURN 1
    
  END SUBROUTINE WriteStringTwoValueLIntg

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted logical firstValue and the the secondString followed by a formatted logical secondValue to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WriteStringTwoValueLL(id,firstString,firstValue,secondString,secondValue,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    LOGICAL, INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    LOGICAL, INTENT(IN) :: secondValue !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString,localString2

!    ENTERS("WriteStringTwoValueLL",err,error,*999)
        
    localString=firstString//LogicalToVString(firstValue,err,error)
    IF(err/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    localString2=localString//secondString
    localString=localString2//LogicalToVString(secondValue,err,error)
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringTwoValueLL")
    RETURN
999 ERRORS("WriteStringTwoValueLL",err,error)
!    EXITS("WriteStringTwoValueLL")
    RETURN 1
    
  END SUBROUTINE WriteStringTwoValueLL

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted logical firstValue and the the secondString followed by a formatted single precision secondValue to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WriteStringTwoValueLSP(id,firstString,firstValue,secondString,secondValue,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    LOGICAL, INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    REAL(SP), INTENT(IN) :: secondValue !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

 !   ENTERS("WriteStringTwoValueLSP",err,error,*999)
        
    localString=firstString//LogicalToVString(firstValue,err,error)
    IF(err/=0) GOTO 999
    localString=localString//secondString//NumberToVString(secondValue,"*",err,error)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringTwoValueLSP")
    RETURN
999 ERRORS("WriteStringTwoValueLSP",err,error)
!    EXITS("WriteStringTwoValueLSP")
    RETURN 1
    
  END SUBROUTINE WriteStringTwoValueLSP

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted logical firstValue and the the secondString followed by a formatted single precision secondValue to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WriteStringTwoValueLVS(id,firstString,firstValue,secondString,secondValue,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    LOGICAL, INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: secondValue !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringTwoValueLVS",err,error,*999)
        
    localString=firstString//LogicalToVString(firstValue,err,error)//secondString//secondValue
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringTwoValueLVS")
    RETURN
999 ERRORS("WriteStringTwoValueLVS",err,error)
!    EXITS("WriteStringTwoValueLVS")
    RETURN 1
    
  END SUBROUTINE WriteStringTwoValueLVS

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted single precision firstValue and the the secondString followed by a formatted character secondValue to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WriteStringTwoValueSPC(id,firstString,firstValue,secondString,secondValue,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    REAL(SP), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondValue !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringTwoValueSPC",err,error,*999)
        
    localString=firstString//NumberToVString(firstValue,"*",err,error)//secondString//secondValue
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringTwoValueSPC")
    RETURN
999 ERRORS("WriteStringTwoValueSPC",err,error)
!    EXITS("WriteStringTwoValueSPC")
    RETURN 1
    
  END SUBROUTINE WriteStringTwoValueSPC

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted single precision firstValue and the the secondString followed by a formatted double precision secondValue to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WriteStringTwoValueSPDP(id,firstString,firstValue,secondString,secondValue,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    REAL(SP), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    REAL(DP), INTENT(IN) :: secondValue !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString,localString2

 !   ENTERS("WriteStringTwoValueSPDP",err,error,*999)
        
    localString=firstString//NumberToVString(firstValue,"*",err,error)
    IF(err/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    localString2=localString//secondString
    localString=localString2//NumberToVString(secondValue,"*",err,error)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringTwoValueSPDP")
    RETURN
999 ERRORS("WriteStringTwoValueSPDP",err,error)
!    EXITS("WriteStringTwoValueSPDP")
    RETURN 1
    
  END SUBROUTINE WriteStringTwoValueSPDP

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted single precision firstValue and the the secondString followed by a formatted integer secondValue to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WriteStringTwoValueSPIntg(id,firstString,firstValue,secondString,secondValue,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    REAL(SP), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    INTEGER(INTG), INTENT(IN) :: secondValue !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString,localString2

!    ENTERS("WriteStringTwoValueSPIntg",err,error,*999)
        
    localString=firstString//NumberToVString(firstValue,"*",err,error)
    IF(err/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    localString2=localString//secondString
    localString=localString2//NumberToVString(secondValue,"*",err,error)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringTwoValueSPIntg")
    RETURN
999 ERRORS("WriteStringTwoValueSPIntg",err,error)
!    EXITS("WriteStringTwoValueSPIntg")
    RETURN 1
    
  END SUBROUTINE WriteStringTwoValueSPIntg

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted single precision firstValue and the the secondString followed by a formatted logical secondValue to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WriteStringTwoValueSPL(id,firstString,firstValue,secondString,secondValue,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    REAL(SP), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    LOGICAL, INTENT(IN) :: secondValue !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringTwoValueSPL",err,error,*999)
        
    localString=firstString//NumberToVString(firstValue,"*",err,error)
    IF(err/=0) GOTO 999
    localString=localString//secondString//LogicalToVString(secondValue,err,error)
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringTwoValueSPL")
    RETURN
999 ERRORS("WriteStringTwoValueSPL",err,error)
!    EXITS("WriteStringTwoValueSPL")
    RETURN 1
    
  END SUBROUTINE WriteStringTwoValueSPL

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted single precision firstValue and the the secondString followed by a formatted single precision secondValue to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WriteStringTwoValueSPSP(id,firstString,firstValue,secondString,secondValue,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    REAL(SP), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    REAL(SP), INTENT(IN) :: secondValue !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString,localString2

 !   ENTERS("WriteStringTwoValueSPSP",err,error,*999)
        
    localString=firstString//NumberToVString(firstValue,"*",err,error)
    IF(err/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    localString2=localString//secondString
    localString=localString2//NumberToVString(secondValue,"*",err,error)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringTwoValueSPSP")
    RETURN
999 ERRORS("WriteStringTwoValueSPSP",err,error)
!    EXITS("WriteStringTwoValueSPSP")
    RETURN 1
    
  END SUBROUTINE WriteStringTwoValueSPSP

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted single precision firstValue and the the secondString followed by a formatted single precision secondValue to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WriteStringTwoValueSPVS(id,firstString,firstValue,secondString,secondValue,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    REAL(SP), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: secondValue !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringTwoValueSPVS",err,error,*999)
        
    localString=firstString//NumberToVString(firstValue,"*",err,error)//secondString//secondValue
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringTwoValueSPVS")
    RETURN
999 ERRORS("WriteStringTwoValueSPVS",err,error)
!    EXITS("WriteStringTwoValueSPVS")
    RETURN 1
    
  END SUBROUTINE WriteStringTwoValueSPVS

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted varying string firstValue and the the secondString followed by a formatted character secondValue to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WriteStringTwoValueVSC(id,firstString,firstValue,secondString,secondValue,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondValue !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringTwoValueVSC",err,error,*999)
        
    localString=firstString//firstValue//secondString//secondValue
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringTwoValueVSC")
    RETURN
999 ERRORS("WriteStringTwoValueVSC",err,error)
!    EXITS("WriteStringTwoValueVSC")
    RETURN 1
    
  END SUBROUTINE WriteStringTwoValueVSC

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted varying string firstValue and the the secondString followed by a formatted double precision secondValue to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WriteStringTwoValueVSDP(id,firstString,firstValue,secondString,secondValue,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    REAL(DP), INTENT(IN) :: secondValue !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

 !   ENTERS("WriteStringTwoValueVSDP",err,error,*999)
        
    localString=firstString//firstValue//secondString//NumberToVString(secondValue,"*",err,error)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringTwoValueVSDP")
    RETURN
999 ERRORS("WriteStringTwoValueVSDP",err,error)
!    EXITS("WriteStringTwoValueVSDP")
    RETURN 1
    
  END SUBROUTINE WriteStringTwoValueVSDP

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted varying string firstValue and the the secondString followed by a formatted integer secondValue to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WriteStringTwoValueVSIntg(id,firstString,firstValue,secondString,secondValue,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    INTEGER(INTG), INTENT(IN) :: secondValue !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringTwoValueVSIntg",err,error,*999)
        
    localString=firstString//firstValue//secondString//NumberToVString(secondValue,"*",err,error)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringTwoValueVSIntg")
    RETURN
999 ERRORS("WriteStringTwoValueVSIntg",err,error)
!    EXITS("WriteStringTwoValueVSIntg")
    RETURN 1
    
  END SUBROUTINE WriteStringTwoValueVSIntg

  !
  !================================================================================================================================
  !
  
  !>Writes the firstString followed by a formatted varying string firstValue and the the secondString followed by a formatted logical secondValue to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WriteStringTwoValueVSL(id,firstString,firstValue,secondString,secondValue,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    LOGICAL, INTENT(IN) :: secondValue !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringTwoValueVSL",err,error,*999)
        
    localString=firstString//firstValue//secondString//LogicalToVString(secondValue,err,error)
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringTwoValueVSL")
    RETURN
999 ERRORS("WriteStringTwoValueVSL",err,error)
!    EXITS("WriteStringTwoValueVSL")
    RETURN 1
    
  END SUBROUTINE WriteStringTwoValueVSL

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted varying string firstValue and the the secondString followed by a formatted single precision secondValue to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WriteStringTwoValueVSSP(id,firstString,firstValue,secondString,secondValue,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    REAL(SP), INTENT(IN) :: secondValue !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

 !   ENTERS("WriteStringTwoValueVSSP",err,error,*999)
        
    localString=firstString//firstValue//secondString//NumberToVString(secondValue,"*",err,error)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringTwoValueVSSP")
    RETURN
999 ERRORS("WriteStringTwoValueVSSP",err,error)
!    EXITS("WriteStringTwoValueVSSP")
    RETURN 1
    
  END SUBROUTINE WriteStringTwoValueVSSP

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted varying string firstValue and the the secondString followed by a formatted varying string secondValue to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WriteStringTwoValueVSVS(id,firstString,firstValue,secondString,secondValue,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: secondValue !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringTwoValueVSVS",err,error,*999)
        
    localString=firstString//firstValue//secondString//secondValue
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringTwoValueVSVS")
    RETURN
999 ERRORS("WriteStringTwoValueVSVS",err,error)
!    EXITS("WriteStringTwoValueVSVS")
    RETURN 1
    
  END SUBROUTINE WriteStringTwoValueVSVS

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted character value to the given output stream specified by ID. formatString is used to format the value.
  SUBROUTINE WriteStringFmtValueC(id,firstString,value,formatString,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    CHARACTER(LEN=*), INTENT(IN) :: value !<The value to be output
    CHARACTER(LEN=*), INTENT(IN) :: formatString !<The format string to be used to format the value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringFmtValueC",err,error,*999)
        
    localString=firstString(1:LEN_TRIM(firstString))//value(1:LEN_TRIM(value))
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtValueC")
    RETURN
999 ERRORS("WriteStringFmtValueC",err,error)
!    EXITS("WriteStringFmtValueC")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtValueC

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted character value to the given output stream specified by ID. formatString is used to format the value.
  SUBROUTINE WriteStringFmtValueDP(id,firstString,value,formatString,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    REAL(DP), INTENT(IN) :: value !<The value to be output
    CHARACTER(LEN=*), INTENT(IN) :: formatString !<The format string to be used to format the value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringFmtValueDP",err,error,*999)
        
    localString=firstString//NumberToVString(value,formatString,err,error,adjust=.FALSE.)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtValueDP")
    RETURN
999 ERRORS("WriteStringFmtValueDP",err,error)
!    EXITS("WriteStringFmtValueDP")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtValueDP
  
  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted character value to the given output stream specified by ID. formatString is used to format the value.
  SUBROUTINE WriteStringFmtValueIntg(id,firstString,value,formatString,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    INTEGER(INTG), INTENT(IN) :: value !<The value to be output
    CHARACTER(LEN=*), INTENT(IN) :: formatString !<The format string to be used to format the value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringFmtValueIntg",err,error,*999)
        
    localString=firstString//NumberToVString(value,formatString,err,error,adjust=.FALSE.)
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtValueIntg")
    RETURN
999 ERRORS("WriteStringFmtValueIntg",err,error)
!    EXITS("WriteStringFmtValueIntg")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtValueIntg

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted character value to the given output stream specified by ID. formatString is used to format the value.
  SUBROUTINE WriteStringFmtValueLIntg(id,firstString,value,formatString,err,error,*)
7
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    INTEGER(LINTG), INTENT(IN) :: value !<The value to be output
    CHARACTER(LEN=*), INTENT(IN) :: formatString !<The format string to be used to format the value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringFmtValueLIntg",err,error,*999)
        
    localString=firstString//NumberToVString(value,formatString,err,error,adjust=.FALSE.)
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtValueLIntg")
    RETURN
999 ERRORS("WriteStringFmtValueLIntg",err,error)
!    EXITS("WriteStringFmtValueLIntg")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtValueLIntg

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted character value to the given output stream specified by ID. formatString is used to format the value.
  SUBROUTINE WriteStringFmtValueL(id,firstString,value,formatString,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    LOGICAL, INTENT(IN) :: value !<The value to be output
    CHARACTER(LEN=*), INTENT(IN) :: formatString !<The format string to be used to format the value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringFmtValueL",err,error,*999)
        
    localString=firstString//LogicalToVString(value,err,error)
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtValueL")
    RETURN
999 ERRORS("WriteStringFmtValueL",err,error)
!    EXITS("WriteStringFmtValueL")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtValueL

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted character value to the given output stream specified by ID. formatString is used to format the value.
  SUBROUTINE WriteStringFmtValueSP(id,firstString,value,formatString,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    REAL(SP), INTENT(IN) :: value !<The value to be output
    CHARACTER(LEN=*), INTENT(IN) :: formatString !<The format string to be used to format the value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringFmtValueSP",err,error,*999)
        
    localString=firstString//NumberToVString(value,formatString,err,error,adjust=.FALSE.)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtValueSP")
    RETURN
999 ERRORS("WriteStringFmtValueSP",err,error)
!    EXITS("WriteStringFmtValueSP")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtValueSP

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted character value to the given output stream specified by ID. formatString is used to format the value.
  SUBROUTINE WriteStringFmtValueVS(id,firstString,value,formatString,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: value !<The value to be output
    CHARACTER(LEN=*), INTENT(IN) :: formatString !<The format string to be used to format the value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringFmtValueVS",err,error,*999)
        
    localString=firstString//value
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtValueVS")
    RETURN
999 ERRORS("WriteStringFmtValueVS",err,error)
!    EXITS("WriteStringFmtValueVS")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtValueVS

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted character firstValue and the the secondString followed by a formatted character secondValue to the given output stream specified by ID. firstFormat is used to format the first value and secondFormat is used to format the second value.
  SUBROUTINE WriteStringFmtTwoValueCC(id,firstString,firstValue,firstFormat,secondString,secondValue,secondFormat,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondValue !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondFormat !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringFmtTwoValueCC",err,error,*999)
        
    localString=firstString//firstValue//secondString//secondValue
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtTwoValueCC")
    RETURN
999 ERRORS("WriteStringFmtTwoValueCC",err,error)
!    EXITS("WriteStringFmtTwoValueCC")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtTwoValueCC

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted character firstValue and the the secondString followed by a formatted double precision secondValue to the given output stream specified by ID. firstFormat is used to format the first value and secondFormat is used to format the second value.
  SUBROUTINE WriteStringFmtTwoValueCDP(id,firstString,firstValue,firstFormat,secondString,secondValue,secondFormat,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    REAL(DP), INTENT(IN) :: secondValue !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondFormat !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

 !   ENTERS("WriteStringFmtTwoValueCDP",err,error,*999)
        
    localString=firstString//firstValue//secondString//NumberToVString(secondValue,secondFormat,err,error,adjust=.FALSE.)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtTwoValueCDP")
    RETURN
999 ERRORS("WriteStringFmtTwoValueCDP",err,error)
!    EXITS("WriteStringFmtTwoValueCDP")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtTwoValueCDP

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted character firstValue and the the secondString followed by a formatted integer secondValue to the given output stream specified by ID. firstFormat is used to format the first value and secondFormat is used to format the second value.
  SUBROUTINE WriteStringFmtTwoValueCIntg(id,firstString,firstValue,firstFormat,secondString,secondValue,secondFormat,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    INTEGER(INTG), INTENT(IN) :: secondValue !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondFormat !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringFmtTwoValueCIntg",err,error,*999)
        
    localString=firstString//firstValue//secondString//NumberToVString(secondValue,secondFormat,err,error,adjust=.FALSE.)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtTwoValueCIntg")
    RETURN
999 ERRORS("WriteStringFmtTwoValueCIntg",err,error)
!    EXITS("WriteStringFmtTwoValueCIntg")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtTwoValueCIntg

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted character firstValue and the the secondString followed by a formatted logical secondValue to the given output stream specified by ID. firstFormat is used to format the first value and secondFormat is used to format the second value.
  SUBROUTINE WriteStringFmtTwoValueCL(id,firstString,firstValue,firstFormat,secondString,secondValue,secondFormat,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    LOGICAL, INTENT(IN) :: secondValue !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondFormat !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringFmtTwoValueCL",err,error,*999)
        
    localString=firstString//firstValue//secondString//LogicalToVString(secondValue,err,error)
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtTwoValueCL")
    RETURN
999 ERRORS("WriteStringFmtTwoValueCL",err,error)
!    EXITS("WriteStringFmtTwoValueCL")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtTwoValueCL

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted character firstValue and the the secondString followed by a formatted single precision secondValue to the given output stream specified by ID. firstFormat is used to format the first value and secondFormat is used to format the second value.
  SUBROUTINE WriteStringFmtTwoValueCSP(id,firstString,firstValue,firstFormat,secondString,secondValue,secondFormat,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    REAL(SP), INTENT(IN) :: secondValue !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondFormat !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

 !   ENTERS("WriteStringFmtTwoValueCSP",err,error,*999)
        
    localString=firstString//firstValue//secondString//NumberToVString(secondValue,secondFormat,err,error,adjust=.FALSE.)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtTwoValueCSP")
    RETURN
999 ERRORS("WriteStringFmtTwoValueCSP",err,error)
!    EXITS("WriteStringFmtTwoValueCSP")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtTwoValueCSP

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted character firstValue and the the secondString followed by a formatted varying string secondValue to the given output stream specified by ID. firstFormat is used to format the first value and secondFormat is used to format the second value.
  SUBROUTINE WriteStringFmtTwoValueCVS(id,firstString,firstValue,firstFormat,secondString,secondValue,secondFormat,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: secondValue !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondFormat !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringFmtTwoValueCVS",err,error,*999)
        
    localString=firstString//firstValue//secondString//secondValue
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtTwoValueCVS")
    RETURN
999 ERRORS("WriteStringFmtTwoValueCVS",err,error)
!    EXITS("WriteStringFmtTwoValueCVS")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtTwoValueCVS

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted double precision firstValue and the the secondString followed by a formatted character secondValue to the given output stream specified by ID. firstFormat is used to format the first value and secondFormat is used to format the second value.
  SUBROUTINE WriteStringFmtTwoValueDPC(id,firstString,firstValue,firstFormat,secondString,secondValue,secondFormat,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    REAL(DP), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondValue !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondFormat !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringFmtTwoValueDPC",err,error,*999)
        
    localString=firstString//NumberToVString(firstValue,firstFormat,err,error,adjust=.FALSE.)//secondString//secondValue
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtTwoValueDPC")
    RETURN
999 ERRORS("WriteStringFmtTwoValueDPC",err,error)
!    EXITS("WriteStringFmtTwoValueDPC")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtTwoValueDPC

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted double precision firstValue and the the secondString followed by a formatted double precision secondValue to the given output stream specified by ID. firstFormat is used to format the first value and secondFormat is used to format the second value.
  SUBROUTINE WriteStringFmtTwoValueDPDP(id,firstString,firstValue,firstFormat,secondString,secondValue,secondFormat,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    REAL(DP), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    REAL(DP), INTENT(IN) :: secondValue !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondFormat !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString,localString2

 !   ENTERS("WriteStringFmtTwoValueDPDP",err,error,*999)
        
    localString=firstString//NumberToVString(firstValue,firstFormat,err,error,adjust=.FALSE.)
    IF(err/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    localString2=localString//secondString
    localString=localString2//NumberToVString(secondValue,secondFormat,err,error,adjust=.FALSE.)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtTwoValueDPDP")
    RETURN
999 ERRORS("WriteStringFmtTwoValueDPDP",err,error)
!    EXITS("WriteStringFmtTwoValueDPDP")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtTwoValueDPDP

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted double precision firstValue and the the secondString followed by a formatted integer secondValue to the given output stream specified by ID. firstFormat is used to format the first value and secondFormat is used to format the second value.
  SUBROUTINE WriteStringFmtTwoValueDPIntg(id,firstString,firstValue,firstFormat,secondString,secondValue,secondFormat,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    REAL(DP), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    INTEGER(INTG), INTENT(IN) :: secondValue !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondFormat !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString,localString2

!    ENTERS("WriteStringFmtTwoValueDPIntg",err,error,*999)
        
    localString=firstString//NumberToVString(firstValue,firstFormat,err,error,adjust=.FALSE.)
    IF(err/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    localString2=localString//secondString
    localString=localString2//NumberToVString(secondValue,secondFormat,err,error,adjust=.FALSE.)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtTwoValueDPIntg")
    RETURN
999 ERRORS("WriteStringFmtTwoValueDPIntg",err,error)
!    EXITS("WriteStringFmtTwoValueDPIntg")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtTwoValueDPIntg

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted double precision firstValue and the the secondString followed by a formatted logical secondValue to the given output stream specified by ID. firstFormat is used to format the first value and secondFormat is used to format the second value.
  SUBROUTINE WriteStringFmtTwoValueDPL(id,firstString,firstValue,firstFormat,secondString,secondValue,secondFormat,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    REAL(DP), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    LOGICAL, INTENT(IN) :: secondValue !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondFormat !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString,localString2

!    ENTERS("WriteStringFmtTwoValueDPL",err,error,*999)
        
    localString=firstString//NumberToVString(firstValue,firstFormat,err,error,adjust=.FALSE.)
    IF(err/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    localString2=localString//secondString
    localString=localString2//LogicalToVString(secondValue,err,error)
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtTwoValueDPL")
    RETURN
999 ERRORS("WriteStringFmtTwoValueDPL",err,error)
!    EXITS("WriteStringFmtTwoValueDPL")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtTwoValueDPL

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted double precision firstValue and the the secondString followed by a formatted single precision secondValue to the given output stream specified by ID. firstFormat is used to format the first value and secondFormat is used to format the second value.
  SUBROUTINE WriteStringFmtTwoValueDPSP(id,firstString,firstValue,firstFormat,secondString,secondValue,secondFormat,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    REAL(DP), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    REAL(SP), INTENT(IN) :: secondValue !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondFormat !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString,localString2

 !   ENTERS("WriteStringFmtTwoValueDPSP",err,error,*999)
        
    localString=firstString//NumberToVString(firstValue,firstFormat,err,error,adjust=.FALSE.)
    IF(err/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    localString2=localString//secondString
    localString=localString2//NumberToVString(secondValue,secondFormat,err,error,adjust=.FALSE.)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtTwoValueDPSP")
    RETURN
999 ERRORS("WriteStringFmtTwoValueDPSP",err,error)
!    EXITS("WriteStringFmtTwoValueDPSP")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtTwoValueDPSP

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted double precision firstValue and the the secondString followed by a formatted single precision secondValue to the given output stream specified by ID. firstFormat is used to format the first value and secondFormat is used to format the second value.
  SUBROUTINE WriteStringFmtTwoValueDPVS(id,firstString,firstValue,firstFormat,secondString,secondValue,secondFormat,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    REAL(DP), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: secondValue !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondFormat !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringFmtTwoValueDPVS",err,error,*999)
        
    localString=firstString//NumberToVString(firstValue,firstFormat,err,error,adjust=.FALSE.)//secondString//secondValue
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtTwoValueDPVS")
    RETURN
999 ERRORS("WriteStringFmtTwoValueDPVS",err,error)
!    EXITS("WriteStringFmtTwoValueDPVS")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtTwoValueDPVS

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted integer firstValue and the the secondString followed by a formatted character secondValue to the given output stream specified by ID. firstFormat is used to format the first value and secondFormat is used to format the second value.
  SUBROUTINE WriteStringFmtTwoValueIntgC(id,firstString,firstValue,firstFormat,secondString,secondValue,secondFormat,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    INTEGER(INTG), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondValue !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondFormat !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringFmtTwoValueIntgC",err,error,*999)
        
    localString=firstString//NumberToVString(firstValue,firstFormat,err,error,adjust=.FALSE.)//secondString//secondValue
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtTwoValueIntgC")
    RETURN
999 ERRORS("WriteStringFmtTwoValueIntgC",err,error)
!    EXITS("WriteStringFmtTwoValueIntgC")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtTwoValueIntgC

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted integer firstValue and the the secondString followed by a formatted double precision secondValue to the given output stream specified by ID. firstFormat is used to format the first value and secondFormat is used to format the second value.
  SUBROUTINE WriteStringFmtTwoValueIntgDP(id,firstString,firstValue,firstFormat,secondString,secondValue,secondFormat,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    INTEGER(INTG), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    REAL(DP), INTENT(IN) :: secondValue !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondFormat !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString,localString2

 !   ENTERS("WriteStringFmtTwoValueIntgDP",err,error,*999)
        
    localString=firstString//NumberToVString(firstValue,firstFormat,err,error,adjust=.FALSE.)
    IF(err/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    localString2=localString//secondString
    localString=localString2//NumberToVString(secondValue,secondFormat,err,error,adjust=.FALSE.)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtTwoValueIntgDP")
    RETURN
999 ERRORS("WriteStringFmtTwoValueIntgDP",err,error)
!    EXITS("WriteStringFmtTwoValueIntgDP")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtTwoValueIntgDP

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted integer firstValue and the the secondString followed by a formatted integer secondValue to the given output stream specified by ID. firstFormat is used to format the first value and secondFormat is used to format the second value.
  SUBROUTINE WriteStringFmtTwoValueIntgIntg(id,firstString,firstValue,firstFormat,secondString,secondValue,secondFormat,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    INTEGER(INTG), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    INTEGER(INTG), INTENT(IN) :: secondValue !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondFormat !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString,localString2

!    ENTERS("WriteStringFmtTwoValueIntgIntg",err,error,*999)
        
    localString=firstString//NumberToVString(firstValue,firstFormat,err,error,adjust=.FALSE.)
    IF(err/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    localString2=localString//secondString
    localString=localString2//NumberToVString(secondValue,secondFormat,err,error,adjust=.FALSE.)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtTwoValueIntgIntg")
    RETURN
999 ERRORS("WriteStringFmtTwoValueIntgIntg",err,error)
!    EXITS("WriteStringFmtTwoValueIntgIntg")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtTwoValueIntgIntg

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted integer firstValue and the the secondString followed by a formatted logical secondValue to the given output stream specified by ID. firstFormat is used to format the first value and secondFormat is used to format the second value.
  SUBROUTINE WriteStringFmtTwoValueIntgL(id,firstString,firstValue,firstFormat,secondString,secondValue,secondFormat,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    INTEGER(INTG), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    LOGICAL, INTENT(IN) :: secondValue !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondFormat !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString,localString2

!    ENTERS("WriteStringFmtTwoValueIntgL",err,error,*999)
        
    localString=firstString//NumberToVString(firstValue,firstFormat,err,error,adjust=.FALSE.)
    IF(err/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    localString2=localString//secondString
    localString=localString2//LogicalToVString(secondValue,err,error)
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtTwoValueIntgL")
    RETURN
999 ERRORS("WriteStringFmtTwoValueIntgL",err,error)
!    EXITS("WriteStringFmtTwoValueIntgL")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtTwoValueIntgL

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted integer firstValue and the the secondString followed by a formatted single precision secondValue to the given output stream specified by ID. firstFormat is used to format the first value and secondFormat is used to format the second value.
  SUBROUTINE WriteStringFmtTwoValueIntgSP(id,firstString,firstValue,firstFormat,secondString,secondValue,secondFormat,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    INTEGER(INTG), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    REAL(SP), INTENT(IN) :: secondValue !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondFormat !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString,localString2

 !   ENTERS("WriteStringFmtTwoValueIntgSP",err,error,*999)
        
    localString=firstString//NumberToVString(firstValue,firstFormat,err,error,adjust=.FALSE.)
    IF(err/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    localString2=localString//secondString
    localString=localString2//NumberToVString(secondValue,secondFormat,err,error,adjust=.FALSE.)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtTwoValueIntgSP")
    RETURN
999 ERRORS("WriteStringFmtTwoValueIntgSP",err,error)
!    EXITS("WriteStringFmtTwoValueIntgSP")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtTwoValueIntgSP

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted integer firstValue and the the secondString followed by a formatted single precision secondValue to the given output stream specified by ID. firstFormat is used to format the first value and secondFormat is used to format the second value.
  SUBROUTINE WriteStringFmtTwoValueIntgVS(id,firstString,firstValue,firstFormat,secondString,secondValue,secondFormat,err,error,*)
   
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    INTEGER(INTG), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: secondValue !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondFormat !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringFmtTwoValueIntgVS",err,error,*999)
        
    localString=firstString//NumberToVString(firstValue,firstFormat,err,error,adjust=.FALSE.)//secondString//secondValue
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtTwoValueIntgVS")
    RETURN
999 ERRORS("WriteStringFmtTwoValueIntgVS",err,error)
!    EXITS("WriteStringFmtTwoValueIntgVS")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtTwoValueIntgVS

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted logical firstValue and the the secondString followed by a formatted character secondValue to the given output stream specified by ID. firstFormat is used to format the first value and secondFormat is used to format the second value.
  SUBROUTINE WriteStringFmtTwoValueLC(id,firstString,firstValue,firstFormat,secondString,secondValue,secondFormat,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    LOGICAL, INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondValue !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondFormat !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringFmtTwoValueLC",err,error,*999)
        
    localString=firstString//LogicalToVString(firstValue,err,error)//secondString//secondValue
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtTwoValueLC")
    RETURN
999 ERRORS("WriteStringFmtTwoValueLC",err,error)
!    EXITS("WriteStringFmtTwoValueLC")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtTwoValueLC

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted logical firstValue and the the secondString followed by a formatted double precision secondValue to the given output stream specified by ID. firstFormat is used to format the first value and secondFormat is used to format the second value.
  SUBROUTINE WriteStringFmtTwoValueLDP(id,firstString,firstValue,firstFormat,secondString,secondValue,secondFormat,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    LOGICAL, INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    REAL(DP), INTENT(IN) :: secondValue !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondFormat !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString,localString2

 !   ENTERS("WriteStringFmtTwoValueLDP",err,error,*999)
        
    localString=firstString//LogicalToVString(firstValue,err,error)
    IF(err/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    localString2=localString//secondString
    localString=localString2//NumberToVString(secondValue,secondFormat,err,error,adjust=.FALSE.)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtTwoValueLDP")
    RETURN
999 ERRORS("WriteStringFmtTwoValueLDP",err,error)
!    EXITS("WriteStringFmtTwoValueLDP")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtTwoValueLDP

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted logical firstValue and the the secondString followed by a formatted integer secondValue to the given output stream specified by ID. firstFormat is used to format the first value and secondFormat is used to format the second value.
  SUBROUTINE WriteStringFmtTwoValueLIntg(id,firstString,firstValue,firstFormat,secondString,secondValue,secondFormat,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    LOGICAL, INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    INTEGER(INTG), INTENT(IN) :: secondValue !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondFormat !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringFmtTwoValueLIntg",err,error,*999)
        
    localString=firstString//LogicalToVString(firstValue,err,error)
    IF(err/=0) GOTO 999
    localString=localString//secondString//NumberToVString(secondValue,secondFormat,err,error,adjust=.FALSE.)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtTwoValueLIntg")
    RETURN
999 ERRORS("WriteStringFmtTwoValueLIntg",err,error)
!    EXITS("WriteStringFmtTwoValueLIntg")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtTwoValueLIntg

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted logical firstValue and the the secondString followed by a formatted logical secondValue to the given output stream specified by ID. firstFormat is used to format the first value and secondFormat is used to format the second value.
  SUBROUTINE WriteStringFmtTwoValueLL(id,firstString,firstValue,firstFormat,secondString,secondValue,secondFormat,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    LOGICAL, INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    LOGICAL, INTENT(IN) :: secondValue !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondFormat !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString,localString2

!    ENTERS("WriteStringFmtTwoValueLL",err,error,*999)
        
    localString=firstString//LogicalToVString(firstValue,err,error)
    IF(err/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    localString2=localString//secondString
    localString=localString2//LogicalToVString(secondValue,err,error)
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WRITE_STRING_TWO_L_INTG_L")
    RETURN
999 ERRORS("WriteStringFmtTwoValueLL",err,error)
!    EXITS("WriteStringFmtTwoValueLL")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtTwoValueLL

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted logical firstValue and the the secondString followed by a formatted single precision secondValue to the given output stream specified by ID. firstFormat is used to format the first value and secondFormat is used to format the second value.
  SUBROUTINE WriteStringFmtTwoValueLSP(id,firstString,firstValue,firstFormat,secondString,secondValue,secondFormat,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    LOGICAL, INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    REAL(SP), INTENT(IN) :: secondValue !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondFormat !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

 !   ENTERS("WriteStringFmtTwoValueLSP",err,error,*999)
        
    localString=firstString//LogicalToVString(firstValue,err,error)
    IF(err/=0) GOTO 999
    localString=localString//secondString//NumberToVString(secondValue,secondFormat,err,error,adjust=.FALSE.)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtTwoValueLSP")
    RETURN
999 ERRORS("WriteStringFmtTwoValueLSP",err,error)
!    EXITS("WriteStringFmtTwoValueLSP")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtTwoValueLSP

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted logical firstValue and the the secondString followed by a formatted single precision secondValue to the given output stream specified by ID. firstFormat is used to format the first value and secondFormat is used to format the second value.
  SUBROUTINE WriteStringFmtTwoValueLVS(id,firstString,firstValue,firstFormat,secondString,secondValue,secondFormat,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    LOGICAL, INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: secondValue !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondFormat !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringFmtTwoValueLVS",err,error,*999)
        
    localString=firstString//LogicalToVString(firstValue,err,error)//secondString//secondValue
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtTwoValueLVS")
    RETURN
999 ERRORS("WriteStringFmtTwoValueLVS",err,error)
!    EXITS("WriteStringFmtTwoValueLVS")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtTwoValueLVS

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted single precision firstValue and the the secondString followed by a formatted character secondValue to the given output stream specified by ID. firstFormat is used to format the first value and secondFormat is used to format the second value.
  SUBROUTINE WriteStringFmtTwoValueSPC(id,firstString,firstValue,firstFormat,secondString,secondValue,secondFormat,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    REAL(SP), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondValue !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondFormat !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringFmtTwoValueSPC",err,error,*999)
        
    localString=firstString//NumberToVString(firstValue,firstFormat,err,error,adjust=.FALSE.)//secondString//secondValue
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtTwoValueSPC")
    RETURN
999 ERRORS("WriteStringFmtTwoValueSPC",err,error)
!    EXITS("WriteStringFmtTwoValueSPC")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtTwoValueSPC

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted single precision firstValue and the the secondString followed by a formatted double precision secondValue to the given output stream specified by ID. firstFormat is used to format the first value and secondFormat is used to format the second value.
  SUBROUTINE WriteStringFmtTwoValueSPDP(id,firstString,firstValue,firstFormat,secondString,secondValue,secondFormat,err,error,*)
   
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    REAL(SP), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    REAL(DP), INTENT(IN) :: secondValue !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondFormat !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString,localString2

 !   ENTERS("WriteStringFmtTwoValueSPDP",err,error,*999)
        
    localString=firstString//NumberToVString(firstValue,firstFormat,err,error,adjust=.FALSE.)
    IF(err/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    localString2=localString//secondString
    localString=localString2//NumberToVString(secondValue,secondFormat,err,error,adjust=.FALSE.)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtTwoValueSPDP")
    RETURN
999 ERRORS("WriteStringFmtTwoValueSPDP",err,error)
!    EXITS("WriteStringFmtTwoValueSPDP")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtTwoValueSPDP

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted single precision firstValue and the the secondString followed by a formatted integer secondValue to the given output stream specified by ID. firstFormat is used to format the first value and secondFormat is used to format the second value.
  SUBROUTINE WriteStringFmtTwoValueSPIntg(id,firstString,firstValue,firstFormat,secondString,secondValue,secondFormat,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    REAL(SP), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    INTEGER(INTG), INTENT(IN) :: secondValue !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondFormat !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString,localString2

!    ENTERS("WriteStringFmtTwoValueSPIntg",err,error,*999)
        
    localString=firstString//NumberToVString(firstValue,firstFormat,err,error,adjust=.FALSE.)
    IF(err/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    localString2=localString//secondString
    localString=localString2//NumberToVString(secondValue,secondFormat,err,error,adjust=.FALSE.)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtTwoValueSPIntg")
    RETURN
999 ERRORS("WriteStringFmtTwoValueSPIntg",err,error)
!    EXITS("WriteStringFmtTwoValueSPIntg")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtTwoValueSPIntg

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted single precision firstValue and the the secondString followed by a formatted logical secondValue to the given output stream specified by ID. firstFormat is used to format the first value and secondFormat is used to format the second value.
  SUBROUTINE WriteStringFmtTwoValueSPL(id,firstString,firstValue,firstFormat,secondString,secondValue,secondFormat,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    REAL(SP), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    LOGICAL, INTENT(IN) :: secondValue !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondFormat !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringFmtTwoValueSPL",err,error,*999)
        
    localString=firstString//NumberToVString(firstValue,firstFormat,err,error,adjust=.FALSE.)
    IF(err/=0) GOTO 999
    localString=localString//secondString//LogicalToVString(secondValue,err,error)
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtTwoValueSPL")
    RETURN
999 ERRORS("WriteStringFmtTwoValueSPL",err,error)
!    EXITS("WriteStringFmtTwoValueSPL")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtTwoValueSPL

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted single precision firstValue and the the secondString followed by a formatted single precision secondValue to the given output stream specified by ID. firstFormat is used to format the first value and secondFormat is used to format the second value.
  SUBROUTINE WriteStringFmtTwoValueSPSP(id,firstString,firstValue,firstFormat,secondString,secondValue,secondFormat,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    REAL(SP), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    REAL(SP), INTENT(IN) :: secondValue !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondFormat !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString,localString2

 !   ENTERS("WriteStringFmtTwoValueSPSP",err,error,*999)
        
    localString=firstString//NumberToVString(firstValue,firstFormat,err,error,adjust=.FALSE.)
    IF(err/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    localString2=localString//secondString
    localString=localString2//NumberToVString(secondValue,secondFormat,err,error,adjust=.FALSE.)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtTwoValueSPSP")
    RETURN
999 ERRORS("WriteStringFmtTwoValueSPSP",err,error)
!    EXITS("WriteStringFmtTwoValueSPSP")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtTwoValueSPSP

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted single precision firstValue and the the secondString followed by a formatted single precision secondValue to the given output stream specified by ID. firstFormat is used to format the first value and secondFormat is used to format the second value.
  SUBROUTINE WriteStringFmtTwoValueSPVS(id,firstString,firstValue,firstFormat,secondString,secondValue,secondFormat,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    REAL(SP), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: secondValue !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondFormat !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringFmtTwoValueSPVS",err,error,*999)
        
    localString=firstString//NumberToVString(firstValue,firstFormat,err,error,adjust=.FALSE.)//secondString//secondValue
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtTwoValueSPVS")
    RETURN
999 ERRORS("WriteStringFmtTwoValueSPVS",err,error)
!    EXITS("WriteStringFmtTwoValueSPVS")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtTwoValueSPVS

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted varying string firstValue and the the secondString followed by a formatted character secondValue to the given output stream specified by ID. firstFormat is used to format the first value and secondFormat is used to format the second value.
  SUBROUTINE WriteStringFmtTwoValueVSC(id,firstString,firstValue,firstFormat,secondString,secondValue,secondFormat,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondValue !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondFormat !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringFmtTwoValueVSC",err,error,*999)
        
    localString=firstString//firstValue//secondString//secondValue
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtTwoValueVSC")
    RETURN
999 ERRORS("WriteStringFmtTwoValueVSC",err,error)
!    EXITS("WriteStringFmtTwoValueVSC")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtTwoValueVSC

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted varying string firstValue and the the secondString followed by a formatted double precision secondValue to the given output stream specified by ID. firstFormat is used to format the first value and secondFormat is used to format the second value.
  SUBROUTINE WriteStringFmtTwoValueVSDP(id,firstString,firstValue,firstFormat,secondString,secondValue,secondFormat,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    REAL(DP), INTENT(IN) :: secondValue !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondFormat !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

 !   ENTERS("WriteStringFmtTwoValueVSDP",err,error,*999)
        
    localString=firstString//firstValue//secondString//NumberToVString(secondValue,secondFormat,err,error,adjust=.FALSE.)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtTwoValueVSDP")
    RETURN
999 ERRORS("WriteStringFmtTwoValueVSDP",err,error)
!    EXITS("WriteStringFmtTwoValueVSDP")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtTwoValueVSDP

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted varying string firstValue and the the secondString followed by a formatted integer secondValue to the given output stream specified by ID. firstFormat is used to format the first value and secondFormat is used to format the second value.
  SUBROUTINE WriteStringFmtTwoValueVSIntg(id,firstString,firstValue,firstFormat,secondString,secondValue,secondFormat,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    INTEGER(INTG), INTENT(IN) :: secondValue !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondFormat !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringFmtTwoValueVSIntg",err,error,*999)
        
    localString=firstString//firstValue//secondString//NumberToVString(secondValue,secondFormat,err,error,adjust=.FALSE.)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtTwoValueVSIntg")
    RETURN
999 ERRORS("WriteStringFmtTwoValueVSIntg",err,error)
!    EXITS("WriteStringFmtTwoValueVSIntg")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtTwoValueVSIntg

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted varying string firstValue and the the secondString followed by a formatted logical secondValue to the given output stream specified by ID. firstFormat is used to format the first value and secondFormat is used to format the second value.
  SUBROUTINE WriteStringFmtTwoValueVSL(id,firstString,firstValue,firstFormat,secondString,secondValue,secondFormat,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    LOGICAL, INTENT(IN) :: secondValue !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondFormat !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringFmtTwoValueVSL",err,error,*999)
        
    localString=firstString//firstValue//secondString//LogicalToVString(secondValue,err,error)
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtTwoValueVSL")
    RETURN
999 ERRORS("WriteStringFmtTwoValueVSL",err,error)
!    EXITS("WriteStringFmtTwoValueVSL")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtTwoValueVSL

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted varying string firstValue and the the secondString followed by a formatted single precision secondValue to the given output stream specified by ID. firstFormat is used to format the first value and secondFormat is used to format the second value.
  SUBROUTINE WriteStringFmtTwoValueVSSP(id,firstString,firstValue,firstFormat,secondString,secondValue,secondFormat,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    REAL(SP), INTENT(IN) :: secondValue !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondFormat !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

 !   ENTERS("WriteStringFmtTwoValueVSSP",err,error,*999)
        
    localString=firstString//firstValue//secondString//NumberToVString(secondValue,secondFormat,err,error,adjust=.FALSE.)
    IF(err/=0) GOTO 999
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtTwoValueVSSP")
    RETURN
999 ERRORS("WriteStringFmtTwoValueVSSP",err,error)
!    EXITS("WriteStringFmtTwoValueVSSP")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtTwoValueVSSP

  !
  !================================================================================================================================
  !

  !>Writes the firstString followed by a formatted varying string firstValue and the the secondString followed by a formatted varying string secondValue to the given output stream specified by ID. firstFormat is used to format the first value and secondFormat is used to format the second value.
  SUBROUTINE WriteStringFmtTwoValueVSVS(id,firstString,firstValue,firstFormat,secondString,secondValue,secondFormat,err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: firstString !<The first string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: firstValue !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: secondString !<The second string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: secondValue !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: secondFormat !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localString

!    ENTERS("WriteStringFmtTwoValueVSVS",err,error,*999)
        
    localString=firstString//firstValue//secondString//secondValue
    WRITE(outputString,'(A)') CHAR(localString)
    CALL WriteStr(id,err,error,*999)
      
!    EXITS("WriteStringFmtTwoValueVSVS")
    RETURN
999 ERRORS("WriteStringFmtTwoValueVSVS",err,error)
!    EXITS("WriteStringFmtTwoValueVSVS")
    RETURN 1
    
  END SUBROUTINE WriteStringFmtTwoValueVSVS

  !
  !================================================================================================================================
  !

  !>Writes the given double precision vector to the given output stream specified by ID. The firstFormat is the format initially used, followed by the repeatFormat which is repeated as many times as necessary. numberFirst is the number of data items in the firstFormat and numberRepeat is the number of data items in the repeatFormat. firstIdx and lastIdx are the extents of the data and delta is the number of indices to skip for each index.
  SUBROUTINE WriteStringVectorDP(id,firstIdx,delta,lastIdx,numberFirst,numberRepeat,vector,firstFormat,repeatFormat,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    INTEGER(INTG), INTENT(IN) :: firstIdx !<The first index of the vector to output
    INTEGER(INTG), INTENT(IN) :: delta !<The delta increment to be used when outputing the first through to the last vector index
    INTEGER(INTG), INTENT(IN) :: lastIdx !<The last index of the vector to output
    INTEGER(INTG), INTENT(IN) :: numberFirst !<The number of vector elements to be output on the first line
    INTEGER(INTG), INTENT(IN) :: numberRepeat !<The number of vector elements to be output on the second and subsequently repeated lines
    REAL(DP), INTENT(IN) :: vector(:) !<The vector to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used for the first line of output
    CHARACTER(LEN=*), INTENT(IN) :: repeatFormat !<The format type to be used for the second and subsequently repeated lines of output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) ::  current,final,count

!    ENTERS("WriteStringVectorDP",err,error,*999)
        
    current=firstIdx
    final=current+(numberFirst-1)*delta
    IF(final>lastIdx) final=lastIdx

    WRITE(outputString,FMT=firstFormat) (vector(count),count=current,final,delta)
    CALL WriteStr(id,err,error,*999)
    DO WHILE(final<lastIdx) !more stuff to do
      current=final+delta
      final=final+numberRepeat*delta
      IF(final>lastIdx) final=lastIdx
      WRITE(outputString,FMT=repeatFormat) (vector(count),count=current,final,delta)
      CALL WriteStr(id,err,error,*999)
    ENDDO !final<lastIdx

!    EXITS("WriteStringVectorDP")
    RETURN
999 ERRORS("WriteStringVectorDP",err,error)
!    EXITS("WriteStringVectorDP")
    RETURN 1
    
  END SUBROUTINE WriteStringVectorDP

  !
  !================================================================================================================================
  !

  !>Writes the given integer vector to the given output stream specified by ID. The firstFormat is the format initially used, followed by the repeatFormat which is repeated as many times as necessary. numberFirst is the number of data items in the firstFormat and numberRepeat is the number of data items in the repeatFormat. firstIdx and lastIdx are the extents of the data and delta is the number of indices to skip for each index.
  SUBROUTINE WriteStringVectorIntg(id,firstIdx,delta,lastIdx,numberFirst,numberRepeat,vector,firstFormat,repeatFormat,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    INTEGER(INTG), INTENT(IN) :: firstIdx !<The first index of the vector to output
    INTEGER(INTG), INTENT(IN) :: delta !<The delta increment to be used when outputing the first through to the last vector index
    INTEGER(INTG), INTENT(IN) :: lastIdx !<The last index of the vector to output
    INTEGER(INTG), INTENT(IN) :: numberFirst !<The number of vector elements to be output on the first line
    INTEGER(INTG), INTENT(IN) :: numberRepeat !<The number of vector elements to be output on the second and subsequently repeated lines
    INTEGER(INTG), INTENT(IN) :: vector(:) !<The vector to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used for the first line of output
    CHARACTER(LEN=*), INTENT(IN) :: repeatFormat !<The format type to be used for the second and subsequently repeated lines of output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) ::  current,final,count

!    ENTERS("WriteStringVectorIntg",err,error,*999)
        
    current=firstIdx
    final=current+(numberFirst-1)*delta
    IF(final>lastIdx) final=lastIdx
    WRITE(outputString,FMT=firstFormat) (vector(count),count=current,final,delta)
    CALL WriteStr(id,err,error,*999)
    DO WHILE(final<lastIdx) !more stuff to do
      current=final+delta
      final=final+numberRepeat*delta
      IF(final>lastIdx) final=lastIdx
      WRITE(outputString,FMT=repeatFormat) (vector(count),count=current,final,delta)
      CALL WriteStr(id,err,error,*999)
    ENDDO !final<lastIdx

!    EXITS("WriteStringVectorIntg")
    RETURN
999 ERRORS("WriteStringVectorIntg",err,error)
!    EXITS("WriteStringVectorIntg")
    RETURN 1
    
  END SUBROUTINE WriteStringVectorIntg

  !
  !================================================================================================================================
  !

  !>Writes the given integer vector to the given output stream specified by ID. The firstFormat is the format initially used, followed by the repeatFormat which is repeated as many times as necessary. numberFirst is the number of data items in the firstFormat and numberRepeat is the number of data items in the repeatFormat. firstIdx and lastIdx are the extents of the data and delta is the number of indices to skip for each index.
  SUBROUTINE WriteStringVectorLIntg(id,firstIdx,delta,lastIdx,numberFirst,numberRepeat,vector,firstFormat,repeatFormat,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    INTEGER(INTG), INTENT(IN) :: firstIdx !<The first index of the vector to output
    INTEGER(INTG), INTENT(IN) :: delta !<The delta increment to be used when outputing the first through to the last vector index
    INTEGER(INTG), INTENT(IN) :: lastIdx !<The last index of the vector to output
    INTEGER(INTG), INTENT(IN) :: numberFirst !<The number of vector elements to be output on the first line
    INTEGER(INTG), INTENT(IN) :: numberRepeat !<The number of vector elements to be output on the second and subsequently repeated lines
    INTEGER(LINTG), INTENT(IN) :: vector(:) !<The vector to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used for the first line of output
    CHARACTER(LEN=*), INTENT(IN) :: repeatFormat !<The format type to be used for the second and subsequently repeated lines of output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) ::  current,final,count

!    ENTERS("WriteStringVectorLIntg",err,error,*999)
        
    current=firstIdx
    final=current+(numberFirst-1)*delta
    IF(final>lastIdx) final=lastIdx
    WRITE(outputString,FMT=firstFormat) (vector(count),count=current,final,delta)
    CALL WriteStr(id,err,error,*999)
    DO WHILE(final<lastIdx) !more stuff to do
      current=final+delta
      final=final+numberRepeat*delta
      IF(final>lastIdx) final=lastIdx
      WRITE(outputString,FMT=repeatFormat) (vector(count),count=current,final,delta)
      CALL WriteStr(id,err,error,*999)
    ENDDO !final<lastIdx

!    EXITS("WriteStringVectorLIntg")
    RETURN
999 ERRORS("WriteStringVectorLIntg",err,error)
!    EXITS("WriteStringVectorLIntg")
    RETURN 1
    
  END SUBROUTINE WriteStringVectorLIntg

  !
  !================================================================================================================================
  !

  !>Writes the given logical vector to the given output stream specified by ID. The firstFormat is the format initially used, followed by the repeatFormat which is repeated as many times as necessary. numberFirst is the number of data items in the firstFormat and numberRepeat is the number of data items in the repeatFormat. firstIdx and lastIdx are the extents of the data and delta is the number of indices to skip for each index.
  SUBROUTINE WriteStringVectorL(id,firstIdx,delta,lastIdx,numberFirst,numberRepeat,vector,firstFormat,repeatFormat,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    INTEGER(INTG), INTENT(IN) :: firstIdx !<The first index of the vector to output
    INTEGER(INTG), INTENT(IN) :: delta !<The delta increment to be used when outputing the first through to the last vector index
    INTEGER(INTG), INTENT(IN) :: lastIdx !<The last index of the vector to output
    INTEGER(INTG), INTENT(IN) :: numberFirst !<The number of vector elements to be output on the first line
    INTEGER(INTG), INTENT(IN) :: numberRepeat !<The number of vector elements to be output on the second and subsequently repeated lines
    LOGICAL, INTENT(IN) :: vector(:) !<The vector to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used for the first line of output
    CHARACTER(LEN=*), INTENT(IN) :: repeatFormat !<The format type to be used for the second and subsequently repeated lines of output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) ::  current,final,count

!    ENTERS("WriteStringVectorL",err,error,*999)
        
    current=firstIdx
    final=current+(numberFirst-1)*delta
    IF(final>lastIdx) final=lastIdx
    WRITE(outputString,FMT=firstFormat) (vector(count),count=current,final,delta)
    CALL WriteStr(id,err,error,*999)
    DO WHILE(final<lastIdx) !more stuff to do
      current=final+delta
      final=final+numberRepeat*delta
      IF(final>lastIdx) final=lastIdx
      WRITE(outputString,FMT=repeatFormat) (vector(count),count=current,final,delta)
      CALL WriteStr(id,err,error,*999)
    ENDDO !final<lastIdx

!    EXITS("WriteStringVectorL")
    RETURN
999 ERRORS("WriteStringVectorL",err,error)
!    EXITS("WriteStringVectorL")
    RETURN 1
    
  END SUBROUTINE WriteStringVectorL

  !
  !================================================================================================================================
  !

  !>Writes the given single precision vector to the given output stream specified by ID. The firstFormat is the format initially used, followed by the repeatFormat which is repeated as many times as necessary. numberFirst is the number of data items in the firstFormat and numberRepeat is the number of data items in the repeatFormat. firstIdx and lastIdx are the extents of the data and delta is the number of indices to skip for each index.
  SUBROUTINE WriteStringVectorSP(id,firstIdx,delta,lastIdx,numberFirst,numberRepeat,vector,firstFormat,repeatFormat,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    INTEGER(INTG), INTENT(IN) :: firstIdx !<The first index of the vector to output
    INTEGER(INTG), INTENT(IN) :: delta !<The delta increment to be used when outputing the first through to the last vector index
    INTEGER(INTG), INTENT(IN) :: lastIdx !<The last index of the vector to output
    INTEGER(INTG), INTENT(IN) :: numberFirst !<The number of vector elements to be output on the first line
    INTEGER(INTG), INTENT(IN) :: numberRepeat !<The number of vector elements to be output on the second and subsequently repeated lines
    REAL(SP), INTENT(IN) :: vector(:) !<The vector to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used for the first line of output
    CHARACTER(LEN=*), INTENT(IN) :: repeatFormat !<The format type to be used for the second and subsequently repeated lines of output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) ::  current,final,count

!    ENTERS("WriteStringVectorSP",err,error,*999)
        
    current=firstIdx
    final=current+(numberFirst-1)*delta
    IF(final>lastIdx) final=lastIdx
    WRITE(outputString,FMT=firstFormat) (vector(count),count=current,final,delta)
    CALL WriteStr(id,err,error,*999)
    DO WHILE(final<lastIdx) !more stuff to do
      current=final+delta
      final=final+numberRepeat*delta
      IF(final>lastIdx) final=lastIdx
      WRITE(outputString,FMT=repeatFormat) (vector(count),count=current,final,delta)
      CALL WriteStr(id,err,error,*999)
    ENDDO !final<lastIdx

!    EXITS("WriteStringVectorSP")
    RETURN
999 ERRORS("WriteStringVectorSP",err,error)
!    EXITS("WriteStringVectorSP")
    RETURN 1
    
  END SUBROUTINE WriteStringVectorSP

  !
  !================================================================================================================================
  !

  !>Writes the given indexed double precision vector to the given output stream specified by ID. numIndices is the number of indices and indices(i) contain the indices of the vector to write. The firstFormat is the format initially used, followed by the repeatFormat which is repeated as many times as necessary. numberFirst is the number of data items in the firstFormat and numberRepeat is the number of data items in the repeatFormat. delta is the number of actual indices to skip for each index.
  SUBROUTINE WriteStringIdxVectorDP(id,numIndices,indices,delta,numberFirst,numberRepeat,vector,firstFormat,repeatFormat, &
    & err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    INTEGER(INTG), INTENT(IN) :: numIndices !<The number of indices of the vector to output
    INTEGER(INTG), INTENT(IN) :: indices(numIndices) !<indices(i). The i'th index of the vector to output
    INTEGER(INTG), INTENT(IN) :: delta !<The delta increment to be used when outputing the first through to the last vector index
    INTEGER(INTG), INTENT(IN) :: numberFirst !<The number of vector elements to be output on the first line
    INTEGER(INTG), INTENT(IN) :: numberRepeat !<The number of vector elements to be output on the second and subsequently repeated lines
    REAL(DP), INTENT(IN) :: vector(:) !<The vector to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used for the first line of output
    CHARACTER(LEN=*), INTENT(IN) :: repeatFormat !<The format type to be used for the second and subsequently repeated lines of output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) ::  current,count,numberToDo

!    ENTERS("WriteStringIdxVectorDP",err,error,*999)
        
    numberToDo=numIndices
    WRITE(outputString,FMT=firstFormat) (vector((indices(count)-1)*delta+1),count=1,MIN(numberFirst,numIndices))
    CALL WriteStr(id,err,error,*999)
    numberToDo=numIndices-numberFirst
    current=numberFirst+1
    DO WHILE(numberToDo>0) !more stuff to do
      WRITE(outputString,FMT=repeatFormat) (vector((indices(count)-1)*delta+1),count=current,MIN(current+numberRepeat-1, &
        & numIndices))
      CALL WriteStr(id,err,error,*999)
      current=current+numberRepeat
      numberToDo=numberToDo-numberRepeat
    ENDDO !numberToDo > 0

!    EXITS("WriteStringIdxVectorDP")
    RETURN
999 ERRORS("WriteStringIdxVectorDP",err,error)
!    EXITS("WriteStringIdxVectorDP")
    RETURN 1
    
  END SUBROUTINE WriteStringIdxVectorDP

  !
  !================================================================================================================================
  !

  !>Writes the given indexed integer vector to the given output stream specified by ID. numIndices is the number of indices and indices(i) contain the indices of the vector to write. The firstFormat is the format initially used, followed by the repeatFormat which is repeated as many times as necessary. numberFirst is the number of data items in the firstFormat and numberRepeat is the number of data items in the repeatFormat. delta is the number of actual indices to skip for each index.
  SUBROUTINE WriteStringIdxVectorIntg(id,numIndices,indices,delta,numberFirst,numberRepeat,vector,firstFormat,repeatFormat, &
    & err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    INTEGER(INTG), INTENT(IN) :: numIndices !<The number of indices of the vector to output
    INTEGER(INTG), INTENT(IN) :: indices(numIndices) !<indices(i). The i'th index of the vector to output
    INTEGER(INTG), INTENT(IN) :: delta !<The delta increment to be used when outputing the first through to the last vector index
    INTEGER(INTG), INTENT(IN) :: numberFirst !<The number of vector elements to be output on the first line
    INTEGER(INTG), INTENT(IN) :: numberRepeat !<The number of vector elements to be output on the second and subsequently repeated lines
    INTEGER(INTG), INTENT(IN) :: vector(:) !<The vector to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used for the first line of output
    CHARACTER(LEN=*), INTENT(IN) :: repeatFormat !<The format type to be used for the second and subsequently repeated lines of output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) ::  current,count,numberToDo

!    ENTERS("WriteStringIdxVectorIntg",err,error,*999)
        
    numberToDo=numIndices
    WRITE(outputString,FMT=firstFormat) (vector((indices(count)-1)*delta+1),count=1,MIN(numberFirst,numIndices))
    CALL WriteStr(id,err,error,*999)
    numberToDo=numIndices-numberFirst
    current=numberFirst+1
    DO WHILE(numberToDo>0) !more stuff to do
      WRITE(outputString,FMT=repeatFormat) (vector((indices(count)-1)*delta+1),count=current,MIN(current+numberRepeat-1, &
        & numIndices))
      CALL WriteStr(id,err,error,*999)
      current=current+numberRepeat
      numberToDo=numberToDo-numberRepeat
    ENDDO !numberToDo > 0

!    EXITS("WriteStringIdxVectorIntg")
    RETURN
999 ERRORS("WriteStringIdxVectorIntg",err,error)
!    EXITS("WriteStringIdxVectorIntg")
    RETURN 1
    
  END SUBROUTINE WriteStringIdxVectorIntg

  !
  !================================================================================================================================
  !

  !>Writes the given indexed integer vector to the given output stream specified by ID. numIndices is the number of indices and indices(i) contain the indices of the vector to write. The firstFormat is the format initially used, followed by the repeatFormat which is repeated as many times as necessary. numberFirst is the number of data items in the firstFormat and numberRepeat is the number of data items in the repeatFormat. delta is the number of actual indices to skip for each index.
  SUBROUTINE WriteStringIdxVectorLIntg(id,numIndices,indices,delta,numberFirst,numberRepeat,vector,firstFormat,repeatFormat, &
    & err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    INTEGER(INTG), INTENT(IN) :: numIndices !<The number of indices of the vector to output
    INTEGER(INTG), INTENT(IN) :: indices(numIndices) !<indices(i). The i'th index of the vector to output
    INTEGER(INTG), INTENT(IN) :: delta !<The delta increment to be used when outputing the first through to the last vector index
    INTEGER(INTG), INTENT(IN) :: numberFirst !<The number of vector elements to be output on the first line
    INTEGER(INTG), INTENT(IN) :: numberRepeat !<The number of vector elements to be output on the second and subsequently repeated lines
    INTEGER(LINTG), INTENT(IN) :: vector(:) !<The vector to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used for the first line of output
    CHARACTER(LEN=*), INTENT(IN) :: repeatFormat !<The format type to be used for the second and subsequently repeated lines of output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) ::  current,count,numberToDo

!    ENTERS("WriteStringIdxVectorLIntg",err,error,*999)
        
    numberToDo=numIndices
    WRITE(outputString,FMT=firstFormat) (vector((indices(count)-1)*delta+1),count=1,MIN(numberFirst,numIndices))
    CALL WriteStr(id,err,error,*999)
    numberToDo=numIndices-numberFirst
    current=numberFirst+1
    DO WHILE(numberToDo>0) !more stuff to do
      WRITE(outputString,FMT=repeatFormat) (vector((indices(count)-1)*delta+1),count=current,MIN(current+numberRepeat-1, &
        & numIndices))
      CALL WriteStr(id,err,error,*999)
      current=current+numberRepeat
      numberToDo=numberToDo-numberRepeat
    ENDDO !numberToDo > 0

!    EXITS("WriteStringIdxVectorLIntg")
    RETURN
999 ERRORS("WriteStringIdxVectorLIntg",err,error)
!    EXITS("WriteStringIdxVectorLIntg")
    RETURN 1
    
  END SUBROUTINE WriteStringIdxVectorLIntg

  !
  !================================================================================================================================
  !

  !>Writes the given indexed logical vector to the given output stream specified by ID. numIndices is the number of indices and indices(i) contain the indices of the vector to write. The firstFormat is the format initially used, followed by the repeatFormat which is repeated as many times as necessary. numberFirst is the number of data items in the firstFormat and numberRepeat is the number of data items in the repeatFormat. delta is the number of actual indices to skip for each index.
  SUBROUTINE WriteStringIdxVectorL(id,numIndices,indices,delta,numberFirst,numberRepeat,vector,firstFormat,repeatFormat, &
    & err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    INTEGER(INTG), INTENT(IN) :: numIndices !<The number of indices of the vector to output
    INTEGER(INTG), INTENT(IN) :: indices(numIndices) !<indices(i). The i'th index of the vector to output
    INTEGER(INTG), INTENT(IN) :: delta !<The delta increment to be used when outputing the first through to the last vector index
    INTEGER(INTG), INTENT(IN) :: numberFirst !<The number of vector elements to be output on the first line
    INTEGER(INTG), INTENT(IN) :: numberRepeat !<The number of vector elements to be output on the second and subsequently repeated lines
    LOGICAL, INTENT(IN) :: vector(:) !<The vector to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used for the first line of output
    CHARACTER(LEN=*), INTENT(IN) :: repeatFormat !<The format type to be used for the second and subsequently repeated lines of output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) ::  current,count,numberToDo

!    ENTERS("WriteStringIdxVectorL",err,error,*999)
        
    numberToDo=numIndices
    WRITE(outputString,FMT=firstFormat) (vector((indices(count)-1)*delta+1),count=1,MIN(numberFirst,numIndices))
    CALL WriteStr(id,err,error,*999)
    numberToDo=numIndices-numberFirst
    current=numberFirst+1
    DO WHILE(numberToDo>0) !more stuff to do
      WRITE(outputString,FMT=repeatFormat) (vector((indices(count)-1)*delta+1),count=current,MIN(current+numberRepeat-1, &
        & numIndices))
      CALL WriteStr(id,err,error,*999)
      current=current+numberRepeat
      numberToDo=numberToDo-numberRepeat
    ENDDO !numberToDo > 0

!    EXITS("WriteStringIdxVectorL")
    RETURN
999 ERRORS("WriteStringIdxVectorL",err,error)
!    EXITS("WriteStringIdxVectorL")
    RETURN 1
    
  END SUBROUTINE WriteStringIdxVectorL

  !
  !================================================================================================================================
  !

  !>Writes the given indexed single precision vector to the given output stream specified by ID. numIndices is the number of indices and indices(i) contain the indices of the vector to write. The firstFormat is the format initially used, followed by the repeatFormat which is repeated as many times as necessary. numberFirst is the number of data items in the firstFormat and numberRepeat is the number of data items in the repeatFormat. delta is the number of actual indices to skip for each index.
  SUBROUTINE WriteStringIdxVectorSP(id,numIndices,indices,delta,numberFirst,numberRepeat,vector,firstFormat,repeatFormat, &
    & err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    INTEGER(INTG), INTENT(IN) :: numIndices !<The number of indices of the vector to output
    INTEGER(INTG), INTENT(IN) :: indices(numIndices) !<indices(i). The i'th index of the vector to output
    INTEGER(INTG), INTENT(IN) :: delta !<The delta increment to be used when outputing the first through to the last vector index
    INTEGER(INTG), INTENT(IN) :: numberFirst !<The number of vector elements to be output on the first line
    INTEGER(INTG), INTENT(IN) :: numberRepeat !<The number of vector elements to be output on the second and subsequently repeated lines
    REAL(SP), INTENT(IN) :: vector(:) !<The vector to be output
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used for the first line of output
    CHARACTER(LEN=*), INTENT(IN) :: repeatFormat !<The format type to be used for the second and subsequently repeated lines of output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) ::  current,count,numberToDo

!    ENTERS("WriteStringIdxVectorSP",err,error,*999)
        
    numberToDo=numIndices
    WRITE(outputString,FMT=firstFormat) (vector((indices(count)-1)*delta+1),count=1,MIN(numberFirst,numIndices))
    CALL WriteStr(id,err,error,*999)
    numberToDo=numIndices-numberFirst
    current=numberFirst+1
    DO WHILE(numberToDo>0) !more stuff to do
      WRITE(outputString,FMT=repeatFormat) (vector((indices(count)-1)*delta+1),count=current,MIN(current+numberRepeat-1, &
        & numIndices))
      CALL WriteStr(id,err,error,*999)
      current=current+numberRepeat
      numberToDo=numberToDo-numberRepeat
    ENDDO !numberToDo > 0

!    EXITS("WriteStringIdxVectorSP")
    RETURN
999 ERRORS("WriteStringIdxVectorSP",err,error)
!    EXITS("WriteStringIdxVectorSP")
    RETURN 1
    
  END SUBROUTINE WriteStringIdxVectorSP

  !
  !================================================================================================================================
  !

  !>Writes the given double precision matrix to the given output stream specified by ID. The basic output is determined by the flag indexFormatType. If indexFormatType is WRITE_STRING_MATRIX_NAME_ONLY then the first line of output for each row is matrixNameFormat concatenated named with the firstFormat. If indexFormatType is WRITE_STRING_MATRIX_NAME_AND_INDICES then the first line of output for each row is matrixNameFormat concatenated with rowIndexFormat and concatenated with firstFormat. Note that with a WRITE_STRING_MATRIX_NAME_AND_INDICES index format type the row number will be supplied to the format before the matrix data. The firstFormat is the format initially used, followed by the repeatFormat which is repeated as many times as necessary. numberFirst is the number of data items in the firstFormat and numberRepeat is the number of data items in the repeatFormat. firstRow/firstColumn and lastRow/lastColumn are the extents of the row/column and deltaRow/deltaColumn is the number of indices to skip for each row/column index.
  SUBROUTINE WriteStringMatrixDP(id,firstRow,deltaRow,lastRow,firstColumn,deltaColumn,lastColumn,numberFirst, &
    & numberRepeat,matrix,indexFormatType,matrixNameFormat,rowIndexFormat,firstFormat,repeatFormat,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    INTEGER(INTG), INTENT(IN) :: firstRow !<The first row of the matrix to be output
    INTEGER(INTG), INTENT(IN) :: deltaRow !<The delta row increment to be used when outputing the first through to the last matrix row
    INTEGER(INTG), INTENT(IN) :: lastRow !<The last row of the matrix to be output
    INTEGER(INTG), INTENT(IN) :: firstColumn !<The first column of the matrix to be output
    INTEGER(INTG), INTENT(IN) :: deltaColumn !<The delta column increate to be used when outputing the first through to the last matrix column
    INTEGER(INTG), INTENT(IN) :: lastColumn !<The last column of the matrix to be output
    INTEGER(INTG), INTENT(IN) :: numberFirst !<The number of matrix elements to be output on the first line
    INTEGER(INTG), INTENT(IN) :: numberRepeat !<The number of matrix elements to be output on the second and subsequently repeated lines
    REAL(DP), INTENT(IN) :: matrix(:,:) !<The matrix to be output
    INTEGER(INTG), INTENT(IN) :: indexFormatType !<The format type to be used for the matrix name and indices \see InputOutput_MatrixNameIndexFormat,InputOutput::MatrixNameIndexFormat
    CHARACTER(LEN=*), INTENT(IN) :: matrixNameFormat !<The format string to be used to format the matrix name
    CHARACTER(LEN=*), INTENT(IN) :: rowIndexFormat !<The format string to be used to format the row indices
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used for the first line of output
    CHARACTER(LEN=*), INTENT(IN) :: repeatFormat !<The format type to be used for the second and subsequently repeated lines of output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) ::  currentRow,currentColumn,finalColumn,count
    CHARACTER(LEN=MAXSTRLEN) :: formatStr

!    ENTERS("WriteStringMatrixDP",err,error,*999)

    IF(indexFormatType==WRITE_STRING_MATRIX_NAME_ONLY) THEN
      formatStr=matrixNameFormat//firstFormat
    ELSE IF(indexFormatType==WRITE_STRING_MATRIX_NAME_AND_INDICES) THEN
      formatStr=matrixNameFormat//rowIndexFormat//firstFormat
    ELSE
      CALL FlagError("Invalid index format type.",err,error,*999)
    ENDIF
    DO currentRow=firstRow,lastRow,deltaRow
      currentColumn=firstColumn
      finalColumn=currentColumn+(numberFirst-1)*deltaColumn
      IF(finalColumn>lastColumn) finalColumn=lastColumn
      IF(indexFormatType==WRITE_STRING_MATRIX_NAME_ONLY) THEN
        WRITE(outputString,FMT=formatStr) (matrix(currentRow,count),count=currentColumn,finalColumn,deltaColumn)
      ELSE IF(indexFormatType==WRITE_STRING_MATRIX_NAME_AND_INDICES) THEN
        WRITE(outputString,FMT=formatStr) currentRow,(matrix(currentRow,count),count=currentColumn,finalColumn,deltaColumn)
      ENDIF
      CALL WriteStr(id,err,error,*999)
      DO WHILE(finalColumn<lastColumn) !more stuff to do
        currentColumn=finalColumn+deltaColumn
        finalColumn=finalColumn+numberRepeat*deltaColumn
        IF(finalColumn>lastColumn) finalColumn=lastColumn
        WRITE(outputString,FMT=repeatFormat) (matrix(currentRow,count),count=currentColumn,finalColumn,deltaColumn)
        CALL WriteStr(id,err,error,*999)
      ENDDO !finalColumnn<lastColumn
    ENDDO !currentRow
    
!    EXITS("WriteStringMatrixDP")
    RETURN
999 ERRORS("WriteStringMatrixDP",err,error)
!    EXITS("WriteStringMatrixDP")
    RETURN 1
    
  END SUBROUTINE WriteStringMatrixDP

  !
  !================================================================================================================================
  !

  !>Writes the given integer matrix to the given output stream specified by ID. The basic output is determined by the flag indexFormatType. If indexFormatType is WRITE_STRING_MATRIX_NAME_ONLY then the first line of output for each row is matrixNameFormat concatenated named with the firstFormat. If indexFormatType is WRITE_STRING_MATRIX_NAME_AND_INDICES then the first line of output for each row is matrixNameFormat concatenated with rowIndexFormat and concatenated with firstFormat. Note that with a WRITE_STRING_MATRIX_NAME_AND_INDICES index format type the row number will be supplied to the format before the matrix data. The firstFormat is the format initially used, followed by the repeatFormat which is repeated as many times as necessary. numberFirst is the number of data items in the firstFormat and numberRepeat is the number of data items in the repeatFormat. firstRow/firstColumn and lastRow/lastColumn are the extents of the row/column and deltaRow/deltaColumn is the number of indices to skip for each row/column index.
  SUBROUTINE WriteStringMatrixIntg(id,firstRow,deltaRow,lastRow,firstColumn,deltaColumn,lastColumn,numberFirst, &
    & numberRepeat,matrix,indexFormatType,matrixNameFormat,rowIndexFormat,firstFormat,repeatFormat,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    INTEGER(INTG), INTENT(IN) :: firstRow !<The first row of the matrix to be output
    INTEGER(INTG), INTENT(IN) :: deltaRow !<The delta row increment to be used when outputing the first through to the last matrix row
    INTEGER(INTG), INTENT(IN) :: lastRow !<The last row of the matrix to be output
    INTEGER(INTG), INTENT(IN) :: firstColumn !<The first column of the matrix to be output
    INTEGER(INTG), INTENT(IN) :: deltaColumn !<The delta column increate to be used when outputing the first through to the last matrix column
    INTEGER(INTG), INTENT(IN) :: lastColumn !<The last column of the matrix to be output
    INTEGER(INTG), INTENT(IN) :: numberFirst !<The number of matrix elements to be output on the first line
    INTEGER(INTG), INTENT(IN) :: numberRepeat !<The number of matrix elements to be output on the second and subsequently repeated lines
    INTEGER(INTG), INTENT(IN) :: matrix(:,:) !<The matrix to be output
    INTEGER(INTG), INTENT(IN) :: indexFormatType !<The format type to be used for the matrix name and indices \see InputOutput_MatrixNameIndexFormat,InputOutput::MatrixNameIndexFormat
    CHARACTER(LEN=*), INTENT(IN) :: matrixNameFormat !<The format string to be used to format the matrix name
    CHARACTER(LEN=*), INTENT(IN) :: rowIndexFormat !<The format string to be used to format the row indices
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used for the first line of output
    CHARACTER(LEN=*), INTENT(IN) :: repeatFormat !<The format type to be used for the second and subsequently repeated lines of output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) ::  currentRow,currentColumn,finalColumn,count
    CHARACTER(LEN=MAXSTRLEN) :: formatStr

!    ENTERS("WriteStringMatrixIntg",err,error,*999)

    IF(indexFormatType==WRITE_STRING_MATRIX_NAME_ONLY) THEN
      formatStr=matrixNameFormat//firstFormat
    ELSE IF(indexFormatType==WRITE_STRING_MATRIX_NAME_AND_INDICES) THEN
      formatStr=matrixNameFormat//rowIndexFormat//firstFormat
    ELSE
      CALL FlagError("Invalid index format type.",err,error,*999)
    ENDIF
    DO currentRow=firstRow,lastRow,deltaRow
      currentColumn=firstColumn
      finalColumn=currentColumn+(numberFirst-1)*deltaColumn
      IF(finalColumn>lastColumn) finalColumn=lastColumn
      IF(indexFormatType==WRITE_STRING_MATRIX_NAME_ONLY) THEN
        WRITE(outputString,FMT=formatStr) (matrix(currentRow,count),count=currentColumn,finalColumn,deltaColumn)
      ELSE IF(indexFormatType==WRITE_STRING_MATRIX_NAME_AND_INDICES) THEN
        WRITE(outputString,FMT=formatStr) currentRow,(matrix(currentRow,count),count=currentColumn,finalColumn,deltaColumn)
      ENDIF
      CALL WriteStr(id,err,error,*999)
      DO WHILE(finalColumn<lastColumn) !more stuff to do
        currentColumn=finalColumn+deltaColumn
        finalColumn=finalColumn+numberRepeat*deltaColumn
        IF(finalColumn>lastColumn) finalColumn=lastColumn
        WRITE(outputString,FMT=repeatFormat) (matrix(currentRow,count),count=currentColumn,finalColumn,deltaColumn)
        CALL WriteStr(id,err,error,*999)
      ENDDO !finalColumnn<lastColumn
    ENDDO !currentRow
    
!    EXITS("WriteStringMatrixIntg")
    RETURN
999 ERRORS("WriteStringMatrixIntg",err,error)
!    EXITS("WriteStringMatrixIntg")
    RETURN 1
    
  END SUBROUTINE WriteStringMatrixIntg

  !
  !================================================================================================================================
  !

  !>Writes the given long integer matrix to the given output stream specified by ID. The basic output is determined by the flag indexFormatType. If indexFormatType is WRITE_STRING_MATRIX_NAME_ONLY then the first line of output for each row is matrixNameFormat concatenated named with the firstFormat. If indexFormatType is WRITE_STRING_MATRIX_NAME_AND_INDICES then the first line of output for each row is matrixNameFormat concatenated with rowIndexFormat and concatenated with firstFormat. Note that with a WRITE_STRING_MATRIX_NAME_AND_INDICES index format type the row number will be supplied to the format before the matrix data. The firstFormat is the format initially used, followed by the repeatFormat which is repeated as many times as necessary. numberFirst is the number of data items in the firstFormat and numberRepeat is the number of data items in the repeatFormat. firstRow/firstColumn and lastRow/lastColumn are the extents of the row/column and deltaRow/deltaColumn is the number of indices to skip for each row/column index.
  SUBROUTINE WriteStringMatrixLIntg(id,firstRow,deltaRow,lastRow,firstColumn,deltaColumn,lastColumn,numberFirst, &
    & numberRepeat,matrix,indexFormatType,matrixNameFormat,rowIndexFormat,firstFormat,repeatFormat,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    INTEGER(INTG), INTENT(IN) :: firstRow !<The first row of the matrix to be output
    INTEGER(INTG), INTENT(IN) :: deltaRow !<The delta row increment to be used when outputing the first through to the last matrix row
    INTEGER(INTG), INTENT(IN) :: lastRow !<The last row of the matrix to be output
    INTEGER(INTG), INTENT(IN) :: firstColumn !<The first column of the matrix to be output
    INTEGER(INTG), INTENT(IN) :: deltaColumn !<The delta column increate to be used when outputing the first through to the last matrix column
    INTEGER(INTG), INTENT(IN) :: lastColumn !<The last column of the matrix to be output
    INTEGER(INTG), INTENT(IN) :: numberFirst !<The number of matrix elements to be output on the first line
    INTEGER(INTG), INTENT(IN) :: numberRepeat !<The number of matrix elements to be output on the second and subsequently repeated lines
    INTEGER(LINTG), INTENT(IN) :: matrix(:,:) !<The matrix to be output
    INTEGER(INTG), INTENT(IN) :: indexFormatType !<The format type to be used for the matrix name and indices \see InputOutput_MatrixNameIndexFormat,InputOutput::MatrixNameIndexFormat
    CHARACTER(LEN=*), INTENT(IN) :: matrixNameFormat !<The format string to be used to format the matrix name
    CHARACTER(LEN=*), INTENT(IN) :: rowIndexFormat !<The format string to be used to format the row indices
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used for the first line of output
    CHARACTER(LEN=*), INTENT(IN) :: repeatFormat !<The format type to be used for the second and subsequently repeated lines of output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) ::  currentRow,currentColumn,finalColumn,count
    CHARACTER(LEN=MAXSTRLEN) :: formatStr

!    ENTERS("WriteStringMatrixLIntg",err,error,*999)

    IF(indexFormatType==WRITE_STRING_MATRIX_NAME_ONLY) THEN
      formatStr=matrixNameFormat//firstFormat
    ELSE IF(indexFormatType==WRITE_STRING_MATRIX_NAME_AND_INDICES) THEN
      formatStr=matrixNameFormat//rowIndexFormat//firstFormat
    ELSE
      CALL FlagError("Invalid index format type.",err,error,*999)
    ENDIF
    DO currentRow=firstRow,lastRow,deltaRow
      currentColumn=firstColumn
      finalColumn=currentColumn+(numberFirst-1)*deltaColumn
      IF(finalColumn>lastColumn) finalColumn=lastColumn
      IF(indexFormatType==WRITE_STRING_MATRIX_NAME_ONLY) THEN
        WRITE(outputString,FMT=formatStr) (matrix(currentRow,count),count=currentColumn,finalColumn,deltaColumn)
      ELSE IF(indexFormatType==WRITE_STRING_MATRIX_NAME_AND_INDICES) THEN
        WRITE(outputString,FMT=formatStr) currentRow,(matrix(currentRow,count),count=currentColumn,finalColumn,deltaColumn)
      ENDIF
      CALL WriteStr(id,err,error,*999)
      DO WHILE(finalColumn<lastColumn) !more stuff to do
        currentColumn=finalColumn+deltaColumn
        finalColumn=finalColumn+numberRepeat*deltaColumn
        IF(finalColumn>lastColumn) finalColumn=lastColumn
        WRITE(outputString,FMT=repeatFormat) (matrix(currentRow,count),count=currentColumn,finalColumn,deltaColumn)
        CALL WriteStr(id,err,error,*999)
      ENDDO !finalColumnn<lastColumn
    ENDDO !currentRow
    
!    EXITS("WriteStringMatrixLIntg")
    RETURN
999 ERRORS("WriteStringMatrixLIntg",err,error)
!    EXITS("WriteStringMatrixLIntg")
    RETURN 1
    
  END SUBROUTINE WriteStringMatrixLIntg

  !
  !================================================================================================================================
  !

  !>Writes the given logical matrix to the given output stream specified by ID. The basic output is determined by the flag indexFormatType. If indexFormatType is WRITE_STRING_MATRIX_NAME_ONLY then the first line of output for each row is matrixNameFormat concatenated named with the firstFormat. If indexFormatType is WRITE_STRING_MATRIX_NAME_AND_INDICES then the first line of output for each row is matrixNameFormat concatenated with rowIndexFormat and concatenated with firstFormat. Note that with a WRITE_STRING_MATRIX_NAME_AND_INDICES index format type the row number will be supplied to the format before the matrix data. The firstFormat is the format initially used, followed by the repeatFormat which is repeated as many times as necessary. numberFirst is the number of data items in the firstFormat and numberRepeat is the number of data items in the repeatFormat. firstRow/firstColumn and lastRow/lastColumn are the extents of the row/column and deltaRow/deltaColumn is the number of indices to skip for each row/column index.
  SUBROUTINE WriteStringMatrixL(id,firstRow,deltaRow,lastRow,firstColumn,deltaColumn,lastColumn,numberFirst, &
    & numberRepeat,matrix,indexFormatType,matrixNameFormat,rowIndexFormat,firstFormat,repeatFormat,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    INTEGER(INTG), INTENT(IN) :: firstRow !<The first row of the matrix to be output
    INTEGER(INTG), INTENT(IN) :: deltaRow !<The delta row increment to be used when outputing the first through to the last matrix row
    INTEGER(INTG), INTENT(IN) :: lastRow !<The last row of the matrix to be output
    INTEGER(INTG), INTENT(IN) :: firstColumn !<The first column of the matrix to be output
    INTEGER(INTG), INTENT(IN) :: deltaColumn !<The delta column increate to be used when outputing the first through to the last matrix column
    INTEGER(INTG), INTENT(IN) :: lastColumn !<The last column of the matrix to be output
    INTEGER(INTG), INTENT(IN) :: numberFirst !<The number of matrix elements to be output on the first line
    INTEGER(INTG), INTENT(IN) :: numberRepeat !<The number of matrix elements to be output on the second and subsequently repeated lines
    LOGICAL, INTENT(IN) :: matrix(:,:) !<The matrix to be output
    INTEGER(INTG), INTENT(IN) :: indexFormatType !<The format type to be used for the matrix name and indices \see InputOutput_MatrixNameIndexFormat,InputOutput::MatrixNameIndexFormat
    CHARACTER(LEN=*), INTENT(IN) :: matrixNameFormat !<The format string to be used to format the matrix name
    CHARACTER(LEN=*), INTENT(IN) :: rowIndexFormat !<The format string to be used to format the row indices
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used for the first line of output
    CHARACTER(LEN=*), INTENT(IN) :: repeatFormat !<The format type to be used for the second and subsequently repeated lines of output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) ::  currentRow,currentColumn,finalColumn,count
    CHARACTER(LEN=MAXSTRLEN) :: formatStr

!    ENTERS("WriteStringMatrixL",err,error,*999)

    IF(indexFormatType==WRITE_STRING_MATRIX_NAME_ONLY) THEN
      formatStr=matrixNameFormat//firstFormat
    ELSE IF(indexFormatType==WRITE_STRING_MATRIX_NAME_AND_INDICES) THEN
      formatStr=matrixNameFormat//rowIndexFormat//firstFormat
    ELSE
      CALL FlagError("Invalid index format type.",err,error,*999)
    ENDIF
    DO currentRow=firstRow,lastRow,deltaRow
      currentColumn=firstColumn
      finalColumn=currentColumn+(numberFirst-1)*deltaColumn
      IF(finalColumn>lastColumn) finalColumn=lastColumn
      IF(indexFormatType==WRITE_STRING_MATRIX_NAME_ONLY) THEN
        WRITE(outputString,FMT=formatStr) (matrix(currentRow,count),count=currentColumn,finalColumn,deltaColumn)
      ELSE IF(indexFormatType==WRITE_STRING_MATRIX_NAME_AND_INDICES) THEN
        WRITE(outputString,FMT=formatStr) currentRow,(matrix(currentRow,count),count=currentColumn,finalColumn,deltaColumn)
      ENDIF
      CALL WriteStr(id,err,error,*999)
      DO WHILE(finalColumn<lastColumn) !more stuff to do
        currentColumn=finalColumn+deltaColumn
        finalColumn=finalColumn+numberRepeat*deltaColumn
        IF(finalColumn>lastColumn) finalColumn=lastColumn
        WRITE(outputString,FMT=repeatFormat) (matrix(currentRow,count),count=currentColumn,finalColumn,deltaColumn)
        CALL WriteStr(id,err,error,*999)
      ENDDO !finalColumnn<lastColumn
    ENDDO !currentRow
    
!    EXITS("WriteStringMatrixL")
    RETURN
999 ERRORS("WriteStringMatrixL",err,error)
!    EXITS("WriteStringMatrixL")
    RETURN 1
    
  END SUBROUTINE WriteStringMatrixL

  !
  !================================================================================================================================
  !

  !>Writes the given single precision matrix to the given output stream specified by ID. The basic output is determined by the flag indexFormatType. If indexFormatType is WRITE_STRING_MATRIX_NAME_ONLY then the first line of output for each row is matrixNameFormat concatenated named with the firstFormat. If indexFormatType is WRITE_STRING_MATRIX_NAME_AND_INDICES then the first line of output for each row is matrixNameFormat concatenated with rowIndexFormat and concatenated with firstFormat. Note that with a WRITE_STRING_MATRIX_NAME_AND_INDICES index format type the row number will be supplied to the format before the matrix data. The firstFormat is the format initially used, followed by the repeatFormat which is repeated as many times as necessary. numberFirst is the number of data items in the firstFormat and numberRepeat is the number of data items in the repeatFormat. firstRow/firstColumn and lastRow/lastColumn are the extents of the row/column and deltaRow/deltaColumn is the number of indices to skip for each row/column index.
  SUBROUTINE WriteStringMatrixSP(id,firstRow,deltaRow,lastRow,firstColumn,deltaColumn,lastColumn,numberFirst, &
    & numberRepeat,matrix,indexFormatType,matrixNameFormat,rowIndexFormat,firstFormat,repeatFormat,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: id !<The ID of the output stream. An ID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    INTEGER(INTG), INTENT(IN) :: firstRow !<The first row of the matrix to be output
    INTEGER(INTG), INTENT(IN) :: deltaRow !<The delta row increment to be used when outputing the first through to the last matrix row
    INTEGER(INTG), INTENT(IN) :: lastRow !<The last row of the matrix to be output
    INTEGER(INTG), INTENT(IN) :: firstColumn !<The first column of the matrix to be output
    INTEGER(INTG), INTENT(IN) :: deltaColumn !<The delta column increate to be used when outputing the first through to the last matrix column
    INTEGER(INTG), INTENT(IN) :: lastColumn !<The last column of the matrix to be output
    INTEGER(INTG), INTENT(IN) :: numberFirst !<The number of matrix elements to be output on the first line
    INTEGER(INTG), INTENT(IN) :: numberRepeat !<The number of matrix elements to be output on the second and subsequently repeated lines
    REAL(SP), INTENT(IN) :: matrix(:,:) !<The matrix to be output
    INTEGER(INTG), INTENT(IN) :: indexFormatType !<The format type to be used for the matrix name and indices \see InputOutput_MatrixNameIndexFormat,InputOutput::MatrixNameIndexFormat
    CHARACTER(LEN=*), INTENT(IN) :: matrixNameFormat !<The format string to be used to format the matrix name
    CHARACTER(LEN=*), INTENT(IN) :: rowIndexFormat !<The format string to be used to format the row indices
    CHARACTER(LEN=*), INTENT(IN) :: firstFormat !<The format string to be used for the first line of output
    CHARACTER(LEN=*), INTENT(IN) :: repeatFormat !<The format type to be used for the second and subsequently repeated lines of output
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) ::  currentRow,currentColumn,finalColumn,count
    CHARACTER(LEN=MAXSTRLEN) :: formatStr

!    ENTERS("WriteStringMatrixSP",err,error,*999)

    IF(indexFormatType==WRITE_STRING_MATRIX_NAME_ONLY) THEN
      formatStr=matrixNameFormat//firstFormat
    ELSE IF(indexFormatType==WRITE_STRING_MATRIX_NAME_AND_INDICES) THEN
      formatStr=matrixNameFormat//rowIndexFormat//firstFormat
    ELSE
      CALL FlagError("Invalid index format type.",err,error,*999)
    ENDIF
    DO currentRow=firstRow,lastRow,deltaRow
      currentColumn=firstColumn
      finalColumn=currentColumn+(numberFirst-1)*deltaColumn
      IF(finalColumn>lastColumn) finalColumn=lastColumn
      IF(indexFormatType==WRITE_STRING_MATRIX_NAME_ONLY) THEN
        WRITE(outputString,FMT=formatStr) (matrix(currentRow,count),count=currentColumn,finalColumn,deltaColumn)
      ELSE IF(indexFormatType==WRITE_STRING_MATRIX_NAME_AND_INDICES) THEN
        WRITE(outputString,FMT=formatStr) currentRow,(matrix(currentRow,count),count=currentColumn,finalColumn,deltaColumn)
      ENDIF
      CALL WriteStr(id,err,error,*999)
      DO WHILE(finalColumn<lastColumn) !more stuff to do
        currentColumn=finalColumn+deltaColumn
        finalColumn=finalColumn+numberRepeat*deltaColumn
        IF(finalColumn>lastColumn) finalColumn=lastColumn
        WRITE(outputString,FMT=repeatFormat) (matrix(currentRow,count),count=currentColumn,finalColumn,deltaColumn)
        CALL WriteStr(id,err,error,*999)
      ENDDO !finalColumnn<lastColumn
    ENDDO !currentRow
    
!    EXITS("WriteStringMatrixSP")
    RETURN
999 ERRORS("WriteStringMatrixSP",err,error)
!    EXITS("WriteStringMatrixSP")
    RETURN 1
    
  END SUBROUTINE WriteStringMatrixSP

  !
  !================================================================================================================================
  !

END MODULE InputOutput
