!> \file
!> \author Chris Bradley
!> \brief This module contains all string manipulation and transformation routines.
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

!>This module contains all string manipulation and transformation routines.
MODULE Strings

  USE BaseRoutines
  USE Constants
  USE Kinds
  USE ISO_VARYING_STRING

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Interfaces

  !>Returns a character string which is the lowercase equivalent of the supplied string.
  INTERFACE CharacterToLowercase
    MODULE PROCEDURE CharacterToLowercaseC
    MODULE PROCEDURE CharacterToLowercaseVS
  END INTERFACE CharacterToLowercase

  !>Returns a character string which is the uppercase equivalent of the supplied string.
  INTERFACE CharacterToUppercase
    MODULE PROCEDURE CharacterToUppercaseC
    MODULE PROCEDURE CharacterToUppercaseVS
  END INTERFACE CharacterToUppercase

  !>Returns .TRUE. if a supplied string is a valid abbreviation of a second supplied string.
  INTERFACE IsAbbreviation
    MODULE PROCEDURE IsAbbreviationCC
    MODULE PROCEDURE IsAbbreviationCVS
    MODULE PROCEDURE IsAbbreviationVSC
    MODULE PROCEDURE IsAbbreviationVSVS
  END INTERFACE IsAbbreviation
  
  !>Converts a list to its equivalent character string representation.
  INTERFACE ListToCharacter
    MODULE PROCEDURE ListToCharacterC
    MODULE PROCEDURE ListToCharacterIntg
    MODULE PROCEDURE ListToCharacterLIntg
    MODULE PROCEDURE ListToCharacterL
    MODULE PROCEDURE ListToCharacterSP
    MODULE PROCEDURE ListToCharacterDP
  END INTERFACE ListToCharacter

  !>Converts a number to its equivalent character string representation.
  INTERFACE NumberToCharacter
    MODULE PROCEDURE NumberToCharacterIntg
    MODULE PROCEDURE NumberToCharacterLIntg
    MODULE PROCEDURE NumberToCharacterSP
    MODULE PROCEDURE NumberToCharacterDP
  END INTERFACE NumberToCharacter

  !Provided to allow conversion to new code style
  INTERFACE NumberToVString
    MODULE PROCEDURE NumberToVStringIntg
    MODULE PROCEDURE NumberToVStringLIntg
    MODULE PROCEDURE NumberToVStringSP
    MODULE PROCEDURE NumberToVStringDP
  END INTERFACE NumberToVString

  !>Converts a string representation of a number to a double precision number.
  INTERFACE StringToDouble
    MODULE PROCEDURE StringToDoubleC
    MODULE PROCEDURE StringToDoubleVS
  END INTERFACE StringToDouble

  !>Converts a string representation of a number to an integer.
  INTERFACE StringToInteger
    MODULE PROCEDURE StringToIntegerC
    MODULE PROCEDURE StringToIntegerVS
  END INTERFACE StringToInteger

  !>Converts a string representation of a number to a long integer.
  INTERFACE StringToLongInteger
    MODULE PROCEDURE StringToLongIntegerC
    MODULE PROCEDURE StringToLongIntegerVS
  END INTERFACE StringToLongInteger

 !>Converts a string representation of a boolean value (TRUE or FALSE) to a logical.
  INTERFACE StringToLogical
    MODULE PROCEDURE StringToLogicalC
    MODULE PROCEDURE StringToLogicalVS
  END INTERFACE StringToLogical

   !>Converts a string representation of a number to a single precision number.
  INTERFACE StingToSingle
    MODULE PROCEDURE StringToSingleC
    MODULE PROCEDURE StringToSingleVS
  END INTERFACE StingToSingle

  !>Returns a varying string which is the lowercase equivalent of the supplied string.
  INTERFACE VStringToLowercase
    MODULE PROCEDURE VStringToLowercaseC
    MODULE PROCEDURE VStringToLowercaseVS
  END INTERFACE VStringToLowercase

  !>Returns a varying string which is the uppercase equivalent of the supplied string.
  INTERFACE VStringToUppercase
    MODULE PROCEDURE VStringToUppercaseC
    MODULE PROCEDURE VStringToUppercaseVS
  END INTERFACE VStringToUppercase

  PUBLIC CharacterToLowercase,CharacterToUppercase
  
  PUBLIC IsAbbreviation,IsDigit,IsLetter,IsWhitespace
  
  PUBLIC ListToCharacter
  
  PUBLIC LogicalToCharacter,LogicalToVString
  
  PUBLIC NumberToCharacter,NumberToVString

  PUBLIC StringToDouble,StringToInteger,StringToLongInteger,StringToLogical,StingToSingle

  PUBLIC VStringToLowercase,VStringToUppercase
  
CONTAINS
  
  !
  !================================================================================================================================
  !

  !>IsAbbreviation returns .TRUE. if the character string short is an abbreviation of the character string long. short must be at least minimumNumberOfCharacters long.
  PURE FUNCTION IsAbbreviationCC(short,long,minimumNumberOfCharacters)

    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: short !<The short form of the string
    CHARACTER(LEN=*), INTENT(IN) :: long !<The long form of the string
    INTEGER(INTG), INTENT(IN) :: minimumNumberOfCharacters !<The minimum number of characters to match
    !Function variable
    LOGICAL :: IsAbbreviationCC !<On exit, .TRUE. if the short string is an abbreviation
    !Local Variables
    INTEGER(INTG) :: characterIdx,numberOfCharacters
    CHARACTER(LEN=LEN(short)) :: upperShort
    CHARACTER(LEN=LEN(long)) :: upperLong
    
    IsAbbreviationCC=.FALSE.
    upperShort=CharacterToUppercase(short)
    upperLong=CharacterToUppercase(long)
    numberOfCharacters=MIN(LEN(long),LEN(short))
    DO characterIdx=minimumNumberOfCharacters,numberOfCharacters
      IF(upperShort==upperLong(:characterIdx)) THEN
          IsAbbreviationCC=.TRUE.
          EXIT
      ENDIF
    ENDDO !characterIdx

    RETURN
  END FUNCTION IsAbbreviationCC

  !
  !================================================================================================================================
  !

  !>IsAbbreviation returns .TRUE. if the character string short is an abbreviation of the varying string long. short must be at least minimumNumberOfCharacters long.
  PURE FUNCTION IsAbbreviationCVS(short,long,minimumNumberOfCharacters)

    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: short !<The short form of the string
    TYPE(VARYING_STRING), INTENT(IN) :: long !<The long form of the string
    INTEGER(INTG), INTENT(IN) :: minimumNumberOfCharacters !<The minimum number of characters to match
    !Function variable
    LOGICAL :: IsAbbreviationCVS !<On exit, .TRUE. if the short string is an abbreviation
    !Local Variables
    INTEGER(INTG) :: characterIdx,numberOfCharacters
    CHARACTER(LEN=LEN(short)) :: upperShort
    TYPE(VARYING_STRING) :: upperLong
    
    IsAbbreviationCVS=.FALSE.
    upperShort=CharacterToUppercase(short)
    upperLong=VStringToUppercase(long)
    numberOfCharacters=MIN(LEN(long),LEN(short))
    DO characterIdx=minimumNumberOfCharacters,numberOfCharacters
      IF(upperShort==EXTRACT(upperLong,1,characterIdx)) THEN
          IsAbbreviationCVS=.TRUE.
          EXIT
      ENDIF
    ENDDO !characterIdx

    RETURN
  END FUNCTION IsAbbreviationCVS

  !
  !================================================================================================================================
  !

  !>IsAbbreviation returns .TRUE. if the varying string short is an abbreviation of the character string long. short must be at least minimumNumberOfCharacters long.
  PURE FUNCTION IsAbbreviationVSC(short,long,minimumNumberOfCharacters)

    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: short !<The short form of the string
    CHARACTER(LEN=*), INTENT(IN) :: long !<The long form of the string
    INTEGER(INTG), INTENT(IN) :: minimumNumberOfCharacters !<The minimum number of characters to match
    !Function variable
    LOGICAL :: IsAbbreviationVSC !<On exit, .TRUE. if the short string is an abbreviation
    !Local Variables
    INTEGER(INTG) :: characterIdx,numberOfCharacters
    TYPE(VARYING_STRING) :: upperShort
    CHARACTER(LEN=LEN(long)) :: upperLong
    
    IsAbbreviationVSC=.FALSE.
    upperShort=VStringToUppercase(short)
    upperLong=CharacterToUppercase(long)
    numberOfCharacters=MIN(LEN(long),LEN(short))
    DO characterIdx=minimumNumberOfCharacters,numberOfCharacters
      IF(upperShort==upperLong(:characterIdx)) THEN
          IsAbbreviationVSC=.TRUE.
          EXIT
      ENDIF
    ENDDO !characterIdx

    RETURN
  END FUNCTION IsAbbreviationVSC

  !
  !================================================================================================================================
  !

  !>IsAbbreviation returns .TRUE. if the varying string short is an abbreviation of the varying string long. short must be at least minimumNumberOfCharacters long.
  PURE FUNCTION IsAbbreviationVSVS(short,long,minimumNumberOfCharacters)

    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: short !<The short form of the string
    TYPE(VARYING_STRING), INTENT(IN) :: long !<The long form of the string
    INTEGER(INTG), INTENT(IN) :: minimumNumberOfCharacters !<The minimum number of characters to match
    !Function variable
    LOGICAL :: IsAbbreviationVSVS !<On exit, .TRUE. if the short string is an abbreviation
    !Local Variables
    INTEGER(INTG) :: characterIdx,numberOfCharacters
    TYPE(VARYING_STRING) :: upperShort,upperLong
    
    IsAbbreviationVSVS=.FALSE.
    upperShort=VStringToUppercase(short)
    upperLong=VStringToUppercase(long)
    numberOfCharacters=MIN(LEN(long),LEN(short))
    DO characterIdx=minimumNumberOfCharacters,numberOfCharacters
      IF(upperShort==EXTRACT(upperLong,1,characterIdx)) THEN
          IsAbbreviationVSVS=.TRUE.
          EXIT
      ENDIF
    ENDDO !characterIdx

    RETURN
  END FUNCTION IsAbbreviationVSVS

  !
  !================================================================================================================================
  !

  !>IsDigit returns .TRUE. if the character charac is a digit character (i.e. 0..9)
  PURE FUNCTION IsDigit(charac)

    !Argument variables
    CHARACTER(LEN=1), INTENT(IN) :: charac !<The character to test if it is a digit
    !Function variable
    LOGICAL :: IsDigit !<On exit, .TRUE. if the character is a digit
    !Local Variables

    IsDigit=(ICHAR(charac)>=ICHAR("0").AND.ICHAR(charac)<=ICHAR("9"))

    RETURN
  END FUNCTION IsDigit

  !
  !================================================================================================================================
  !

  !>IsLetter returns .TRUE. if the character charac is a letter character (i.e. A..Z or a..z)
  PURE FUNCTION IsLetter(charac)

    !Argument variables
    CHARACTER(LEN=1), INTENT(IN) :: charac !<The character to test if it is a letter
    !Function variable
    LOGICAL :: IsLetter !<On exit, .TRUE. if the character is a letter
    !Local Variables

    IsLetter=((ICHAR(charac)>=ICHAR("A").AND.ICHAR(charac)<=ICHAR("Z")).OR.&
	    & (ICHAR(charac)>=ICHAR("a").AND.ICHAR(charac)<=ICHAR("z")))

    RETURN
  END FUNCTION IsLetter

  !
  !================================================================================================================================
  !

  !>Returns .TRUE. if the supplied character is a lowercase character.
  PURE FUNCTION IsLowercase(charac)

    !Argument variables
    CHARACTER(LEN=1), INTENT(IN) :: charac !<The character to test if it is lowercase
    !Function variable
    LOGICAL :: IsLowercase !<On exit, .TRUE. if the character is lowercase
    !Local Variables

    IF(LGE(charac,"a").AND.LLE(charac,"z")) THEN
      IsLowercase=.TRUE.
    ELSE
      IsLowercase=.FALSE.
    ENDIF

    RETURN
  END FUNCTION IsLowercase

  !
  !================================================================================================================================
  !

  !>Returns .TRUE. if the supplied character is an uppercase character.
  PURE FUNCTION IsUppercase(charac)

    !Argument variables
    CHARACTER(LEN=1), INTENT(IN) :: charac !<The character to test if it is uppercase
    !Function variable
    LOGICAL :: IsUppercase !<On exit, .TRUE. if the character is uppercase
    !Local Variables

    IF(LGE(charac,"A").AND.LLE(charac,"Z")) THEN
      IsUppercase=.TRUE.
    ELSE
      IsUppercase=.FALSE.
    ENDIF

    RETURN
  END FUNCTION IsUppercase

  !
  !================================================================================================================================
  !

  !>IsWhitespace returns .TRUE. if the character charac is a whitespace character (i.e. space, tabs, etc.)
  PURE FUNCTION IsWhitespace(charac)

    !Argument variables
    CHARACTER(LEN=1), INTENT(IN) :: charac !<The character to test if it is whitespace
    !Function variable
    LOGICAL :: IsWhitespace !<On exit, .TRUE. if the character is whitespace
    !Local Variables
    
    !!WARNING: Assumes ASCII encoding
    IsWhitespace=(charac==CHAR(32).OR.charac==CHAR(9))

    RETURN
  END FUNCTION IsWhitespace

  !
  !================================================================================================================================
  !

  !>Converts a character list to its equivalent character string representation as determined by the supplied format. If present, the optional argument listLengths is used for the lengths of each list elements length otherwise the trimmed length is used. NOTE: The format is ignored for this child FUNCTION.
  FUNCTION ListToCharacterC(numberInList,list,format,err,error,listLengths)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: numberInList !<The number of items in the list
    CHARACTER(LEN=*), INTENT(IN) :: list(numberInList) !<list(listItemIdx). The listItemIdx'th item in the list
    CHARACTER(LEN=*), INTENT(IN) :: format !<The format to use. Ignored for character lists.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    INTEGER(INTG), OPTIONAL, INTENT(IN) :: listLengths(numberInList) !<listLengths(listItemidx). Optional, The length of the listItemIdx'th list item.
    !Function variable
    CHARACTER(LEN=MAXSTRLEN) :: ListToCharacterC !<On exit, the character equivalent of the list
    !Local variables
    INTEGER(INTG) :: length,listItemIdx,position
    
    ENTERS("ListToCharacterC",err,error,*999)

    ListToCharacterC=""
    IF(numberInList>0) THEN
      IF(PRESENT(listLengths)) THEN
        length=listLengths(1)
        ListToCharacterC=list(1)(1:length)
        DO listItemIdx=2,numberInList
          IF(length+listLengths(listItemIdx)+1<=MAXSTRLEN) THEN
            ListToCharacterC=ListToCharacterC(1:length)//","//list(listItemIdx)(1:listLengths(listItemIdx))
            length=length+listLengths(listItemIdx)+1
          ELSE IF(length+5<=MAXSTRLEN) THEN
            ListToCharacterC=ListToCharacterC(1:length)//",...."
            EXIT
          ELSE
            position=INDEX(ListToCharacterC(1:MAXSTRLEN-4),",",.TRUE.)
            IF(position/=0) THEN
              ListToCharacterC=ListToCharacterC(1:position)//"...."
            ELSE
              ListToCharacterC=ListToCharacterC(1:MAXSTRLEN-5)//",...."
            ENDIF
            EXIT
          ENDIF
        ENDDO !listItemIdx
      ELSE
        ListToCharacterC=list(1)(1:LEN_TRIM(list(1)))
        DO listItemIdx=2,numberInList
          IF(LEN_TRIM(ListToCharacterC)+LEN_TRIM(list(listItemIdx))+1<=MAXSTRLEN) THEN
            ListToCharacterC=ListToCharacterC(1:LEN_TRIM(ListToCharacterC))//","//list(listItemIdx)(1:LEN_TRIM(list(listItemIdx)))
          ELSE IF(LEN_TRIM(ListToCharacterC)+5<=MAXSTRLEN) THEN
            ListToCharacterC=ListToCharacterC(1:LEN_TRIM(ListToCharacterC))//",...."
            EXIT
          ELSE
            position=INDEX(ListToCharacterC(1:MAXSTRLEN-4),",",.TRUE.)
            IF(position/=0) THEN
              ListToCharacterC=ListToCharacterC(1:position)//"...."
            ELSE
              ListToCharacterC=ListToCharacterC(1:MAXSTRLEN-5)//",...."
            ENDIF
            EXIT
          ENDIF
        ENDDO !listItemIdx
      ENDIF
    ENDIF

    EXITS("ListToCharacterC")
    RETURN
999 ERRORSEXITS("ListToCharacterC",err,error)
    RETURN
    
  END FUNCTION ListToCharacterC
  
  !
  !================================================================================================================================
  !
  
  !>Converts an integer list to its equivalent character string representation as determined by the supplied format. 
  FUNCTION ListToCharacterIntg(numberInList,list,format,err,error)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: numberInList !<The number of items in the list
    INTEGER(INTG), INTENT(IN) :: list(numberInList) !<list(listItemIdx). The listItemIdx'th item in the list
    CHARACTER(LEN=*), INTENT(IN) :: format !<The format to use for the conversion
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    CHARACTER(LEN=MAXSTRLEN) :: ListToCharacterIntg !<On exit, the character equivalent of the list
    !Local variables
    INTEGER(INTG) :: listItemIdx,position
    CHARACTER(LEN=MAXSTRLEN) :: listValue
    
    ENTERS("ListToCharacterIntg",err,error,*999)

    ListToCharacterIntg=""
    IF(numberInList>0) THEN
      ListToCharacterIntg=NumberToCharacterIntg(list(1),format,err,error)
      IF(err/=0) GOTO 999
      DO listItemIdx=2,numberInList
        listValue=NumberToCharacterIntg(list(listItemIdx),format,err,error)
        IF(err/=0) GOTO 999
        IF(LEN_TRIM(ListToCharacterIntg)+LEN_TRIM(listValue)+1<=MAXSTRLEN) THEN
          ListToCharacterIntg=ListToCharacterIntg(1:LEN_TRIM(ListToCharacterIntg))//","// &
            & listValue(1:LEN_TRIM(listValue))
        ELSE IF(LEN_TRIM(ListToCharacterIntg)+5<=MAXSTRLEN) THEN
          ListToCharacterIntg=ListToCharacterIntg(1:LEN_TRIM(ListToCharacterIntg))//",...."
          EXIT
        ELSE
          position=INDEX(ListToCharacterIntg(1:MAXSTRLEN-4),",",.TRUE.)
          IF(position/=0) THEN
            ListToCharacterIntg=ListToCharacterIntg(1:position)//"...."
          ELSE
            ListToCharacterIntg=ListToCharacterIntg(1:MAXSTRLEN-5)//",...."
          ENDIF
          EXIT
        ENDIF
      ENDDO !listItemIdx
    ENDIF

    EXITS("ListToCharacterIntg")
    RETURN
999 ERRORSEXITS("ListToCharacterIntg",err,error)
    RETURN
    
  END FUNCTION ListToCharacterIntg
  
  !
  !================================================================================================================================
  !

  !>Converts an long integer list to its equivalent character string representation as determined by the supplied format. 
  FUNCTION ListToCharacterLIntg(numberInList,list,format,err,error)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: numberInList !<The number of items in the list
    INTEGER(LINTG), INTENT(IN) :: list(numberInList) !<list(listItemIdx). The listItemIdx'th item in the list
    CHARACTER(LEN=*), INTENT(IN) :: format !<The format to use for the conversion
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    CHARACTER(LEN=MAXSTRLEN) :: ListToCharacterLIntg !<On exit, the character equivalent of the list
    !Local variables
    INTEGER(INTG) :: listItemIdx,position
    CHARACTER(LEN=MAXSTRLEN) :: listValue
    
    ENTERS("ListToCharacterLIntg",err,error,*999)

    ListToCharacterLIntg=""
    IF(numberInList>0) THEN
      ListToCharacterLIntg=NumberToCharacterLIntg(list(1),format,err,error)
      IF(err/=0) GOTO 999
      DO listItemIdx=2,numberInList
        listValue=NumberToCharacterLIntg(list(listItemIdx),format,err,error)
        IF(err/=0) GOTO 999
        IF(LEN_TRIM(ListToCharacterLIntg)+LEN_TRIM(listValue)+1<=MAXSTRLEN) THEN
          ListToCharacterLIntg=ListToCharacterLIntg(1:LEN_TRIM(ListToCharacterLIntg))//","// &
            & listValue(1:LEN_TRIM(listValue))
        ELSE IF(LEN_TRIM(ListToCharacterLIntg)+5<=MAXSTRLEN) THEN
          ListToCharacterLIntg=ListToCharacterLIntg(1:LEN_TRIM(ListToCharacterLIntg))//",...."
          EXIT
        ELSE
          position=INDEX(ListToCharacterLIntg(1:MAXSTRLEN-4),",",.TRUE.)
          IF(position/=0) THEN
            ListToCharacterLIntg=ListToCharacterLIntg(1:position)//"...."
          ELSE
            ListToCharacterLIntg=ListToCharacterLIntg(1:MAXSTRLEN-5)//",...."
          ENDIF
          EXIT
        ENDIF
      ENDDO !listItemIdx
    ENDIF

    EXITS("ListToCharacterLIntg")
    RETURN
999 ERRORSEXITS("ListToCharacterLIntg",err,error)
    RETURN
    
  END FUNCTION ListToCharacterLIntg
  
  !
  !================================================================================================================================
  !

  !>Converts a logical list to its equivalent character string representation as determined by the supplied format string. The format is ignored for this child FUNCTION.
  FUNCTION ListToCharacterL(numberInList,list,format,err,error)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: numberInList !<The number of items in the list
    LOGICAL, INTENT(IN) :: list(numberInList) !<list(listItemIdx). The listItemIdx'th item in the list
    CHARACTER(LEN=*), INTENT(IN) :: format !<The format to use. Ignored for logical lists.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    CHARACTER(LEN=MAXSTRLEN) :: ListToCharacterL !<On exit, the character equivalent of the list
    !Local variables
    INTEGER(INTG) :: listItemIdx,position
    CHARACTER(LEN=MAXSTRLEN) :: listValue
    
    ENTERS("ListToCharacterL",err,error,*999)

    ListToCharacterL=""
    IF(numberInList>0) THEN
      ListToCharacterL=LogicalToCharacter(list(1),err,error)
      IF(err/=0) GOTO 999
      DO listItemIdx=2,numberInList
        listValue=LogicalToCharacter(list(listItemIdx),err,error)
        IF(err/=0) GOTO 999
        IF(LEN_TRIM(ListToCharacterL)+LEN_TRIM(listValue)+1<=MAXSTRLEN) THEN
          ListToCharacterL=ListToCharacterL(1:LEN_TRIM(ListToCharacterL))//","//listValue(1:LEN_TRIM(listValue))
        ELSE IF(LEN_TRIM(ListToCharacterL)+5<=MAXSTRLEN) THEN
          ListToCharacterL=ListToCharacterL(1:LEN_TRIM(ListToCharacterL))//",...."
          EXIT
        ELSE
          position=INDEX(ListToCharacterL(1:MAXSTRLEN-4),",",.TRUE.)
          IF(position/=0) THEN
            ListToCharacterL=ListToCharacterL(1:position)//"...."
          ELSE
            ListToCharacterL=ListToCharacterL(1:MAXSTRLEN-5)//",...."
          ENDIF
          EXIT
        ENDIF
      ENDDO !listItemIdx
    ENDIF

    EXITS("ListToCharacterL")
    RETURN
999 ERRORSEXITS("ListToCharacterL",err,error)
    RETURN
    
  END FUNCTION ListToCharacterL
  
  !
  !================================================================================================================================
  !

  !>Converts a single precision list to its equivalent character string representation as determined by the supplied format string.
  FUNCTION ListToCharacterSP(numberInList,list,format,err,error)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: numberInList !<The number of items in the list
    REAL(SP), INTENT(IN) :: list(numberInList) !<list(listItemIdx). The listItemIdx'th item in the list 
    CHARACTER(LEN=*), INTENT(IN) :: format !<The format to use for the conversion
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    CHARACTER(LEN=MAXSTRLEN) :: ListToCharacterSP !<On exit, the character equivalent of the list
    !Local variables
    INTEGER(INTG) :: listItemIdx,position
    CHARACTER(LEN=MAXSTRLEN) :: listValue
    
    ENTERS("ListToCharacterSP",err,error,*999)

    ListToCharacterSP=""
    IF(numberInList>0) THEN
      ListToCharacterSP=NumberToCharacterSP(list(1),format,err,error)
      IF(err/=0) GOTO 999
      DO listItemIdx=2,numberInList
        listValue=NumberToCharacterSP(list(listItemIdx),format,err,error)
        IF(err/=0) GOTO 999
        IF(LEN_TRIM(ListToCharacterSP)+LEN_TRIM(listValue)+1<=MAXSTRLEN) THEN
          ListToCharacterSP=ListToCharacterSP(1:LEN_TRIM(ListToCharacterSP))//","//listValue(1:LEN_TRIM(listValue))
        ELSE IF(LEN_TRIM(ListToCharacterSP)+5<=MAXSTRLEN) THEN
          ListToCharacterSP=ListToCharacterSP(1:LEN_TRIM(ListToCharacterSP))//",...."
          EXIT
        ELSE
          position=INDEX(ListToCharacterSP(1:MAXSTRLEN-4),",",.TRUE.)
          IF(position/=0) THEN
            ListToCharacterSP=ListToCharacterSP(1:position)//"...."
          ELSE
            ListToCharacterSP=ListToCharacterSP(1:MAXSTRLEN-5)//",...."
          ENDIF
          EXIT
        ENDIF
      ENDDO !listItemIdx
    ENDIF

    EXITS("ListToCharacterSP")
    RETURN
999 ERRORSEXITS("ListToCharacterSP",err,error)
    RETURN
    
  END FUNCTION ListToCharacterSP
  
  !
  !================================================================================================================================
  !

  !>Converts a double precision list to its equivalent character string representation as determined by the supplied format string.
  FUNCTION ListToCharacterDP(numberInList,list,format,err,error)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: numberInList !<The number of items in the list
    REAL(DP), INTENT(IN) :: list(numberInList) !<list(listItemIdx). The listItemIdx'th item in the list
    CHARACTER(LEN=*), INTENT(IN) :: format !<The format to use for the conversion
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    CHARACTER(LEN=MAXSTRLEN) :: ListToCharacterDP !<On exit, the character equivalent of the list
    !Local variables
    INTEGER(INTG) :: listItemIdx,position
    CHARACTER(LEN=MAXSTRLEN) :: listValue
    
    ENTERS("ListToCharacterDP",err,error,*999)

    ListToCharacterDP=""
    IF(numberInList>0) THEN
      ListToCharacterDP=NumberToCharacterDP(list(1),format,err,error)
      IF(err/=0) GOTO 999
      DO listItemIdx=2,numberInList
        listValue=NumberToCharacterDP(list(listItemIdx),format,err,error)
        IF(err/=0) GOTO 999
        IF(LEN_TRIM(ListToCharacterDP)+LEN_TRIM(listValue)+1<=MAXSTRLEN) THEN
          ListToCharacterDP=ListToCharacterDP(1:LEN_TRIM(ListToCharacterDP))//","//listValue(1:LEN_TRIM(listValue))
        ELSE IF(LEN_TRIM(ListToCharacterDP)+5<=MAXSTRLEN) THEN
          ListToCharacterDP=ListToCharacterDP(1:LEN_TRIM(ListToCharacterDP))//",...."
          EXIT
        ELSE
          position=INDEX(ListToCharacterDP(1:MAXSTRLEN-4),",",.TRUE.)
          IF(position/=0) THEN
            ListToCharacterDP=ListToCharacterDP(1:position)//"...."
          ELSE
            ListToCharacterDP=ListToCharacterDP(1:MAXSTRLEN-5)//",...."
          ENDIF
          EXIT
        ENDIF
      ENDDO !listItemIdx
    ENDIF

    EXITS("ListToCharacterDP")
    RETURN
999 ERRORSEXITS("ListToCharacterDP",err,error)
    RETURN
    
  END FUNCTION ListToCharacterDP
  
  !
  !================================================================================================================================
  !
  
  !>Converts a logical value to either a "TRUE" or "FALSE" character string.
  FUNCTION LogicalToCharacter(logicalValue,err,error)
  
    !Argument variables
    LOGICAL, INTENT(IN) :: logicalValue !<The logical value to convert
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    CHARACTER(LEN=MAXSTRLEN) :: LogicalToCharacter !<On exit, the character equivalent value
    !Local variables
    
    !ENTERS("LogicalToCharacter",err,error,*999)

    IF(logicalValue) THEN
      LogicalToCharacter="TRUE"
    ELSE
      LogicalToCharacter="FALSE"
    ENDIF

    !EXITS("LogicalToCharacter")
    RETURN
!999 ERRORSEXITS("LogicalToCharacter",err,error)
999 ERRORS("LogicalToCharacter",err,error)
    RETURN
    
  END FUNCTION LogicalToCharacter
  
  !
  !================================================================================================================================
  !

  !>Converts a logical value to either a "TRUE" or "FALSE" varying string.
  FUNCTION LogicalToVString(logicalValue,err,error)
    
    !Argument variables
    LOGICAL, INTENT(IN) :: logicalValue !<The logical value to convert
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    TYPE(VARYING_STRING) :: LogicalToVString !<On exit, the varying string equivalent value
    !Local variables
    
    !ENTERS("LogicalToVString",err,error,*999)

    IF(logicalValue) THEN
      LogicalToVString="TRUE"
    ELSE
      LogicalToVString="FALSE"
    ENDIF

    !EXITS("LogicalToVString")
    RETURN
!999 ERRORSEXITS("LogicalToVString",err,error)
999 ERRORS("LogicalToVString",err,error)
    RETURN
    
  END FUNCTION LogicalToVString
  
  !
  !================================================================================================================================
  !

  !>Converts an integer number to its equivalent character string representation as determined by the supplied format. The format is of the form of a standard Fortran integer format e.g. "I3".
  FUNCTION NumberToCharacterIntg(number,format,err,error,adjust)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: number !<The number to convert
    CHARACTER(LEN=*), INTENT(IN) :: format !<The format to use in the conversion
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    LOGICAL, OPTIONAL, INTENT(IN) :: adjust !<Optional argument. If .TRUE. (default) the leading white space is stripped. 
    !Function variable
    CHARACTER(LEN=MAXSTRLEN) :: NumberToCharacterIntg !<On exit, the character equivalent of the number
    !Local variables
    CHARACTER(LEN=MAXSTRLEN) :: localFormat
    LOGICAL :: adjustLeft
    
    !ENTERS("NumberToCharacterIntg",err,error,*999)

    IF(PRESENT(adjust)) THEN
      adjustLeft=adjust
    ELSE
      adjustLeft=.TRUE.
    ENDIF
    
    IF(format(1:1)=="*") THEN
      localFormat="(I12)"
    ELSE
      localFormat="("//format(1:LEN_TRIM(format))//")"
    ENDIF
    WRITE(NumberToCharacterIntg,localFormat,ERR=999) number

    IF(adjustLeft) THEN
      !Trim leading blanks
      NumberToCharacterIntg=ADJUSTL(NumberToCharacterIntg)
    ENDIF

    !EXITS("NumberToCharacterIntg")
    RETURN
999 CALL FlagError("Error converting an integer to a character string.",err,error,*998)
!998 ERRORSEXITS("NumberToCharacterIntg",err,error)
998 ERRORS("NumberToCharacterIntg",err,error)
    RETURN
    
  END FUNCTION NumberToCharacterIntg
  
  !
  !================================================================================================================================
  !

  !>Converts a long integer number to its equivalent character string representation as determined by the supplied format. The format is of the form of a standard Fortran integer format e.g. "I3".
  FUNCTION NumberToCharacterLIntg(number,format,err,error,adjust)
  
    !Argument variables
    INTEGER(LINTG), INTENT(IN) :: number !<The number to convert
    CHARACTER(LEN=*), INTENT(IN) :: format !<The format to use in the conversion
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    LOGICAL, OPTIONAL, INTENT(IN) :: adjust !<Optional argument. If .TRUE. (default) the leading white space is stripped. 
    !Function variable
    CHARACTER(LEN=MAXSTRLEN) :: NumberToCharacterLIntg !<On exit, the character equivalent of the number
    !Local variables
    CHARACTER(LEN=MAXSTRLEN) :: localFormat
    LOGICAL :: adjustLeft
   
    !ENTERS("NumberToCharacterLIntg",err,error,*999)

    IF(PRESENT(adjust)) THEN
      adjustLeft=adjust
    ELSE
      adjustLeft=.TRUE.
    ENDIF

    IF(format(1:1)=="*") THEN
      localFormat="(I18)"
    ELSE
      localFormat="("//format(1:LEN_TRIM(format))//")"
    ENDIF
    WRITE(NumberToCharacterLIntg,localFormat,ERR=999) number

    IF(adjustLeft) THEN
      !Trim leading blanks
      NumberToCharacterLIntg=ADJUSTL(NumberToCharacterLIntg)
    ENDIF
      
    !EXITS("NumberToCharacterLIntg")
    RETURN
999 CALL FlagError("Error converting a long integer to a character string.",err,error,*998)
!998 ERRORSEXITS("NumberToCharacterLIntg",err,error)
998 ERRORS("NumberToCharacterLIntg",err,error)
    RETURN
    
  END FUNCTION NumberToCharacterLIntg
  
  !
  !================================================================================================================================
  !

  !>Converts a single precision number to its equivalent character string representation as determined by the supplied format string. NOTE: If format is an asterisk followed by a number between 1 and 32 the format will be chosen to maximise the number of significant digits, e.g., format="*8" will return a string of 8 characters representing the supplied number in either F8.? or E8.? format depending on its magnitude.
  FUNCTION NumberToCharacterSP(number,format,err,error,adjust)
  
    !Argument variables
    REAL(SP), INTENT(IN) :: number !<The number to convert
    CHARACTER(LEN=*), INTENT(IN) :: format !<The format to use in the conversion
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    LOGICAL, OPTIONAL, INTENT(IN) :: adjust !<Optional argument. If .TRUE. (default) the leading white space is stripped. 
    !Function variable
    CHARACTER(LEN=MAXSTRLEN) :: NumberToCharacterSP !<On exit, the character equivalent of the number
    !Local variables
    INTEGER(INTG) :: asteriskPosition,i0,i1,length
    CHARACTER(LEN=MAXSTRLEN) :: ci0,ci1
    CHARACTER(LEN=MAXSTRLEN) :: localFormat
    LOGICAL :: adjustLeft
    
    !ENTERS("NumberToCharacterSP",err,error,*999)

    IF(PRESENT(adjust)) THEN
      adjustLeft=adjust
    ELSE
      adjustLeft=.TRUE.
    ENDIF

    asteriskPosition=INDEX(format,"*")
    length=LEN_TRIM(format)
    IF(asteriskPosition==1.AND.length==1) THEN !Free format
      WRITE(NumberToCharacterSP,*,ERR=999) number      
    ELSE IF(asteriskPosition>0) THEN !Adjustable format
      ci0=FORMAT(asteriskPosition+1:LEN_TRIM(format))
      READ(ci0,'(BN,I2)') i0
      IF(i0<=MAXSTRLEN) THEN
        IF(number>=0.0_SP) THEN
          IF((number<10.0_SP**(i0-1)).AND.(number>=0.1_SP**(MIN(i0-2,5)))) THEN
            IF(number>1.0_SP) THEN
              i1=i0-2-FLOOR(LOG10(number))
              localFormat="(I2)"
              WRITE(ci1,localFormat) i1
              localFormat="(F"//ci0(1:LEN_TRIM(ci0))//"."//ci1(1:LEN_TRIM(ci1))//")"
              WRITE(NumberToCharacterSP,localFormat,ERR=999) number
            ELSE
              localFormat="(I2)"
              WRITE(ci1,localFormat) i0-2
              localFormat="(F"//ci0(1:LEN_TRIM(ci0))//"."//ci1(1:LEN_TRIM(ci1))//")"
              WRITE(NumberToCharacterSP,localFormat,ERR=999) number
            ENDIF
          ELSE
            localFormat="(I2)"
            WRITE(ci1,localFormat) i0-6
            localFormat="(E"//ci0(1:LEN_TRIM(ci0))//"."//ci1(1:LEN_TRIM(ci1))//")"
            WRITE(NumberToCharacterSP,localFormat,ERR=999) number
          ENDIF
        ELSE
          IF((-number<10.0_SP**(i0-2)).AND.(-number>=0.01_SP**(MIN(i0-2,5)))) THEN
            IF(-number>=1.0_SP) THEN
              i1=i0-3-FLOOR(LOG10(number))
              localFormat="(I2)"
              WRITE(ci1,'(I2)') i1
              localFormat="(F"//ci0(1:LEN_TRIM(ci0))//"."//ci1(1:LEN_TRIM(ci1))//")"
              WRITE(NumberToCharacterSP,localFormat,ERR=999) number
            ELSE
              localFormat="(I2)"
              WRITE(ci1,localFormat) i0-2
              localFormat="(F"//ci0(1:LEN_TRIM(ci0))//"."//ci1(1:LEN_TRIM(ci1))//")"
              WRITE(NumberToCharacterSP,localFormat,ERR=999) number
            ENDIF
          ELSE
            localFormat="(I2)"
            WRITE(ci1,localFormat) i0-6
            localFormat="(E"//ci0(1:LEN_TRIM(ci0))//"."//ci1(1:LEN_TRIM(ci1))//")"
            WRITE(NumberToCharacterSP,localFormat,ERR=999) number
          ENDIF
        ENDIF
      ELSE
        CALL FlagError("Invalid format.",err,error,*999)
        GOTO 999
      ENDIF
    ELSE
      localFormat='('//format(1:LEN_TRIM(format))//')'
      WRITE(NumberToCharacterSP,localFormat,ERR=999) number
    ENDIF

    !Add an extra zero if required
    IF(NumberToCharacterSP(LEN_TRIM(NumberToCharacterSP):LEN_TRIM(NumberToCharacterSP))==".") &
      & NumberToCharacterSP=NumberToCharacterSP(1:LEN_TRIM(NumberToCharacterSP))//"0"

    IF(adjustLeft) THEN
      !Trim leading blanks
      NumberToCharacterSP=ADJUSTL(NumberToCharacterSP)
    ENDIF

    !EXITS("NumberToCharacterSP")
    RETURN
999 CALL FlagError("Error converting a single precision number to a character string.",err,error,*998)
!998 ERRORSEXITS("NumberToCharacterSP",err,error)
998 ERRORS("NumberToCharacterSP",err,error)
    RETURN
    
  END FUNCTION NumberToCharacterSP
  
  !
  !================================================================================================================================
  !

  !>Converts a double precision number to its equivalent character string representation as determined by the supplied format string. Note If format is an asterisk followed by a number between 1 and 32 the format will be chosen to maximise the number of significant digits, e.g., format="*8" will return a string of 8 characters representing the supplied number in either F8.? or E8.? format depending on its magnitude.
  FUNCTION NumberToCharacterDP(number,format,err,error,adjust)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: number !<The number to convert
    CHARACTER(LEN=*), INTENT(IN) :: format !<The format to use in the conversion
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    LOGICAL, OPTIONAL, INTENT(IN) :: adjust !<Optional argument. If .TRUE. (default) the leading white space is stripped. 
    !Function variable
    CHARACTER(LEN=MAXSTRLEN) :: NumberToCharacterDP !<On exit, the character equivalent of the number
    !Local variables
    INTEGER(INTG) :: asteriskPosition,i0,i1,length
    CHARACTER(LEN=2) :: ci0,ci1
    CHARACTER(LEN=MAXSTRLEN) :: localFormat
    LOGICAL :: adjustLeft
    
    !ENTERS("NumberToCharacterDP",err,error,*999)

    IF(PRESENT(adjust)) THEN
      adjustLeft=adjust
    ELSE
      adjustLeft=.TRUE.
    ENDIF
    
    asteriskPosition=INDEX(format,"*")
    length=LEN_TRIM(format)
    IF(asteriskPosition==1.AND.length==1) THEN !Free format
      WRITE(NumberToCharacterDP,*,ERR=999) number      
    ELSE IF(asteriskPosition>0) THEN !Adjustable format
      ci0=FORMAT(asteriskPosition+1:LEN_TRIM(format))
      READ(ci0,'(BN,I2)') i0
      IF(i0<=MAXSTRLEN) THEN
        IF(number>=0.0_DP) THEN
          IF((number<10.0_DP**(i0-1)).AND.(number>=0.1_DP**(MIN(i0-2,5)))) THEN
            IF(number>1.0_DP) THEN
              i1=i0-2-FLOOR(LOG10(number))
              localFormat="(I2)"
              WRITE(ci1,localFormat) i1
              localFormat="(F"//ci0(1:LEN_TRIM(ci0))//"."//ci1(1:LEN_TRIM(ci1))//")"
              WRITE(NumberToCharacterDP,localFormat,ERR=999) number
            ELSE
              localFormat="(I2)"
              WRITE(ci1,localFormat) i0-2
              localFormat="(F"//ci0(1:LEN_TRIM(ci0))//"."//ci1(1:LEN_TRIM(ci1))//")"
              WRITE(NumberToCharacterDP,localFormat,ERR=999) number
            ENDIF
          ELSE
            localFormat="(I2)"
            WRITE(ci1,localFormat) i0-6
            localFormat="(E"//ci0(1:LEN_TRIM(ci0))//"."//ci1(1:LEN_TRIM(ci1))//")"
            WRITE(NumberToCharacterDP,localFormat,ERR=999) number
          ENDIF
        ELSE
          IF((-number<10.0_DP**(i0-2)).AND.(-number>=0.01_DP**(MIN(i0-2,5)))) THEN
            IF(-number>=1.0_DP) THEN
              i1=i0-3-FLOOR(LOG10(number))
              localFormat="(I2)"
              WRITE(ci1,'(I2)') i1
              localFormat="(F"//ci0(1:LEN_TRIM(ci0))//"."//ci1(1:LEN_TRIM(ci1))//")"
              WRITE(NumberToCharacterDP,localFormat,ERR=999) number
            ELSE
              localFormat="(I2)"
              WRITE(ci1,localFormat) i0-2
              localFormat="(F"//ci0(1:LEN_TRIM(ci0))//"."//ci1(1:LEN_TRIM(ci1))//")"
              WRITE(NumberToCharacterDP,localFormat,ERR=999) number
            ENDIF
          ELSE
            localFormat="(I2)"
            WRITE(ci1,localFormat) i0-6
            localFormat="(E"//ci0(1:LEN_TRIM(ci0))//"."//ci1(1:LEN_TRIM(ci1))//")"
            WRITE(NumberToCharacterDP,localFormat,ERR=999) number
          ENDIF
        ENDIF
      ELSE
        CALL FlagError("Invalid format.",err,error,*999)
      ENDIF
    ELSE
      localFormat='('//format(1:LEN_TRIM(format))//')'
      WRITE(NumberToCharacterDP,localFormat,ERR=999) number
    ENDIF

    !Add an extra zero if required
    IF(NumberToCharacterDP(LEN_TRIM(NumberToCharacterDP):LEN_TRIM(NumberToCharacterDP))==".") &
      & NumberToCharacterDP=NumberToCharacterDP(1:LEN_TRIM(NumberToCharacterDP))//"0"

    IF(adjustLeft) THEN
      !Trim leading blanks
      NumberToCharacterDP=ADJUSTL(NumberToCharacterDP)
    ENDIF

    !EXITS("NumberToCharacterDP")
    RETURN
999 CALL FlagError("Error converting double precision number to a character string.",err,error,*998)
!998 ERRORSEXITS("NumberToCharacterDP",err,error)
998 ERRORS("NumberToCharacterDP",err,error)
    RETURN
    
  END FUNCTION NumberToCharacterDP
  
  !
  !================================================================================================================================
  !

  !>Converts an integer number to its equivalent varying string representation as determined by the supplied format. The format is of the form of a standard Fortran integer format e.g. "I3".
  FUNCTION NumberToVStringIntg(number,format,err,error,adjust)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: number !<The number to convert
    CHARACTER(LEN=*), INTENT(IN) :: format !<The format to use in the conversion
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    LOGICAL, OPTIONAL, INTENT(IN) :: adjust !<Optional argument. If .TRUE. (default) the leading white space is stripped. 
    !Function variable
    TYPE(VARYING_STRING) :: NumberToVStringIntg !<On exit, the varying string equivalent of the number
    !Local variables
    LOGICAL :: adjustLeft
    CHARACTER(LEN=MAXSTRLEN) :: localFormat,localString
    
    !ENTERS("NumberToVStringIntg",err,error,*999)

!!TODO: remove dependance on localString
    
    NumberToVStringIntg=""
    
    IF(PRESENT(adjust)) THEN
      adjustLeft=adjust
    ELSE
      adjustLeft=.TRUE.
    ENDIF
    
    IF(format(1:1)=="*") THEN
      localFormat="(I12)"
    ELSE
      localFormat="("//format(1:LEN_TRIM(format))//")"
    ENDIF
    WRITE(localString,localFormat,ERR=999) number

    IF(adjustLeft) THEN
      !Trim leading blanks
      NumberToVStringIntg=ADJUSTL(localString(1:LEN_TRIM(localString)))
    ELSE
      NumberToVStringIntg=localString(1:LEN_TRIM(localString))
    ENDIF

    !EXITS("NumberToVStringIntg")
    RETURN
999 CALL FlagError("Error converting an integer to a varying string.",err,error,*998)
998 ERRORS("NumberToVStringIntg",err,error)
    !EXITS("NumberToVStringIntg")
    RETURN
    
  END FUNCTION NumberToVStringIntg
  
  !
  !================================================================================================================================
  !

  !>Converts a long integer number to its equivalent varying string representation as determined by the supplied format. The format is of the form of a standard Fortran integer format e.g., "I3".
  FUNCTION NumberToVStringLIntg(number,format,err,error,adjust)
  
    !Argument variables
    INTEGER(LINTG), INTENT(IN) :: number !<The number to convert
    CHARACTER(LEN=*), INTENT(IN) :: format !<The format to use in the conversion
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    LOGICAL, OPTIONAL, INTENT(IN) :: adjust !<Optional argument. If .TRUE. (default) the leading white space is stripped. 
    !Function variable
    TYPE(VARYING_STRING) :: NumberToVStringLIntg !<On exit, the varying string equivalent of the number
    !Local variables
    CHARACTER(LEN=MAXSTRLEN) :: localFormat,localString
    LOGICAL :: adjustLeft   

    !ENTERS("NumberToVStringLIntg",err,error,*999)

!!TODO: remove dependance on localString
    
    NumberToVStringLIntg=""
    
    IF(PRESENT(adjust)) THEN
      adjustLeft=adjust
    ELSE
      adjustLeft=.TRUE.
    ENDIF
    
    IF(format(1:1)=="*") THEN
      localFormat="(I18)"
    ELSE
      localFormat="("//format(1:LEN_TRIM(format))//")"
    ENDIF
    WRITE(localString,localFormat,ERR=999) number

    IF(adjustLeft) THEN
      !Trim leading blanks
      NumberToVStringLIntg=ADJUSTL(localString(1:LEN_TRIM(localString)))
    ELSE
      NumberToVStringLIntg=localString(1:LEN_TRIM(localString))
    ENDIF

    !EXITS("NumberToVStringLIntg")
    RETURN
999 CALL FlagError("Error converting a long integer to a varying string.",err,error,*998)
998 ERRORS("NumberToVStringLIntg",err,error)
    !EXITS("NumberToVStringLIntg")
    RETURN
    
  END FUNCTION NumberToVStringLIntg
  
  !
  !================================================================================================================================
  !

  !>Converts a single precision number to its equivalent varying string representation as determined by the supplied format string. Note If format is an asterisk followed by a number between 1 and 32 the format will be chosen to maximise the number of significant digits, e.g., format="*8" will return a string of 8 characters representing the supplied number in either F8.? or E8.? format depending on its magnitude.
  FUNCTION NumberToVStringSP(number,format,err,error,adjust)
  
    !Argument variables
    REAL(SP), INTENT(IN) :: number !<The number to convert
    CHARACTER(LEN=*), INTENT(IN) :: format !<The format to use in the conversion
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    LOGICAL, OPTIONAL, INTENT(IN) :: adjust !<Optional argument. If .TRUE. (default) the leading white space is stripped. 
    !Function variable
    TYPE(VARYING_STRING) :: NumberToVStringSP !<On exit, the varying string equivalent of the number
    !Local variables
    INTEGER(INTG) :: asteriskPosition,i0,i1,length
    CHARACTER(LEN=MAXSTRLEN) :: ci0,ci1
    CHARACTER(LEN=MAXSTRLEN) :: localFormat,localString
    LOGICAL :: adjustLeft
    
    !ENTERS("NumberToVStringSP",err,error,*999)

!!TODO: remove dependance on localString
    
    NumberToVStringSP=""    

    IF(PRESENT(adjust)) THEN
      adjustLeft=adjust
    ELSE
      adjustLeft=.TRUE.
    ENDIF

    asteriskPosition=INDEX(format,"*")
    length=LEN_TRIM(format)
    IF(asteriskPosition==1.AND.length==1) THEN !Free format
      WRITE(localString,*,ERR=999) number      
    ELSE IF(asteriskPosition>0) THEN !Adjustable format
      ci0=FORMAT(asteriskPosition+1:LEN_TRIM(format))
      READ(ci0,'(BN,I2)') i0
      IF(i0<=MAXSTRLEN) THEN
        IF(number>=0.0_SP) THEN
          IF((number<10.0_SP**(i0-1)).AND.(number>=0.1_SP**(MIN(i0-2,5)))) THEN
            IF(number>1.0_SP) THEN
              i1=i0-2-FLOOR(LOG10(number))
              localFormat="(I2)"
              WRITE(ci1,localFormat) i1
              localFormat="(F"//ci0(1:LEN_TRIM(ci0))//"."//ci1(1:LEN_TRIM(ci1))//")"
              WRITE(localString,localFormat,ERR=999) number
            ELSE
              localFormat="(I2)"
              WRITE(ci1,localFormat) i0-2
              localFormat="(F"//ci0(1:LEN_TRIM(ci0))//"."//ci1(1:LEN_TRIM(ci1))//")"
              WRITE(localString,localFormat,ERR=999) number
            ENDIF
          ELSE
            localFormat="(I2)"
            WRITE(ci1,localFormat) i0-6
            localFormat="(E"//ci0(1:LEN_TRIM(ci0))//"."//ci1(1:LEN_TRIM(ci1))//")"
            WRITE(localString,localFormat,ERR=999) number
          ENDIF
        ELSE
          IF((-number<10.0_SP**(i0-2)).AND.(-number>=0.01_SP**(MIN(i0-2,5)))) THEN
            IF(-number>=1.0_SP) THEN
              i1=i0-3-FLOOR(LOG10(number))
              localFormat="(I2)"
              WRITE(ci1,'(I2)') i1
              localFormat="(F"//ci0(1:LEN_TRIM(ci0))//"."//ci1(1:LEN_TRIM(ci1))//")"
              WRITE(localString,localFormat,ERR=999) number
            ELSE
              localFormat="(I2)"
              WRITE(ci1,localFormat) i0-2
              localFormat="(F"//ci0(1:LEN_TRIM(ci0))//"."//ci1(1:LEN_TRIM(ci1))//")"
              WRITE(localString,localFormat,ERR=999) number
            ENDIF
          ELSE
            localFormat="(I2)"
            WRITE(ci1,localFormat) i0-6
            localFormat="(E"//ci0(1:LEN_TRIM(ci0))//"."//ci1(1:LEN_TRIM(ci1))//")"
            WRITE(localString,localFormat,ERR=999) number
          ENDIF
        ENDIF
      ELSE
        CALL FlagError("Invalid format.",err,error,*999)
        GOTO 999
      ENDIF
    ELSE
      localFormat='('//format(1:LEN_TRIM(format))//')'
      WRITE(localString,localFormat,ERR=999) number
    ENDIF

    !Add an extra zero if required
    IF(localString(LEN_TRIM(localString):LEN_TRIM(localString))==".") localString=localString(1:LEN_TRIM(localString))//"0"
    
    IF(adjustLeft) THEN
      !Trim leading blanks
      NumberToVStringSP=ADJUSTL(localString(1:LEN_TRIM(localString)))
    ELSE
      NumberToVStringSP=localString(1:LEN_TRIM(localString))
    ENDIF

    !EXITS("NumberToVStringSP")
    RETURN
999 CALL FlagError("Error converting a single precision number to a varying string.",err,error,*998)
998 ERRORS("NumberToVStringSP",err,error)
    !EXITS("NumberToVStringSP")
    RETURN
    
  END FUNCTION NumberToVStringSP
  
  !
  !================================================================================================================================
  !

  !>Converts a double precision number to its equivalent varying string representation as determined by the supplied format string. Note If format is an asterisk followed by a number between 1 and 32 the format will be chosen to maximise the number of significant digits, e.g., format="*8" will return a string of 8 characters representing the supplied number in either F8.? or E8.? format depending on its magnitude.
  FUNCTION NumberToVStringDP(number,format,err,error,adjust)
      
    !Argument variables
    REAL(DP), INTENT(IN) :: number !<The number to convert
    CHARACTER(LEN=*), INTENT(IN) :: format !<The format to use in the conversion
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    LOGICAL, OPTIONAL, INTENT(IN) :: adjust !<Optional argument. If .TRUE. (default) the leading white space is stripped. 
    !Function variable
    TYPE(VARYING_STRING) :: NumberToVStringDP !<On exit, the varying string equivalent of the number
    !Local variables
    INTEGER(INTG) :: asteriskPosition,i0,i1,length
    CHARACTER(LEN=2) :: ci0,ci1
    CHARACTER(LEN=MAXSTRLEN) :: localFormat,localString
    LOGICAL :: adjustLeft
     
    !ENTERS("NumberToVStringDP",err,error,*999)

!!TODO: remove dependance on localString
    
    NumberToVStringDP=""    

    IF(PRESENT(adjust)) THEN
      adjustLeft=adjust
    ELSE
      adjustLeft=.TRUE.
    ENDIF
    
    asteriskPosition=INDEX(format,"*")
    length=LEN_TRIM(format)
    IF(asteriskPosition==1.AND.length==1) THEN !Free format
      WRITE(localString,*,ERR=999) number      
    ELSE IF(asteriskPosition>0) THEN !Adjustable format
      ci0=FORMAT(asteriskPosition+1:LEN_TRIM(format))
      READ(ci0,'(BN,I2)') i0
      IF(i0<=MAXSTRLEN) THEN
        IF(number>=0.0_DP) THEN
          IF((number<10.0_DP**(i0-1)).AND.(number>=0.1_DP**(MIN(i0-2,5)))) THEN
            IF(number>1.0_DP) THEN
              i1=i0-2-FLOOR(LOG10(number))
              localFormat="(I2)"
              WRITE(ci1,localFormat) i1
              localFormat="(F"//ci0(1:LEN_TRIM(ci0))//"."//ci1(1:LEN_TRIM(ci1))//")"
              WRITE(localString,localFormat,ERR=999) number
            ELSE
              localFormat="(I2)"
              WRITE(ci1,localFormat) i0-2
              localFormat="(F"//ci0(1:LEN_TRIM(ci0))//"."//ci1(1:LEN_TRIM(ci1))//")"
              WRITE(localString,localFormat,ERR=999) number
            ENDIF
          ELSE
            localFormat="(I2)"
            WRITE(ci1,localFormat) i0-6
            localFormat="(E"//ci0(1:LEN_TRIM(ci0))//"."//ci1(1:LEN_TRIM(ci1))//")"
            WRITE(localString,localFormat,ERR=999) number
          ENDIF
        ELSE
          IF((-number<10.0_DP**(i0-2)).AND.(-number>=0.01_DP**(MIN(i0-2,5)))) THEN
            IF(-number>=1.0_DP) THEN
              i1=i0-3-FLOOR(LOG10(number))
              localFormat="(I2)"
              WRITE(ci1,'(I2)') i1
              localFormat="(F"//ci0(1:LEN_TRIM(ci0))//"."//ci1(1:LEN_TRIM(ci1))//")"
              WRITE(localString,localFormat,ERR=999) number
            ELSE
              localFormat="(I2)"
              WRITE(ci1,localFormat) i0-2
              localFormat="(F"//ci0(1:LEN_TRIM(ci0))//"."//ci1(1:LEN_TRIM(ci1))//")"
              WRITE(localString,localFormat,ERR=999) number
            ENDIF
          ELSE
            localFormat="(I2)"
            WRITE(ci1,localFormat) i0-6
            localFormat="(E"//ci0(1:LEN_TRIM(ci0))//"."//ci1(1:LEN_TRIM(ci1))//")"
            WRITE(localString,localFormat,ERR=999) number
          ENDIF
        ENDIF
      ELSE
        CALL FlagError("Invalid format.",err,error,*999)
      ENDIF
    ELSE
      localFormat='('//format(1:LEN_TRIM(format))//')'
      WRITE(localString,localFormat,ERR=999) number
    ENDIF

    !Add an extra zero if required
    IF(localString(LEN_TRIM(localString):LEN_TRIM(localString))==".") localString=localString(1:LEN_TRIM(localString))//"0"

    IF(adjustLeft) THEN
!!Do you really want to do this???
      !Trim leading blanks
      !NumberToVStringDP=ADJUSTL(localString(1:LEN_TRIM(localString)))
      NumberToVStringDP=localString(1:LEN_TRIM(localString))
    ELSE
      NumberToVStringDP=localString(1:LEN_TRIM(localString))
    ENDIF

    !EXITS("NumberToVStringDP")
    RETURN
999 CALL FlagError("Error converting double precision number to a varying string.",err,error,*998)
998 ERRORS("NumberToVStringDP",err,error)
    !EXITS("NumberToVStringDP")
    RETURN
    
  END FUNCTION NumberToVStringDP
  
  !
  !================================================================================================================================
  !

  !>Converts a character string representation of a number to a double precision number.
  FUNCTION StringToDoubleC(string,err,error)
  
    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: string !<The string to convert
    INTEGER(INTG), INTENT(OUT) :: err !<The error code !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string !<The error string
    !Function variable
    REAL(DP) :: StringToDoubleC !<On exit, the double precision equivalent of the string
    !Local variables
    
    !ENTERS("StringToDoubleC",err,error,*999)

    READ(string,*,IOSTAT=err,ERR=999) StringToDoubleC

    !EXITS("StringToDoubleC")
    RETURN
999 CALL FlagError("Cannot convert '"//STRING(1:LEN_TRIM(string))//"' to a double real.",err,error,*998)
!998 ERRORSEXITS("StringToDoubleC",err,error)
998 ERRORS("StringToDoubleC",err,error)
    RETURN
    
  END FUNCTION StringToDoubleC
  
  !
  !================================================================================================================================
  !

  !>Converts a varying string representation of a number to a double precision number.
  FUNCTION StringToDoubleVS(string,err,error)
  
    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: string !<The string to convert
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    REAL(DP) :: StringToDoubleVS !<On exit, the double precision equivalent of the string
    !Local variables
    CHARACTER(LEN=MAXSTRLEN) :: localString
    
    !ENTERS("StringToDoubleVS",err,error,*999)

!!TODO: remove dependance on localString

    localString=CHAR(string)
    READ(localString,*,IOSTAT=err,ERR=999) StringToDoubleVS

    !EXITS("StringToDoubleVS")
    RETURN
999 CALL FlagError("Cannot convert '"//CHAR(string)//"' to a double real.",err,error,*998)
!998 ERRORSEXITS("StringToDoubleVS",err,error)
998 ERRORS("StringToDoubleVS",err,error)
    RETURN
    
  END FUNCTION StringToDoubleVS
  
  !
  !================================================================================================================================
  !

  !>Converts a character string representation of a number to an integer.
  FUNCTION StringToIntegerC(string,err,error)
  
    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: string !<The string to convert
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    INTEGER(INTG) :: StringToIntegerC !<On exit, the integer equivalent of the string
    !Local variables
    
    !ENTERS("StringToIntegerC",err,error,*999)

    READ(string,*,IOSTAT=err,ERR=999) StringToIntegerC

    !EXITS("StringToIntegerC")
    RETURN
999 CALL FlagError("Cannot convert '"//STRING(1:LEN_TRIM(string))//"' to an integer.",err,error,*998)
!998 ERRORSEXITS("StringToIntegerC",err,error)
998 ERRORS("StringToIntegerC",err,error)
    RETURN
    
  END FUNCTION StringToIntegerC
  
  !
  !================================================================================================================================
  !

  !>Converts a varying string representation of a number to an integer.
  FUNCTION StringToIntegerVS(string,err,error)
  
    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: string !<The string to convert
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    INTEGER(INTG) :: StringToIntegerVS !<On exit, the integer equivalent of the string
    !Local variables
    CHARACTER(LEN=MAXSTRLEN) :: localString
    
    !ENTERS("StringToIntegerVS",err,error,*999)

!!TODO: remove dependance on localString

    localString=CHAR(string)
    READ(localString,*,IOSTAT=err,ERR=999) StringToIntegerVS

    !EXITS("StringToIntegerVS")
    RETURN
999 CALL FlagError("Cannot convert '"//CHAR(string)//"' to an integer.",err,error,*998)
!998 ERRORSEXITS("StringToIntegerVS",err,error)
998 ERRORS("StringToIntegerVS",err,error)
    RETURN
    
  END FUNCTION StringToIntegerVS
  
  !
  !================================================================================================================================
  !

  !>Converts a character string representation of a number to a long integer.
  FUNCTION StringToLongIntegerC(string,err,error)
  
    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: string !<The string to convert
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    INTEGER(LINTG) :: StringToLongIntegerC !<On exit, the long integer equivalent of the string
    !Local variables
    
    !ENTERS("StringToLongIntegerC",err,error,*999)

    READ(string,*,IOSTAT=err,ERR=999) StringToLongIntegerC

    !EXITS("StringToLongIntegerC")
    RETURN
999 CALL FlagError("Cannot convert '"//STRING(1:LEN_TRIM(string))//"' to a long integer.",err,error,*998)
!998 ERRORSEXITS("StringToLongIntegerC",err,error)
998 ERRORS("StringToLongIntegerC",err,error)
    RETURN
    
  END FUNCTION StringToLongIntegerC
  
  !
  !================================================================================================================================
  !

  !>Converts a varying string representation of a number to a long integer.
  FUNCTION StringToLongIntegerVS(string,err,error)
  
    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: string !<The string to convert
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    INTEGER(LINTG) :: StringToLongIntegerVS !<On exit, the long integer equivalent of the string
    !Local variables
    CHARACTER(LEN=MAXSTRLEN) :: localString
    
    !ENTERS("StringToLongIntegerVS",err,error,*999)

!!TODO: remove dependance on localString

    localString=CHAR(string)
    READ(localString,*,IOSTAT=err,ERR=999) StringToLongIntegerVS

    !EXITS("StringToLongIntegerVS")
    RETURN
999 CALL FlagError("Cannot convert '"//CHAR(string)//"' to a long integer.",err,error,*998)
!998 ERRORSEXITS("StringToLongIntegerVS",err,error)
998 ERRORS("StringToLongIntegerVS",err,error)
    RETURN
    
  END FUNCTION StringToLongIntegerVS
  
  !
  !================================================================================================================================
  !

  !>Converts a character string representation of a boolean (TRUE or FALSE) to a logical.
  FUNCTION StringToLogicalC(string,err,error)
  
    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: string !<The string to convert
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    LOGICAL :: StringToLogicalC !<On exit, the logical equivalent of the string
    !Local variables
    
    !ENTERS("StringToLogicalC",err,error,*999)

    READ(string,*,IOSTAT=err,ERR=999) StringToLogicalC

    !EXITS("StringToLogicalC")
    RETURN
999 CALL FlagError("Cannot convert '"//STRING(1:LEN_TRIM(string))//"' to a logical.",err,error,*998)
!998 ERRORSEXITS("StringToLogicalC",err,error)
998 ERRORS("StringToLogicalC",err,error)
    RETURN
    
  END FUNCTION StringToLogicalC
  
  !
  !================================================================================================================================
  !

  !>Converts a varying string representation of a boolean (TRUE or FALSE) to a logical.
  FUNCTION StringToLogicalVS(string,err,error)
  
    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: string !<The string to convert
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    LOGICAL :: StringToLogicalVS !<On exit, the logical equivalent of the string
    !Local variables
    CHARACTER(LEN=MAXSTRLEN) :: localString
    
    !ENTERS("StringToLogicalVS",err,error,*999)

    localString=CHAR(string)
    READ(localString,*,IOSTAT=err,ERR=999) StringToLogicalVS

    !EXITS("StringToLogicalVS")
    RETURN
999 CALL FlagError("Cannot convert '"//CHAR(string)//"' to a logical.",err,error,*998)
!998 ERRORSEXITS("StringToLogicalVS",err,error)
998 ERRORS("StringToLogicalVS",err,error)
    RETURN
    
  END FUNCTION StringToLogicalVS
  
  !
  !================================================================================================================================
  !

  !>Converts a character string representation of a number to a single precision number.
  FUNCTION StringToSingleC(string,err,error)
  
    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: string !<The string to convert
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    REAL(SP) :: StringToSingleC !<On exit, the single precision equivalent of the string
    !Local variables
    
    !ENTERS("StringToSingleC",err,error,*999)

    READ(string,*,IOSTAT=err,ERR=999) StringToSingleC
    
    !EXITS("StringToSingleC")
    RETURN
999 CALL FlagError("Cannot convert '"//STRING(1:LEN_TRIM(string))//"' to a single real.",err,error,*998)
!998 ERRORSEXITS("StringToSingleC",err,error)
998 ERRORS("StringToSingleC",err,error)
    RETURN
    
  END FUNCTION StringToSingleC
    
  !
  !================================================================================================================================
  !

  !>Converts a varying string representation of a number to a single precision number.
  FUNCTION StringToSingleVS(string,err,error)
  
    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: string !<The string to convert
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    REAL(SP) :: StringToSingleVS !<On exit, the single precision equivalent of the string
    !Local variables
    CHARACTER(LEN=MAXSTRLEN) :: localString
    
    !ENTERS("StringToSingleVS",err,error,*999)

!!TODO: remove dependance on localString

    localString=CHAR(string)
    READ(localString,*,IOSTAT=err,ERR=999) StringToSingleVS
    
    !EXITS("StringToSingleVS")
    RETURN
999 CALL FlagError("Cannot convert '"//CHAR(string)//"' to a single real.",err,error,*998)
!998 ERRORSEXITS("StringToSingleVS",err,error)
998 ERRORS("StringToSingleVS",err,error)
    RETURN
    
  END FUNCTION StringToSingleVS
    
  !
  !================================================================================================================================
  !
    
  PURE FUNCTION CharacterToLowercaseC(string)

    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: string !<The string to convert to lowercase
    !Function variable
    CHARACTER(LEN=LEN(string)) :: CharacterToLowercaseC !<On exit, the lowercase equivalent of the string
    !Local Variables
    INTEGER(INTG), PARAMETER :: offset=(ICHAR("a")-ICHAR("A"))
    INTEGER(INTG) :: characterIdx

    CharacterToLowercaseC=string
    DO characterIdx=1,LEN(string)
      IF(IsUppercase(string(characterIdx:characterIdx))) THEN
        CharacterToLowercaseC(characterIdx:characterIdx)=CHAR(ICHAR(string(characterIdx:characterIdx))+offset)
      ENDIF
    ENDDO !characterIdx

    RETURN
  END FUNCTION CharacterToLowercaseC

  !
  !================================================================================================================================
  !

  !>Returns a character string that is the lowercase equivalent of the supplied varying string.
  PURE FUNCTION CharacterToLowercaseVS(string)

    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: string !<The string to convert
    !Function variable
    CHARACTER(LEN=LEN(string)) :: CharacterToLowercaseVS !<On exit, the lowercase equivalent of the string
    !Local Variables
    INTEGER(INTG), PARAMETER :: offset=(ICHAR("a")-ICHAR("A"))
    INTEGER(INTG) :: characterIdx

    CharacterToLowercaseVS=CHAR(string)
    DO characterIdx=1,LEN(string)
      IF(IsUppercase(CHAR(EXTRACT(string,characterIdx,characterIdx)))) THEN
        CharacterToLowercaseVS(characterIdx:characterIdx)=CHAR(ICHAR(EXTRACT(string,characterIdx,characterIdx))+offset)
      ENDIF
    ENDDO !characterIdx

    RETURN
  END FUNCTION CharacterToLowercaseVS

  !
  !================================================================================================================================
  !

  !>Returns a varying string that is the lowercase equivalent of the supplied character string.
  PURE FUNCTION VStringToLowercaseC(string)

    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: string !<The string to convert
    !Function variable
    TYPE(VARYING_STRING) :: VStringToLowercaseC !<On exit, the lowercase equivalent of the string
    !Local Variables
    INTEGER(INTG), PARAMETER :: offset=(ICHAR("a")-ICHAR("A"))
    INTEGER(INTG) :: characterIdx

    VStringToLowercaseC=string
    DO characterIdx=1,LEN(string)
      IF(IsUppercase(string(characterIdx:characterIdx))) THEN
        VStringToLowercaseC=INSERT(VStringToLowercaseC,characterIdx,CHAR(ICHAR(string(characterIdx:characterIdx))+offset))
      ENDIF
    ENDDO !characterIdx

    RETURN
  END FUNCTION VStringToLowercaseC

  !
  !================================================================================================================================
  !

  !>Returns a varying string that is the lowercase equivalent of the supplied varying string.
  PURE FUNCTION VStringToLowercaseVS(string)

    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: string !<The string to convert
    !Function variable
    TYPE(VARYING_STRING) :: VStringToLowercaseVS !<On exit, the lowercase equivalent of the string
    !Local Variables
    INTEGER(INTG), PARAMETER :: offset=(ICHAR("a")-ICHAR("A"))
    INTEGER(INTG) :: characterIdx

    VStringToLowercaseVS=string
    DO characterIdx=1,LEN(string)
      IF(IsUppercase(CHAR(EXTRACT(string,characterIdx,characterIdx)))) THEN
        VStringToLowercaseVS=INSERT(VStringToLowercaseVS,characterIdx,CHAR(ICHAR(EXTRACT(string,characterIdx,characterIdx))+offset))
      ENDIF
    ENDDO !charaterIdx

    RETURN
  END FUNCTION VStringToLowercaseVS

  !
  !================================================================================================================================
  !

  !>Returns a character string which is uppercase equivalent of the supplied character string.
  PURE FUNCTION CharacterToUppercaseC(string)

    !Argument variables 
    CHARACTER(LEN=*), INTENT(IN) :: string !<The string to convert
    !Function variable
    CHARACTER(LEN=LEN(string)) :: CharacterToUppercaseC !<On exit, the uppercase equivalent of the string
    !Local Variables
    INTEGER(INTG), PARAMETER :: offset=(ICHAR("A")-ICHAR("a"))
    INTEGER(INTG) :: characterIdx

    CharacterToUppercaseC=string
    DO characterIdx=1,LEN(string)
      IF(IsLowercase(string(characterIdx:characterIdx))) THEN        
        CharacterToUppercaseC(characterIdx:characterIdx)=CHAR(ICHAR(string(characterIdx:characterIdx))+offset)
      ENDIF
    ENDDO !i

    RETURN
  END FUNCTION CharacterToUppercaseC

  !
  !================================================================================================================================
  !

  !>Returns a character string which is uppercase equivalent of the supplied varying string.
  PURE FUNCTION CharacterToUppercaseVS(string)

    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: string !<The string to convert
    !Function variable
    CHARACTER(LEN=LEN(string)) :: CharacterToUppercaseVS !<On exit, the uppercase equivalent of the string
    !Local Variables
    INTEGER(INTG), PARAMETER :: offset=(ICHAR("A")-ICHAR("a"))
    INTEGER(INTG) :: characterIdx

    CharacterToUppercaseVS=CHAR(string)
    DO characterIdx=1,LEN(string)
      IF(IsLowercase(CHAR(EXTRACT(string,characterIdx,characterIdx)))) THEN        
        CharacterToUppercaseVS(characterIdx:characterIdx)=CHAR(ICHAR(EXTRACT(string,characterIdx,characterIdx))+offset)
      ENDIF
    ENDDO !characterIdx

    RETURN
  END FUNCTION CharacterToUppercaseVS

  !
  !================================================================================================================================
  !

  !>Returns a varying string which is uppercase equivalent of the supplied character string.
  PURE FUNCTION VStringToUppercaseC(string)

    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: string !<The string to convert
    !Function variable
    TYPE(VARYING_STRING) :: VStringToUppercaseC !<On exit, the uppercase equivalent of the string
    !Local Variables
    INTEGER(INTG), PARAMETER :: offset=(ICHAR("A")-ICHAR("a"))
    INTEGER(INTG) :: characterIdx

    VStringToUppercaseC=string
    DO characterIdx=1,LEN(string)
      IF(IsLowercase(string(characterIdx:characterIdx))) THEN
        VStringToUppercaseC=INSERT(VStringToUppercaseC,characterIdx,CHAR(ICHAR(string(characterIdx:characterIdx))+offset))
      ENDIF
    ENDDO !characterIdx

    RETURN
  END FUNCTION VStringToUppercaseC

  !
  !================================================================================================================================
  !

  !>Returns a varying string which is uppercase equivalent of the supplied varying string.
  PURE FUNCTION VStringToUppercaseVS(string)

    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: string !<The string to convert
    !Function variable
    TYPE(VARYING_STRING) :: VStringToUppercaseVS !<On exit, the uppercase equivalent of the string
    !Local Variables
    INTEGER(INTG), PARAMETER :: offset=(ICHAR("A")-ICHAR("a"))
    INTEGER(INTG) :: characterIdx

    VStringToUppercaseVS=string
    DO characterIdx=1,LEN(string)
      IF(IsLowercase(CHAR(EXTRACT(string,characterIdx,characterIdx)))) THEN
        VStringToUppercaseVS=INSERT(VStringToUppercaseVS,characterIdx,CHAR(ICHAR(EXTRACT(string,characterIdx,characterIdx))+offset))
      ENDIF
    ENDDO !characterIdx

    RETURN
  END FUNCTION VStringToUppercaseVS

  !
  !================================================================================================================================
  !
  
END MODULE Strings
