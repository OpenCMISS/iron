!> \file
!> \author Chris Bradley
!> \brief This module contains all interface condition access method routines.
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

!> This module contains all interface condition access method routines.
MODULE InterfaceConditionAccessRoutines
  
  USE BaseRoutines
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

  INTERFACE INTERFACE_CONDITION_EQUATIONS_GET
    MODULE PROCEDURE InterfaceCondition_EquationsGet
  END INTERFACE INTERFACE_CONDITION_EQUATIONS_GET

  INTERFACE InterfaceCondition_LabelGet
    MODULE PROCEDURE InterfaceCondition_LabelGetC
    MODULE PROCEDURE InterfaceCondition_LabelGetVS
  END INTERFACE InterfaceCondition_LabelGet

  INTERFACE InterfaceCondition_LabelSet
    MODULE PROCEDURE InterfaceCondition_LabelSetC
    MODULE PROCEDURE InterfaceCondition_LabelSetVS
  END INTERFACE InterfaceCondition_LabelSet

  INTERFACE INTERFACE_CONDITION_USER_NUMBER_FIND
    MODULE PROCEDURE InterfaceCondition_UserNumberFind
  END INTERFACE INTERFACE_CONDITION_USER_NUMBER_FIND

  PUBLIC InterfaceCondition_EquationsGet

  PUBLIC INTERFACE_CONDITION_EQUATIONS_GET

  PUBLIC InterfaceCondition_GeometricFieldGet

  PUBLIC InterfaceCondition_InterfaceGet

  PUBLIC InterfaceCondition_LabelGet,InterfaceCondition_LabelSet

  PUBLIC InterfaceCondition_UserNumberFind

  PUBLIC INTERFACE_CONDITION_USER_NUMBER_FIND

CONTAINS

  !
  !================================================================================================================================
  !

  !>Gets the interface equations for an interface conditions.
  SUBROUTINE InterfaceCondition_EquationsGet(interfaceCondition,interfaceEquations,err,error,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: interfaceCondition !<A pointer to the interface conditions to get the interface equations for
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: interfaceEquations !<On exit, a pointer to the interface equations in the specified interface condition. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceCondition_EquationsGet",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)
    IF(.NOT.interfaceCondition%INTERFACE_CONDITION_FINISHED) &
      & CALL FlagError("Interface condition has not been finished.",err,error,*999)
    IF(ASSOCIATED(interfaceEquations)) CALL FlagError("Interface equations is already associated.",err,error,*999)

    interfaceEquations=>interfaceCondition%INTERFACE_EQUATIONS
    IF(.NOT.ASSOCIATED(interfaceEquations)) &
      & CALL FlagError("Interface condition equations is not associated.",err,error,*999)
       
    EXITS("InterfaceCondition_EquationsGet")
    RETURN
999 ERRORSEXITS("InterfaceCondition_EquationsGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_EquationsGet

  !
  !================================================================================================================================
  !

  !>Gets the geometric field for an interface condition.
  SUBROUTINE InterfaceCondition_GeometricFieldGet(interfaceCondition,geometricField,err,error,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: interfaceCondition !<A pointer to the interface condition to get the geometric field for
    TYPE(FIELD_TYPE), POINTER :: geometricField !<On exit, a pointer to the geometric field in the specified interface condition. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("InterfaceCondition_GeometricFieldGet",err,error,*998)

    IF(ASSOCIATED(geometricField)) CALL FlagError("Geometric field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("InterfaceCondition is not associated.",err,error,*999)

    geometricField=>interfaceCondition%geometry%GEOMETRIC_FIELD
    IF(.NOT.ASSOCIATED(geometricField)) THEN
      localError="Geometric field is not associated for interface condition number "// &
      & TRIM(NumberToVString(interfaceCondition%USER_NUMBER,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("InterfaceCondition_GeometricFieldGet")
    RETURN
999 NULLIFY(geometricField)
998 ERRORSEXITS("InterfaceCondition_GeometricFieldGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_GeometricFieldGet

  !
  !================================================================================================================================
  !

  !>Gets the interface for an interface condition.
  SUBROUTINE InterfaceCondition_InterfaceGet(interfaceCondition,interface,err,error,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: interfaceCondition !<A pointer to the interface condition to get the interface for
    TYPE(INTERFACE_TYPE), POINTER :: interface !<On exit, a pointer to the interface in the specified interface condition. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceCondition_InterfaceGet",err,error,*998)

    IF(ASSOCIATED(interface)) CALL FlagError("Interface is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)
    IF(.NOT.interfaceCondition%INTERFACE_CONDITION_FINISHED) &
      & CALL FlagError("Interface condition has not been finished.",err,error,*999)
 
    INTERFACE=>interfaceCondition%INTERFACE
    IF(.NOT.ASSOCIATED(interface)) CALL FlagError("Interface condition interface is not associated.",err,error,*999)
       
    EXITS("InterfaceCondition_InterfaceGet")
    RETURN
999 NULLIFY(interface)
998 ERRORSEXITS("InterfaceCondition_InterfaceGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_InterfaceGet

  !
  !================================================================================================================================
  !

  !>Returns the label of an interface condition into a character string. \see OpenCMISS::cmfe_InterfaceCondition_LabelGet
  SUBROUTINE InterfaceCondition_LabelGetC(interfaceCondition,label,err,error,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: interfaceCondition !<A pointer to the interface condition to get the label for
    CHARACTER(LEN=*), INTENT(OUT) :: label !<On return, the interface condition label.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: cLength,vsLength

    ENTERS("InterfaceCondition_LabelGetC",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)
    
    cLength=LEN(label)
    vsLength=LEN_TRIM(interfaceCondition%label)
    IF(cLength>vsLength) THEN
      label=CHAR(interfaceCondition%label,vsLength)
    ELSE
      label=CHAR(interfaceCondition%label,cLength)
    ENDIF
    
    EXITS("InterfaceCondition_LabelGetC")
    RETURN
999 ERRORSEXITS("InterfaceCondition_LabelGetC",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_LabelGetC

   !
  !================================================================================================================================
  !

  !>Returns the label of a interface condition into a varying string. \see OpenCMISS::cmfe_InterfaceCondition_LabelGet
  SUBROUTINE InterfaceCondition_LabelGetVS(interfaceCondition,label,err,error,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: interfaceCondition !<A pointer to the interface condition to get the label for
    TYPE(VARYING_STRING), INTENT(OUT) :: label !<On return, the interface condition label.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("InterfaceCondition_LabelGetVS",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)
    
    label=VAR_STR(CHAR(interfaceCondition%label))
          
    EXITS("InterfaceCondition_LabelGetVS")
    RETURN
999 ERRORSEXITS("InterfaceCondition_LabelGetVS",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_LabelGetVS

  !
  !================================================================================================================================
  !

  !>Sets the label of an interface condition from a character string. \see OpenCMISS::cmfe_InterfaceCondition_LabelSet
  SUBROUTINE InterfaceCondition_LabelSetC(interfaceCondition,label,err,error,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: interfaceCondition !<A pointer to the interface condition to set the label for 
    CHARACTER(LEN=*), INTENT(IN) :: label !<The label to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("InterfaceCondition_LabelSetC",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)
    IF(interfaceCondition%INTERFACE_CONDITION_FINISHED) CALL FlagError("Interface condition has been finished.",err,error,*999)
    
    interfaceCondition%label=label
        
    EXITS("InterfaceCondition_LabelSetC")
    RETURN
999 ERRORSEXITS("InterfaceCondition_LabelSetC",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_LabelSetC

  !
  !================================================================================================================================
  !

  !>Sets the label of an interface condition from a varying string. \see OpenCMISS::cmfe_InterfaceCondition_LabelSet
  SUBROUTINE InterfaceCondition_LabelSetVS(interfaceCondition,label,err,error,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: interfaceCondition !<A pointer to the interface condition to set the label for 
    TYPE(VARYING_STRING), INTENT(IN) :: label !<The label to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("InterfaceCondition_LabelSetVS",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)
    IF(interfaceCondition%INTERFACE_CONDITION_FINISHED) CALL FlagError("Interface condition has been finished.",err,error,*999)
    
    interfaceCondition%label=label
    
    EXITS("InterfaceCondition_LabelSetVS")
    RETURN
999 ERRORSEXITS("InterfaceCondition_LabelSetVS",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_LabelSetVS

  !
  !================================================================================================================================
  !

  !>Finds and returns a pointer to the interface condition identified by user number in the given interface. If no interface condition with that user number exists the interface condition is left nullified.
  SUBROUTINE InterfaceCondition_UserNumberFind(userNumber,interface,interfaceCondition,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number to find.
    TYPE(INTERFACE_TYPE), POINTER :: interface !<The interface to find the interface condition in.
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: interfaceCondition !<On return a pointer to the interface condition with the given user number. If no interface condition with that user number exists then the pointer is returned as NULL. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: interfaceConditionIdx
    TYPE(VARYING_STRING) :: localError

    ENTERS("InterfaceCondition_UserNumberFind",err,error,*999)

    IF(.NOT.ASSOCIATED(interface)) CALL FlagError("Interface is not associated.",err,error,*999)
    IF(ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(interface%INTERFACE_CONDITIONS)) THEN
      localError="The interface interface conditions are not associated for interface number "// &
        & TRIM(NumberToVString(interface%USER_NUMBER,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    NULLIFY(interfaceCondition)
    IF(ASSOCIATED(INTERFACE%INTERFACE_CONDITIONS%INTERFACE_CONDITIONS)) THEN
      DO interfaceConditionIdx=1,INTERFACE%INTERFACE_CONDITIONS%NUMBER_OF_INTERFACE_CONDITIONS
        IF(ASSOCIATED(interface%INTERFACE_CONDITIONS%INTERFACE_CONDITIONS(interfaceConditionIdx)%PTR)) THEN
          IF(interface%INTERFACE_CONDITIONS%INTERFACE_CONDITIONS(interfaceConditionIdx)%PTR%USER_NUMBER==userNumber) THEN
            interfaceCondition=>interface%INTERFACE_CONDITIONS%INTERFACE_CONDITIONS(interfaceConditionIdx)%PTR
            EXIT
          ENDIF
        ELSE
          localError="The interface condition pointer in interface conditions is not associated for interface condition index "// &
            & TRIM(NumberToVString(interfaceConditionIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDDO !interfaceConditionIdx
    ENDIF
    
    EXITS("InterfaceCondition_UserNumberFind")
    RETURN
999 ERRORSEXITS("InterfaceCondition_UserNumberFind",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_UserNumberFind

  !
  !================================================================================================================================
  !

  
END MODULE InterfaceConditionAccessRoutines
