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
  USE ISO_VARYING_STRING
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
    MODULE PROCEDURE InterfaceCondition_InterfaceEquationsGet
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

  PUBLIC InterfaceCondition_AssertIsFinished,InterfaceCondition_AssertNotFinished

  PUBLIC InterfaceCondition_InterfaceDependentGet

  PUBLIC InterfaceCondition_InterfaceEquationsGet

  PUBLIC INTERFACE_CONDITION_EQUATIONS_GET

  PUBLIC InterfaceCondition_GeometricFieldGet

  PUBLIC InterfaceCondition_InterfaceGet

  PUBLIC InterfaceCondition_LabelGet,InterfaceCondition_LabelSet

  PUBLIC InterfaceCondition_InterfaceLagrangeGet
  
  PUBLIC InterfaceCondition_LagrangeFieldGet

  PUBLIC InterfaceCondition_PenaltyFieldGet

  PUBLIC InterfaceCondition_UserNumberFind

  PUBLIC INTERFACE_CONDITION_USER_NUMBER_FIND

  PUBLIC InterfaceDependent_DependentVariableGet
  
  PUBLIC InterfaceDependent_EquationsSetGet

  PUBLIC InterfaceLagrange_AssertIsFinished,InterfaceLagrange_AssertNotFinished

  PUBLIC InterfaceLagrange_LagrangeFieldGet

CONTAINS

  !
  !================================================================================================================================
  !

  !>Assert that an interface condition has been finished
  SUBROUTINE InterfaceCondition_AssertIsFinished(interfaceCondition,err,error,*)

    !Argument Variables
    TYPE(InterfaceConditionType), POINTER, INTENT(IN) :: interfaceCondition !<The interface condition to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("InterfaceCondition_AssertIsFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)

    IF(.NOT.interfaceCondition%interfaceConditionFinished) THEN
      localError="Interface condition number "//TRIM(NumberToVString(interfaceCondition%userNumber,"*",err,error))// &
        & " has not been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("InterfaceCondition_AssertIsFinished")
    RETURN
999 ERRORSEXITS("InterfaceCondition_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_AssertIsFinished

  !
  !================================================================================================================================
  !

  !>Assert that an interface condition has not been finished
  SUBROUTINE InterfaceCondition_AssertNotFinished(interfaceCondition,err,error,*)

    !Argument Variables
    TYPE(InterfaceConditionType), POINTER, INTENT(IN) :: interfaceCondition !<The interface condition to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("InterfaceCondition_AssertNotFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)

    IF(interfaceCondition%interfaceConditionFinished) THEN
      localError="Interface condition number "//TRIM(NumberToVString(interfaceCondition%userNumber,"*",err,error))// &
        & " has already been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("InterfaceCondition_AssertNotFinished")
    RETURN
999 ERRORSEXITS("InterfaceCondition_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_AssertNotFinished

  !
  !================================================================================================================================
  !

  !>Gets the interface dependent for an interface conditions.
  SUBROUTINE InterfaceCondition_InterfaceDependentGet(interfaceCondition,interfaceDependent,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface conditions to get the interface dependent for
    TYPE(InterfaceDependentType), POINTER :: interfaceDependent !<On exit, a pointer to the interface dependent in the specified interface condition. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceCondition_InterfaceDependentGet",err,error,*998)

    IF(ASSOCIATED(interfaceDependent)) CALL FlagError("Interface dependent is already associated.",err,error,*998)
    CALL InterfaceCondition_AssertIsFinished(interfaceCondition,err,error,*999)
 
    interfaceDependent=>interfaceCondition%dependent
    IF(.NOT.ASSOCIATED(interfaceDependent)) &
      & CALL FlagError("Interface condition dependent is not associated.",err,error,*999)
       
    EXITS("InterfaceCondition_InterfaceDependentGet")
    RETURN
999 NULLIFY(interfaceDependent)
998 ERRORSEXITS("InterfaceCondition_InterfaceDependentGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_InterfaceDependentGet

  !
  !================================================================================================================================
  !

  !>Gets the interface equations for an interface conditions.
  SUBROUTINE InterfaceCondition_InterfaceEquationsGet(interfaceCondition,interfaceEquations,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface conditions to get the interface equations for
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations !<On exit, a pointer to the interface equations in the specified interface condition. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceCondition_InterfaceEquationsGet",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)
    IF(.NOT.interfaceCondition%interfaceConditionFinished) &
      & CALL FlagError("Interface condition has not been finished.",err,error,*999)
    IF(ASSOCIATED(interfaceEquations)) CALL FlagError("Interface equations is already associated.",err,error,*999)

    interfaceEquations=>interfaceCondition%interfaceEquations
    IF(.NOT.ASSOCIATED(interfaceEquations)) &
      & CALL FlagError("Interface condition equations is not associated.",err,error,*999)
       
    EXITS("InterfaceCondition_InterfaceEquationsGet")
    RETURN
999 ERRORSEXITS("InterfaceCondition_InterfaceEquationsGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_InterfaceEquationsGet

  !
  !================================================================================================================================
  !

  !>Gets the geometric field for an interface condition.
  SUBROUTINE InterfaceCondition_GeometricFieldGet(interfaceCondition,geometricField,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition to get the geometric field for
    TYPE(FieldType), POINTER :: geometricField !<On exit, a pointer to the geometric field in the specified interface condition. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("InterfaceCondition_GeometricFieldGet",err,error,*998)

    IF(ASSOCIATED(geometricField)) CALL FlagError("Geometric field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("InterfaceCondition is not associated.",err,error,*999)

    geometricField=>interfaceCondition%geometry%geometricField
    IF(.NOT.ASSOCIATED(geometricField)) THEN
      localError="Geometric field is not associated for interface condition number "// &
      & TRIM(NumberToVString(interfaceCondition%userNumber,"*",err,error))//"."
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
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition to get the interface for
    TYPE(InterfaceType), POINTER :: interface !<On exit, a pointer to the interface in the specified interface condition. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceCondition_InterfaceGet",err,error,*998)

    IF(ASSOCIATED(interface)) CALL FlagError("Interface is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)
    IF(.NOT.interfaceCondition%interfaceConditionFinished) &
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

  !>Returns the label of an interface condition into a character string. \see OpenCMISS::Iron::cmfe_InterfaceCondition_LabelGet
  SUBROUTINE InterfaceCondition_LabelGetC(interfaceCondition,label,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition to get the label for
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

  !>Returns the label of a interface condition into a varying string. \see OpenCMISS::Iron::cmfe_InterfaceCondition_LabelGet
  SUBROUTINE InterfaceCondition_LabelGetVS(interfaceCondition,label,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition to get the label for
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

  !>Sets the label of an interface condition from a character string. \see OpenCMISS::Iron::cmfe_InterfaceCondition_LabelSet
  SUBROUTINE InterfaceCondition_LabelSetC(interfaceCondition,label,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition to set the label for 
    CHARACTER(LEN=*), INTENT(IN) :: label !<The label to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("InterfaceCondition_LabelSetC",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)
    IF(interfaceCondition%interfaceConditionFinished) CALL FlagError("Interface condition has been finished.",err,error,*999)
    
    interfaceCondition%label=label
        
    EXITS("InterfaceCondition_LabelSetC")
    RETURN
999 ERRORSEXITS("InterfaceCondition_LabelSetC",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_LabelSetC

  !
  !================================================================================================================================
  !

  !>Sets the label of an interface condition from a varying string. \see OpenCMISS::Iron::cmfe_InterfaceCondition_LabelSet
  SUBROUTINE InterfaceCondition_LabelSetVS(interfaceCondition,label,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition to set the label for 
    TYPE(VARYING_STRING), INTENT(IN) :: label !<The label to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("InterfaceCondition_LabelSetVS",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)
    IF(interfaceCondition%interfaceConditionFinished) CALL FlagError("Interface condition has been finished.",err,error,*999)
    
    interfaceCondition%label=label
    
    EXITS("InterfaceCondition_LabelSetVS")
    RETURN
999 ERRORSEXITS("InterfaceCondition_LabelSetVS",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_LabelSetVS

  !
  !================================================================================================================================
  !

  !>Gets the interface lagrange for an interface conditions.
  SUBROUTINE InterfaceCondition_InterfaceLagrangeGet(interfaceCondition,interfaceLagrange,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface conditions to get the interface lagrange for
    TYPE(InterfaceLagrangeType), POINTER :: interfaceLagrange !<On exit, a pointer to the interface lagrange in the specified interface condition. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceCondition_InterfaceLagrangeGet",err,error,*998)

    IF(ASSOCIATED(interfaceLagrange)) CALL FlagError("Interface lagrange is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)

    interfaceLagrange=>interfaceCondition%lagrange
    IF(.NOT.ASSOCIATED(interfaceLagrange)) &
      & CALL FlagError("Interface condition Lagrange is not associated.",err,error,*999)
       
    EXITS("InterfaceCondition_InterfaceLagrangeGet")
    RETURN
999 NULLIFY(interfaceLagrange)    
998 ERRORSEXITS("InterfaceCondition_InterfaceLagrangeGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_InterfaceLagrangeGet

  !
  !================================================================================================================================
  !

  !>Gets the Lagrange field for an interface condition.
  SUBROUTINE InterfaceCondition_LagrangeFieldGet(interfaceCondition,lagrangeField,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition to get the Lagrange field for
    TYPE(FieldType), POINTER :: lagrangeField !<On exit, a pointer to the Lagrange field in the specified interface condition. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("InterfaceCondition_LagrangeFieldGet",err,error,*998)

    IF(ASSOCIATED(lagrangeField)) CALL FlagError("Lagrange field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceCondition%lagrange)) THEN
      localError="Lagrange is not associated for interface condition number "// &
      & TRIM(NumberToVString(interfaceCondition%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    lagrangeField=>interfaceCondition%lagrange%lagrangeField
    IF(.NOT.ASSOCIATED(lagrangeField)) THEN
      localError="Lagrange field is not associated for interface condition number "// &
      & TRIM(NumberToVString(interfaceCondition%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("InterfaceCondition_LagrangeFieldGet")
    RETURN
999 NULLIFY(lagrangeField)
998 ERRORSEXITS("InterfaceCondition_LagrangeFieldGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_LagrangeFieldGet

  !
  !================================================================================================================================
  !

  !>Gets the penalty field for an interface condition.
  SUBROUTINE InterfaceCondition_PenaltyFieldGet(interfaceCondition,penaltyField,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition to get the penalty field for
    TYPE(FieldType), POINTER :: penaltyField !<On exit, a pointer to the penalty field in the specified interface condition. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("InterfaceCondition_PenaltyFieldGet",err,error,*998)

    IF(ASSOCIATED(penaltyField)) CALL FlagError("Penalty field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceCondition%penalty)) THEN
      localError="Penalty is not associated for interface condition number "// &
      & TRIM(NumberToVString(interfaceCondition%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    penaltyField=>interfaceCondition%penalty%penaltyField
    IF(.NOT.ASSOCIATED(penaltyField)) THEN
      localError="Penalty field is not associated for interface condition number "// &
      & TRIM(NumberToVString(interfaceCondition%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("InterfaceCondition_PenaltyFieldGet")
    RETURN
999 NULLIFY(penaltyField)
998 ERRORSEXITS("InterfaceCondition_PenaltyFieldGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_PenaltyFieldGet

  !
  !================================================================================================================================
  !

  !>Finds and returns a pointer to the interface condition identified by user number in the given interface. If no interface condition with that user number exists the interface condition is left nullified.
  SUBROUTINE InterfaceCondition_UserNumberFind(userNumber,interface,interfaceCondition,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number to find.
    TYPE(InterfaceType), POINTER :: interface !<The interface to find the interface condition in.
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<On return a pointer to the interface condition with the given user number. If no interface condition with that user number exists then the pointer is returned as NULL. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: interfaceConditionIdx
    TYPE(VARYING_STRING) :: localError

    ENTERS("InterfaceCondition_UserNumberFind",err,error,*999)

    IF(.NOT.ASSOCIATED(interface)) CALL FlagError("Interface is not associated.",err,error,*999)
    IF(ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(interface%interfaceConditions)) THEN
      localError="The interface interface conditions are not associated for interface number "// &
        & TRIM(NumberToVString(interface%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    NULLIFY(interfaceCondition)
    IF(ASSOCIATED(INTERFACE%interfaceConditions%interfaceConditions)) THEN
      DO interfaceConditionIdx=1,INTERFACE%interfaceConditions%numberOfInterfaceConditions
        IF(ASSOCIATED(interface%interfaceConditions%interfaceConditions(interfaceConditionIdx)%PTR)) THEN
          IF(interface%interfaceConditions%interfaceConditions(interfaceConditionIdx)%PTR%userNumber==userNumber) THEN
            interfaceCondition=>interface%interfaceConditions%interfaceConditions(interfaceConditionIdx)%PTR
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

  !>Gets an equations set for an interface dependent.
  SUBROUTINE InterfaceDependent_EquationsSetGet(interfaceDependent,variableIdx,equationsSet,err,error,*)

    !Argument variables
    TYPE(InterfaceDependentType), POINTER :: interfaceDependent !<A pointer to the interface dependent to get the equations set for
    INTEGER(INTG), INTENT(IN) :: variableIdx !<The index of the dependent variable to get the equations set for
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<On exit, a pointer to the specified equations set in the interface dependent. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("InterfaceDependent_EquationsSetGet",err,error,*998)

    IF(ASSOCIATED(equationsSet)) CALL FlagError("Equations set is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceDependent)) CALL FlagError("Interface dependent is not associated.",err,error,*999)
    IF(variableIdx<0.OR.variableIdx>interfaceDependent%numberOfDependentVariables) THEN
      localError="The specified variable index of "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
        & " is invalid. The variable index must be >= 1 and <= "// &
        & TRIM(NumberToVString(interfaceDependent%numberOfDependentVariables,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ASSOCIATED(interfaceDependent%equationsSets)) &
      & CALL FlagError("Interface dependent equations sets is not associated.",err,error,*999)
 
    equationsSet=>interfaceDependent%equationsSets(variableIdx)%ptr      
    IF(.NOT.ASSOCIATED(equationsSet)) THEN
      localError="The equations set for variable index "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)      
    ENDIF
       
    EXITS("InterfaceDependent_EquationsSetGet")
    RETURN
999 NULLIFY(equationsSet)
998 ERRORSEXITS("InterfaceDependent_EquationsSetGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceDependent_EquationsSetGet

  !
  !================================================================================================================================
  !

  !>Gets a dependent field variable for an interface dependent.
  SUBROUTINE InterfaceDependent_DependentVariableGet(interfaceDependent,variableIdx,dependentVariable,err,error,*)

    !Argument variables
    TYPE(InterfaceDependentType), POINTER :: interfaceDependent !<A pointer to the interface dependent to get the dependent variable for
    INTEGER(INTG), INTENT(IN) :: variableIdx !<The index of the dependent variable to get the field variable for
    TYPE(FieldVariableType), POINTER :: dependentVariable !<On exit, a pointer to the specified dependent variable in the interface dependent. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("InterfaceDependent_DependentVariableGet",err,error,*998)

    IF(ASSOCIATED(dependentVariable)) CALL FlagError("Dependent variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceDependent)) CALL FlagError("Interface dependent is not associated.",err,error,*999)
    IF(variableIdx<0.OR.variableIdx>interfaceDependent%numberOfDependentVariables) THEN
      localError="The specified variable index of "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
        & " is invalid. The variable index must be >= 1 and <= "// &
        & TRIM(NumberToVString(interfaceDependent%numberOfDependentVariables,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ASSOCIATED(interfaceDependent%fieldVariables)) &
      & CALL FlagError("Interface dependent field variables is not associated.",err,error,*999)
 
    dependentVariable=>interfaceDependent%fieldVariables(variableIdx)%ptr      
    IF(.NOT.ASSOCIATED(dependentVariable)) THEN
      localError="The dependent field variable for variable index "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)      
    ENDIF
       
    EXITS("InterfaceDependent_DependentVariableGet")
    RETURN
999 NULLIFY(dependentVariable)
998 ERRORSEXITS("InterfaceDependent_DependentVariableGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceDependent_DependentVariableGet

  !
  !================================================================================================================================
  !

  !>Assert that an interface lagrange has been finished
  SUBROUTINE InterfaceLagrange_AssertIsFinished(interfaceLagrange,err,error,*)

    !Argument Variables
    TYPE(InterfaceLagrangeType), POINTER, INTENT(IN) :: interfaceLagrange !<The interface lagrange to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceLagrange_AssertIsFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceLagrange)) CALL FlagError("Interface Lagrange is not associated.",err,error,*999)

    IF(.NOT.interfaceLagrange%lagrangeFinished) &
      & CALL FlagError("Interface Lagrange has not been finished.",err,error,*999)
    
    EXITS("InterfaceLagrange_AssertIsFinished")
    RETURN
999 ERRORSEXITS("InterfaceLagrange_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceLagrange_AssertIsFinished

  !
  !================================================================================================================================
  !

  !>Assert that an interface lagrange has not been finished
  SUBROUTINE InterfaceLagrange_AssertNotFinished(interfaceLagrange,err,error,*)

    !Argument Variables
    TYPE(InterfaceLagrangeType), POINTER, INTENT(IN) :: interfaceLagrange !<The interface lagrange to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceLagrange_AssertNotFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceLagrange)) CALL FlagError("Interface Lagrange is not associated.",err,error,*999)

    IF(interfaceLagrange%lagrangeFinished)  &
      & CALL FlagError("Interface Lagrange has already been finished.",err,error,*999)
    
    EXITS("InterfaceLagrange_AssertNotFinished")
    RETURN
999 ERRORSEXITS("InterfaceLagrange_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceLagrange_AssertNotFinished

  !
  !================================================================================================================================
  !

  !>Gets the Lagrange field for an interface Lagrange.
  SUBROUTINE InterfaceLagrange_LagrangeFieldGet(interfaceLagrange,lagrangeField,err,error,*)

    !Argument variables
    TYPE(InterfaceLagrangeType), POINTER :: interfaceLagrange !<A pointer to the interface Lagrange to get the Lagrange field for
    TYPE(FieldType), POINTER :: lagrangeField !<On exit, a pointer to the Lagrange field in the specified interface Lagrange. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("InterfaceLagrange_LagrangeFieldGet",err,error,*998)

    IF(ASSOCIATED(lagrangeField)) CALL FlagError("Lagrange field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceLagrange)) CALL FlagError("Interface Lagrange is not associated.",err,error,*999)

    lagrangeField=>interfaceLagrange%lagrangeField
    IF(.NOT.ASSOCIATED(lagrangeField)) CALL FlagError("The interface Lagrange Lagrange field is not associated.",err,error,*999)
       
    EXITS("InterfaceLagrange_LagrangeFieldGet")
    RETURN
999 NULLIFY(lagrangeField)
998 ERRORSEXITS("InterfaceLagrange_LagrangeFieldGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceLagrange_LagrangeFieldGet

  !
  !================================================================================================================================
  !
  
END MODULE InterfaceConditionAccessRoutines
