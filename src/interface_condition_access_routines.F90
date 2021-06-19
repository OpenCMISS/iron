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


  !> \addtogroup InterfaceConditions_Methods InterfaceConditions::Constants::Methods
  !> \brief Interface condition methods.
  !> \see InterfaceConditions
  !>@{
  INTEGER(INTG), PARAMETER :: INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD=1 !<Lagrange multipliers interface condition method. \see InterfaceConditions_Methods,InterfaceConditions
  INTEGER(INTG), PARAMETER :: INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD=2 !<Augmented Lagrange multiplers interface condition method. \see InterfaceConditions_Methods,InterfaceConditions
  INTEGER(INTG), PARAMETER :: INTERFACE_CONDITION_PENALTY_METHOD=3 !<Penalty interface condition method. \see InterfaceConditions_Methods,InterfaceConditions
  INTEGER(INTG), PARAMETER :: INTERFACE_CONDITION_POINT_TO_POINT_METHOD=4 !<Point to point interface condition method. \see InterfaceConditions_Methods,InterfaceConditions
  !>@}

  !> \addtogroup InterfaceConditions_Operators InterfaceConditions::Constants::Operators
  !> \brief Interface condition operators.
  !> \see InterfaceConditions
  !>@{
  INTEGER(INTG), PARAMETER :: INTERFACE_CONDITION_FIELD_CONTINUITY_OPERATOR=1 !<Continuous field operator, i.e., lambda.(u_1-u_2). \see InterfaceConditions_Operators,InterfaceConditions 
  INTEGER(INTG), PARAMETER :: INTERFACE_CONDITION_FIELD_NORMAL_CONTINUITY_OPERATOR=2 !<Continuous field normal operator, i.e., lambda(u_1.n_1-u_2.n_2). \see InterfaceConditions_Operators,InterfaceConditions 
  INTEGER(INTG), PARAMETER :: INTERFACE_CONDITION_FLS_CONTACT_OPERATOR=3 !<Frictionless contact operator, i.e., lambda.(x_1.n-x_2.n). \see InterfaceConditions_Operators,InterfaceConditions 
  INTEGER(INTG), PARAMETER :: INTERFACE_CONDITION_FLS_CONTACT_REPROJECT_OPERATOR=4 !<Frictionless contact operator, reproject at each newton iteration and has geometric linearisation terms i.e., lambda.(x_1.n-x_2.n). \see InterfaceConditions_Operators,InterfaceConditions 
  INTEGER(INTG), PARAMETER :: INTERFACE_CONDITION_SOLID_FLUID_OPERATOR=5 !<Solid fluid operator, i.e., lambda.(v_f-du_s/dt). \see InterfaceConditions_Operators,InterfaceConditions 
  INTEGER(INTG), PARAMETER :: INTERFACE_CONDITION_SOLID_FLUID_NORMAL_OPERATOR=6 !<Solid fluid normal operator, i.e., lambda(v_f.n_f-du_s/dt.n_s). \see InterfaceConditions_Operators,InterfaceConditions 
  !>@}


  !> \addtogroup InterfaceConditions_OutputTypes InterfaceConditions::Constants::OutputTypes
  !> \brief The interface conditions output types
  !> \see InterfaceConditions,OpenCMISS_InterfaceConditionsConstants
  !>@{
  INTEGER(INTG), PARAMETER :: INTERFACE_CONDITION_NO_OUTPUT=0 !<No output. \see InterfaceConditions_OutputTypes,InterfaceConditions
  INTEGER(INTG), PARAMETER :: INTERFACE_CONDITION_PROGRESS_OUTPUT=1 !<Progress information output. \see InterfaceConditions_OutputTypes,InterfaceConditions
  !>@}

  !> \addtogroup InterfaceConditions_IntegrationType InterfaceConditions::Constants::IntegrationType
  !> \brief Interface condition IntegrationType.
  !> \see InterfaceConditions
  !>@{
  INTEGER(INTG), PARAMETER :: INTERFACE_CONDITION_GAUSS_INTEGRATION=1 !<Gauss points integration type, i.e. Loop over element Gauss points and sum up their contribution. \see InterfaceConditions_IntegrationType,InterfaceConditions 
  INTEGER(INTG), PARAMETER :: INTERFACE_CONDITION_DATA_POINTS_INTEGRATION=2 !< Data points integration type i.e. Loop over data points and  sum up their contribution.\see InterfaceConditions_IntegrationType,InterfaceConditions 
  !>@}
  
  !Module types

  !Module variables

  !Interfaces

  INTERFACE InterfaceCondition_LabelGet
    MODULE PROCEDURE InterfaceCondition_LabelGetC
    MODULE PROCEDURE InterfaceCondition_LabelGetVS
  END INTERFACE InterfaceCondition_LabelGet

  INTERFACE InterfaceCondition_LabelSet
    MODULE PROCEDURE InterfaceCondition_LabelSetC
    MODULE PROCEDURE InterfaceCondition_LabelSetVS
  END INTERFACE InterfaceCondition_LabelSet

  PUBLIC INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD, &
    & INTERFACE_CONDITION_PENALTY_METHOD,INTERFACE_CONDITION_POINT_TO_POINT_METHOD

  PUBLIC INTERFACE_CONDITION_FIELD_CONTINUITY_OPERATOR,INTERFACE_CONDITION_FIELD_NORMAL_CONTINUITY_OPERATOR, &
    & INTERFACE_CONDITION_FLS_CONTACT_OPERATOR,INTERFACE_CONDITION_FLS_CONTACT_REPROJECT_OPERATOR, &
    & INTERFACE_CONDITION_SOLID_FLUID_OPERATOR,INTERFACE_CONDITION_SOLID_FLUID_NORMAL_OPERATOR

  PUBLIC INTERFACE_CONDITION_NO_OUTPUT,INTERFACE_CONDITION_PROGRESS_OUTPUT

  PUBLIC INTERFACE_CONDITION_GAUSS_INTEGRATION,INTERFACE_CONDITION_DATA_POINTS_INTEGRATION

  PUBLIC InterfaceCondition_AssertIsFinished,InterfaceCondition_AssertNotFinished

  PUBLIC InterfaceCondition_InterfaceDependentGet

  PUBLIC InterfaceCondition_InterfaceEquationsGet

  PUBLIC InterfaceCondition_GeometricFieldGet

  PUBLIC InterfaceCondition_GlobalNumberGet

  PUBLIC InterfaceCondition_InterfaceGet

  PUBLIC InterfaceCondition_LabelGet,InterfaceCondition_LabelSet

  PUBLIC InterfaceCondition_IntegrationTypeGet

  PUBLIC InterfaceCondition_InterfaceLagrangeGet
  
  PUBLIC InterfaceCondition_LagrangeFieldGet

  PUBLIC InterfaceCondition_MethodGet

  PUBLIC InterfaceCondition_OperatorGet

  PUBLIC InterfaceCondition_OutputTypeGet

  PUBLIC InterfaceCondition_InterfacePenaltyGet

  PUBLIC InterfaceCondition_PenaltyFieldExists

  PUBLIC InterfaceCondition_PenaltyFieldGet

  PUBLIC InterfaceCondition_UserNumberFind

  PUBLIC InterfaceCondition_UserNumberGet

  PUBLIC InterfaceDependent_DependentVariableGet
  
  PUBLIC InterfaceDependent_EquationsSetGet

  PUBLIC InterfaceDependent_InterfaceConditionGet

  PUBLIC InterfaceDependent_NumberOfDependentVariablesGet

  PUBLIC InterfaceDependent_VariableMeshIndexGet

  PUBLIC InterfaceLagrange_AssertIsFinished,InterfaceLagrange_AssertNotFinished

  PUBLIC InterfaceLagrange_LagrangeFieldGet

  PUBLIC InterfacePenalty_AssertIsFinished,InterfacePenalty_AssertNotFinished

  PUBLIC InterfacePenalty_PenaltyFieldGet

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

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)
#endif    

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

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)
#endif    

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

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(interfaceDependent)) CALL FlagError("Interface dependent is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)
#endif    
 
    interfaceDependent=>interfaceCondition%dependent

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(interfaceDependent)) &
      & CALL FlagError("Interface condition dependent is not associated.",err,error,*999)
#endif    
       
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
 
    ENTERS("InterfaceCondition_InterfaceEquationsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(interfaceEquations)) CALL FlagError("Interface equations is already associated.",err,error,*998)
    CALL InterfaceCondition_AssertIsFinished(interfaceCondition,err,error,*999)
#endif    
 
    interfaceEquations=>interfaceCondition%interfaceEquations

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(interfaceEquations)) &
      & CALL FlagError("Interface condition equations is not associated.",err,error,*999)
#endif    
       
    EXITS("InterfaceCondition_InterfaceEquationsGet")
    RETURN
999 NULLIFY(interfaceEquations)
998 ERRORSEXITS("InterfaceCondition_InterfaceEquationsGet",err,error)
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
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("InterfaceCondition_GeometricFieldGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(geometricField)) CALL FlagError("Geometric field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("InterfaceCondition is not associated.",err,error,*999)
#endif    

    geometricField=>interfaceCondition%geometry%geometricField

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(geometricField)) THEN
      localError="Geometric field is not associated for interface condition number "// &
      & TRIM(NumberToVString(interfaceCondition%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
    EXITS("InterfaceCondition_GeometricFieldGet")
    RETURN
999 NULLIFY(geometricField)
998 ERRORSEXITS("InterfaceCondition_GeometricFieldGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_GeometricFieldGet

  !
  !================================================================================================================================
  !

  !>Gets the global number for an interface condition.
  SUBROUTINE InterfaceCondition_GlobalNumberGet(interfaceCondition,globalNumber,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition to get the global number for
    INTEGER(INTG), INTENT(OUT) :: globalNumber !<On exit, the global number for the interface condition
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceCondition_GlobalNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)
#endif    

    globalNumber=interfaceCondition%globalNumber

    EXITS("InterfaceCondition_GlobalNumberGet")
    RETURN
999 ERRORSEXITS("InterfaceCondition_GlobalNumberGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_GlobalNumberGet

  !
  !================================================================================================================================
  !

  !>Gets the integration type for an interface condition.
  SUBROUTINE InterfaceCondition_IntegrationTypeGet(interfaceCondition,integrationType,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition to get the integration type for
    INTEGER(INTG), INTENT(OUT) :: integrationType !<On exit, the intergration type for the interface condition \see InterfaceConditions_IntegrationType
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceCondition_IntegrationTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)
#endif    

    integrationType=interfaceCondition%integrationType

    EXITS("InterfaceCondition_IntegrationTypeGet")
    RETURN
999 ERRORSEXITS("InterfaceCondition_IntegrationTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_IntegrationTypeGet

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

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(INTERFACE)) CALL FlagError("Interface is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)
#endif    
 
    INTERFACE=>interfaceCondition%INTERFACE

#ifdef WITH_POSTCHECKS      
    IF(.NOT.ASSOCIATED(INTERFACE)) CALL FlagError("Interface condition interface is not associated.",err,error,*999)
#endif      
       
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

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)
#endif    
    
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

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)
#endif    
    
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
    
#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)
#endif
    
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

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)
#endif    
    
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

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(interfaceLagrange)) CALL FlagError("Interface lagrange is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)
#endif    

    interfaceLagrange=>interfaceCondition%lagrange

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(interfaceLagrange)) &
      & CALL FlagError("Interface condition Lagrange is not associated.",err,error,*999)
#endif    
       
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
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("InterfaceCondition_LagrangeFieldGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(lagrangeField)) CALL FlagError("Lagrange field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(interfaceCondition%lagrange)) THEN
      localError="Lagrange is not associated for interface condition number "// &
      & TRIM(NumberToVString(interfaceCondition%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    lagrangeField=>interfaceCondition%lagrange%lagrangeField

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(lagrangeField)) THEN
      localError="Lagrange field is not associated for interface condition number "// &
      & TRIM(NumberToVString(interfaceCondition%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
    EXITS("InterfaceCondition_LagrangeFieldGet")
    RETURN
999 NULLIFY(lagrangeField)
998 ERRORSEXITS("InterfaceCondition_LagrangeFieldGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_LagrangeFieldGet

  !
  !================================================================================================================================
  !

  !>Gets the method for an interface condition.
  SUBROUTINE InterfaceCondition_MethodGet(interfaceCondition,method,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition to get the method for
    INTEGER(INTG), INTENT(OUT) :: method !<On exit, the method for the interface condition \see InterfaceConditions_Methods
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceCondition_MethodGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)
#endif    

    method=interfaceCondition%method

    EXITS("InterfaceCondition_MethodGet")
    RETURN
999 ERRORSEXITS("InterfaceCondition_MethodGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_MethodGet

  !
  !================================================================================================================================
  !

  !>Gets the operator for an interface condition.
  SUBROUTINE InterfaceCondition_OperatorGet(interfaceCondition,operator,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition to get the operator for
    INTEGER(INTG), INTENT(OUT) :: operator !<On exit, the operator for the interface condition \see InterfaceConditions_Operators
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceCondition_OperatorGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)
#endif    

    operator=interfaceCondition%operator

    EXITS("InterfaceCondition_OperatorGet")
    RETURN
999 ERRORSEXITS("InterfaceCondition_OperatorGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_OperatorGet

  !
  !================================================================================================================================
  !

  !>Gets the output type for an interface condition.
  SUBROUTINE InterfaceCondition_OutputTypeGet(interfaceCondition,outputType,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition to get the output type for
    INTEGER(INTG), INTENT(OUT) :: outputType !<On exit, the output type for the interface condition \see InterfaceConditions_OutputTypes
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceCondition_OutputTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)
#endif    

    outputType=interfaceCondition%outputType

    EXITS("InterfaceCondition_OutputTypeGet")
    RETURN
999 ERRORSEXITS("InterfaceCondition_OutputTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_OutputTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the interface penalty for an interface conditions.
  SUBROUTINE InterfaceCondition_InterfacePenaltyGet(interfaceCondition,interfacePenalty,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface conditions to get the interface penalty for
    TYPE(InterfacePenaltyType), POINTER :: interfacePenalty !<On exit, a pointer to the interface penalty in the specified interface condition. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceCondition_InterfacePenaltyGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(interfacePenalty)) CALL FlagError("Interface penalty is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)
#endif    

    interfacePenalty=>interfaceCondition%penalty

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(interfacePenalty)) &
      & CALL FlagError("Interface condition penalty is not associated.",err,error,*999)
#endif    
       
    EXITS("InterfaceCondition_InterfacePenaltyGet")
    RETURN
999 NULLIFY(interfacePenalty)    
998 ERRORSEXITS("InterfaceCondition_InterfacePenaltyGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_InterfacePenaltyGet

  !
  !================================================================================================================================
  !

  !>Checks the penalty field exists for an interface condition.
  SUBROUTINE InterfaceCondition_PenaltyFieldExists(interfaceCondition,penaltyField,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition to check the penalty field for
    TYPE(FieldType), POINTER :: penaltyField !<On exit, a pointer to the penalty field in the specified interface condition if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif
    
    ENTERS("InterfaceCondition_PenaltyFieldExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(penaltyField)) CALL FlagError("Penalty field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)
#endif    
    
    IF(ASSOCIATED(interfaceCondition%penalty)) THEN
      penaltyField=>interfaceCondition%penalty%penaltyField
    ELSE
      NULLIFY(penaltyField)
    ENDIF
       
    EXITS("InterfaceCondition_PenaltyFieldExists")
    RETURN
999 NULLIFY(penaltyField)
998 ERRORSEXITS("InterfaceCondition_PenaltyFieldExists",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_PenaltyFieldExists

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
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("InterfaceCondition_PenaltyFieldGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(penaltyField)) CALL FlagError("Penalty field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(interfaceCondition%penalty)) THEN
      localError="Penalty is not associated for interface condition number "// &
      & TRIM(NumberToVString(interfaceCondition%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    penaltyField=>interfaceCondition%penalty%penaltyField

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(penaltyField)) THEN
      localError="Penalty field is not associated for interface condition number "// &
      & TRIM(NumberToVString(interfaceCondition%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
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
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    

    ENTERS("InterfaceCondition_UserNumberFind",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(interface)) CALL FlagError("Interface is not associated.",err,error,*999)
    IF(ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(INTERFACE%interfaceConditions)) THEN
      localError="The interface interface conditions are not associated for interface number "// &
        & TRIM(NumberToVString(interface%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    IF(ASSOCIATED(INTERFACE%interfaceConditions%interfaceConditions)) THEN
      DO interfaceConditionIdx=1,INTERFACE%interfaceConditions%numberOfInterfaceConditions
#ifdef WITH_PRECHECKS        
        IF(.NOT.ASSOCIATED(INTERFACE%interfaceConditions%interfaceConditions(interfaceConditionIdx)%ptr)) THEN
          localError="The interface condition pointer in interface conditions is not associated for interface condition index "// &
            & TRIM(NumberToVString(interfaceConditionIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
#endif        
        IF(INTERFACE%interfaceConditions%interfaceConditions(interfaceConditionIdx)%PTR%userNumber==userNumber) THEN
          interfaceCondition=>interface%interfaceConditions%interfaceConditions(interfaceConditionIdx)%ptr
          EXIT
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

  !>Gets the user number for an interface condition.
  SUBROUTINE InterfaceCondition_UserNumberGet(interfaceCondition,userNumber,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition to get the user number for
    INTEGER(INTG), INTENT(OUT) :: userNumber !<On exit, the user number for the interface condition
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceCondition_UserNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)
#endif    

    userNumber=interfaceCondition%userNumber

    EXITS("InterfaceCondition_UserNumberGet")
    RETURN
999 ERRORSEXITS("InterfaceCondition_UserNumberGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceCondition_UserNumberGet

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

#ifdef WITH_PRECHECKS    
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
#endif    
 
    dependentVariable=>interfaceDependent%fieldVariables(variableIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(dependentVariable)) THEN
      localError="The dependent field variable for variable index "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)      
    ENDIF
#endif    
       
    EXITS("InterfaceDependent_DependentVariableGet")
    RETURN
999 NULLIFY(dependentVariable)
998 ERRORSEXITS("InterfaceDependent_DependentVariableGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceDependent_DependentVariableGet

  !
  !================================================================================================================================
  !

  !>Gets an equations set for an interface dependent.
  SUBROUTINE InterfaceDependent_EquationsSetGet(interfaceDependent,variableIdx,equationsSet,err,error,*)

    !Argument variables
    TYPE(InterfaceDependentType), POINTER :: interfaceDependent !<A pointer to the interface dependent to get the equations set for
    INTEGER(INTG), INTENT(IN) :: variableIdx !<The index of the dependent variable to get the equations set for
    TYPE(EquationsSetType), POINTER :: equationsSet !<On exit, a pointer to the specified equations set in the interface dependent. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("InterfaceDependent_EquationsSetGet",err,error,*998)

#ifdef WITH_PRECHECKS 
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
#endif    
 
    equationsSet=>interfaceDependent%equationsSets(variableIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(equationsSet)) THEN
      localError="The equations set for variable index "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)      
    ENDIF
#endif    
       
    EXITS("InterfaceDependent_EquationsSetGet")
    RETURN
999 NULLIFY(equationsSet)
998 ERRORSEXITS("InterfaceDependent_EquationsSetGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceDependent_EquationsSetGet

  !
  !================================================================================================================================
  !

  !>Gets the interface condition for an interface dependent.
  SUBROUTINE InterfaceDependent_InterfaceConditionGet(interfaceDependent,interfaceCondition,err,error,*)

    !Argument variables
    TYPE(InterfaceDependentType), POINTER :: interfaceDependent !<A pointer to the interface dependent to get the interface condition for
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<On exit, a pointer to the interface condition for the interface dependent. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("InterfaceDependent_InterfaceConditionGet",err,error,*998)

#ifdef WITH_PRECHECKS 
    IF(ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceDependent)) CALL FlagError("Interface dependent is not associated.",err,error,*999)
#endif    
 
    interfaceCondition=>interfaceDependent%interfaceCondition

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(interfaceCondition)) &
      & CALL FlagError("The interface condition is not associated for the interface dependent.",err,error,*999)      
#endif    
       
    EXITS("InterfaceDependent_InterfaceConditionGet")
    RETURN
999 NULLIFY(interfaceCondition)
998 ERRORSEXITS("InterfaceDependent_InterfaceConditionGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceDependent_InterfaceConditionGet

  !
  !================================================================================================================================
  !

  !>Gets the number of dependent variables for an interface dependent.
  SUBROUTINE InterfaceDependent_NumberOfDependentVariablesGet(interfaceDependent,numberOfDependentVariables,err,error,*)

    !Argument variables
    TYPE(InterfaceDependentType), POINTER :: interfaceDependent !<A pointer to the interface dependent to get the number of dependent variables for
    INTEGER(INTG), INTENT(OUT) :: numberOfDependentVariables !<On exit, the number of dependent variables for the interface dependent.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("InterfaceDependent_NumberOfDependentVariablesGet",err,error,*999)

#ifdef WITH_PRECHECKS 
    IF(.NOT.ASSOCIATED(interfaceDependent)) CALL FlagError("Interface dependent is not associated.",err,error,*999)
#endif    
 
    numberOfDependentVariables=interfaceDependent%numberOfDependentVariables
       
    EXITS("InterfaceDependent_NumberOfDependentVariablesGet")
    RETURN
999 ERRORSEXITS("InterfaceDependent_NumberOfDependentVariablesGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceDependent_NumberOfDependentVariablesGet

  !
  !================================================================================================================================
  !

  !>Gets the variable mesh index for an interface condition dependent information.
  SUBROUTINE InterfaceDependent_VariableMeshIndexGet(interfaceDependent,variableIdx,variableMeshIndex,err,error,*)

    !Argument variables
    TYPE(InterfaceDependentType), POINTER :: interfaceDependent !<A pointer to the interface dependent to get the variable mesh index for
    INTEGER(INTG), INTENT(IN) :: variableIdx !<The variable index to get the variable mesh index for
    INTEGER(INTG), INTENT(OUT) :: variableMeshIndex !<On exit, the variable mesh index for the interface dependent variable
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("InterfaceCondition_VariableMeshIndexGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(interfaceDependent)) CALL FlagError("Interface dependent is not associated.",err,error,*999)
    IF(variableIdx<1.OR.variableIdx>interfaceDependent%numberOfDependentVariables) THEN
      localError="The specified variable index of "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
        & " is invalid. The variable index should be >= 1 and <= "// &
        & TRIM(NumberToVString(interfaceDependent%numberOfDependentVariables,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ASSOCIATED(interfaceDependent%variableMeshIndices)) &
      & CALL FlagError("The variable mesh indices is not associated for the interface dependent.",err,error,*999)
#endif    

    variableMeshIndex=interfaceDependent%variableMeshIndices(variableIdx)

    EXITS("InterfaceCondition_VariableMeshIndexGet")
    RETURN
999 ERRORSEXITS("InterfaceCondition_VariableMeshIndexGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceDependent_VariableMeshIndexGet

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

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(interfaceLagrange)) CALL FlagError("Interface Lagrange is not associated.",err,error,*999)
#endif    

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

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(interfaceLagrange)) CALL FlagError("Interface Lagrange is not associated.",err,error,*999)
#endif    

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

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(lagrangeField)) CALL FlagError("Lagrange field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceLagrange)) CALL FlagError("Interface Lagrange is not associated.",err,error,*999)
#endif    

    lagrangeField=>interfaceLagrange%lagrangeField

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(lagrangeField)) CALL FlagError("The interface Lagrange Lagrange field is not associated.",err,error,*999)
#endif    
       
    EXITS("InterfaceLagrange_LagrangeFieldGet")
    RETURN
999 NULLIFY(lagrangeField)
998 ERRORSEXITS("InterfaceLagrange_LagrangeFieldGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceLagrange_LagrangeFieldGet

  !
  !================================================================================================================================
  !

  !>Assert that an interface penalty has been finished
  SUBROUTINE InterfacePenalty_AssertIsFinished(interfacePenalty,err,error,*)

    !Argument Variables
    TYPE(InterfacePenaltyType), POINTER, INTENT(IN) :: interfacePenalty !<The interface penalty to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfacePenalty_AssertIsFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(interfacePenalty)) CALL FlagError("Interface penalty is not associated.",err,error,*999)
#endif    

    IF(.NOT.interfacePenalty%penaltyFinished) &
      & CALL FlagError("Interface Penalty has not been finished.",err,error,*999)
    
    EXITS("InterfacePenalty_AssertIsFinished")
    RETURN
999 ERRORSEXITS("InterfacePenalty_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE InterfacePenalty_AssertIsFinished

  !
  !================================================================================================================================
  !

  !>Assert that an interface penalty has not been finished
  SUBROUTINE InterfacePenalty_AssertNotFinished(interfacePenalty,err,error,*)

    !Argument Variables
    TYPE(InterfacePenaltyType), POINTER, INTENT(IN) :: interfacePenalty !<The interface penalty to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfacePenalty_AssertNotFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(interfacePenalty)) CALL FlagError("Interface penalty is not associated.",err,error,*999)
#endif    

    IF(interfacePenalty%penaltyFinished)  &
      & CALL FlagError("Interface penalty has already been finished.",err,error,*999)
    
    EXITS("InterfacePenalty_AssertNotFinished")
    RETURN
999 ERRORSEXITS("InterfacePenalty_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE InterfacePenalty_AssertNotFinished

  !
  !================================================================================================================================
  !

  !>Gets the Penalty field for an interface Penalty.
  SUBROUTINE InterfacePenalty_PenaltyFieldGet(interfacePenalty,penaltyField,err,error,*)

    !Argument variables
    TYPE(InterfacePenaltyType), POINTER :: interfacePenalty !<A pointer to the interface penalty to get the penalty field for
    TYPE(FieldType), POINTER :: penaltyField !<On exit, a pointer to the penalty field in the specified interface penalty. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("InterfacePenalty_PenaltyFieldGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(penaltyField)) CALL FlagError("Penalty field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfacePenalty)) CALL FlagError("Interface penalty is not associated.",err,error,*999)
#endif    

    penaltyField=>interfacePenalty%penaltyField

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(penaltyField)) CALL FlagError("The interface penalty penalty field is not associated.",err,error,*999)
#endif    
       
    EXITS("InterfacePenalty_PenaltyFieldGet")
    RETURN
999 NULLIFY(penaltyField)
998 ERRORSEXITS("InterfacePenalty_PenaltyFieldGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfacePenalty_PenaltyFieldGet

  !
  !================================================================================================================================
  !
  
END MODULE InterfaceConditionAccessRoutines
