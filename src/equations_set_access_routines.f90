!> \file
!> \author Chris Bradley
!> \brief This module contains all equations set access method routines.
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

!> This module contains all equations set access method routines.
MODULE EquationsSetAccessRoutines
  
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
 
  INTERFACE EQUATIONS_SET_ANALYTIC_TIME_GET
    MODULE PROCEDURE EquationsSet_AnalyticTimeGet
  END INTERFACE EQUATIONS_SET_ANALYTIC_TIME_GET
  
  INTERFACE EQUATIONS_SET_ANALYTIC_TIME_SET
    MODULE PROCEDURE EquationsSet_AnalyticTimeSet
  END INTERFACE EQUATIONS_SET_ANALYTIC_TIME_SET
  
  INTERFACE EQUATIONS_SET_EQUATIONS_GET
    MODULE PROCEDURE EquationsSet_EquationsGet
  END INTERFACE EQUATIONS_SET_EQUATIONS_GET
  
  INTERFACE EquationsSet_LabelGet
    MODULE PROCEDURE EquationsSet_LabelGetC
    MODULE PROCEDURE EquationsSet_LabelGetVS
  END INTERFACE EquationsSet_LabelGet

   INTERFACE EquationsSet_LabelSet
    MODULE PROCEDURE EquationsSet_LabelSetC
    MODULE PROCEDURE EquationsSet_LabelSetVS
  END INTERFACE EquationsSet_LabelSet

  INTERFACE EQUATIONS_SET_USER_NUMBER_FIND
    MODULE PROCEDURE EquationsSet_UserNumberFind
  END INTERFACE EQUATIONS_SET_USER_NUMBER_FIND

  PUBLIC EquationsSet_AnalyticCreated
  
  PUBLIC EquationsSet_AnalyticFieldExists
  
  PUBLIC EquationsSet_AnalyticFieldGet

  PUBLIC EquationsSet_AnalyticTimeGet,EquationsSet_AnalyticTimeSet

  PUBLIC EQUATIONS_SET_ANALYTIC_TIME_GET,EQUATIONS_SET_ANALYTIC_TIME_SET
  
  PUBLIC EquationsSet_CoordinateSystemGet

  PUBLIC EquationsSet_DerivedCreated
  
  PUBLIC EquationsSet_DerivedFieldExists
  
  PUBLIC EquationsSet_DerivedFieldGet
  
  PUBLIC EquationsSet_DependentFieldGet

  PUBLIC EquationsSet_EquationsGet

  PUBLIC EQUATIONS_SET_EQUATIONS_GET
  
  PUBLIC EquationsSet_EquationsSetFieldGet
  
  PUBLIC EquationsSet_FibreFieldExists
  
  PUBLIC EquationsSet_GeometricFieldGet

  PUBLIC EquationsSet_IndependentCreated
  
  PUBLIC EquationsSet_IndependentFieldExists
  
  PUBLIC EquationsSet_IndependentFieldGet
  
  PUBLIC EquationsSet_LabelGet,EquationsSet_LabelSet

  PUBLIC EquationsSet_MaterialsCreated

  PUBLIC EquationsSet_MaterialsFieldExists

  PUBLIC EquationsSet_MaterialsFieldGet

  PUBLIC EquationsSet_SourceCreated
  
  PUBLIC EquationsSet_SourceFieldExists
  
  PUBLIC EquationsSet_SourceFieldGet
  
  PUBLIC EquationsSet_RegionGet
  
  PUBLIC EquationsSet_UserNumberFind

  PUBLIC EQUATIONS_SET_USER_NUMBER_FIND

CONTAINS

  !
  !================================================================================================================================
  !

  !>Asserts that analytic has been created for an equations set.
  SUBROUTINE EquationsSet_AnalyticCreated(equationsSet,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to check that the analytic has been created for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_AnalyticCreated",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet%analytic)) THEN
      localError="Equations set analytic has not been created for equations set number "// &
        & TRIM(NumberToVString(equationsSet%USER_NUMBER,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
     
    EXITS("EquationsSet_AnalyticCreated")
    RETURN
999 ERRORSEXITS("EquationsSet_AnalyticCreated",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AnalyticCreated

  !
  !================================================================================================================================
  !

  !>Gets the analytic field for an equations set if it exists.
  SUBROUTINE EquationsSet_AnalyticFieldExists(equationsSet,analyticField,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to get the analytic field for
    TYPE(FIELD_TYPE), POINTER :: analyticField !<On exit, a pointer to the analytic field in the specified equations set if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsSet_AnalyticFieldExists",err,error,*998)

    IF(ASSOCIATED(analyticField)) CALL FlagError("Analytic field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    IF(ASSOCIATED(equationsSet%analytic)) THEN
      analyticField=>equationsSet%analytic%ANALYTIC_FIELD
    ELSE
      NULLIFY(analyticField)
    ENDIF
     
    EXITS("EquationsSet_AnalyticFieldExists")
    RETURN
999 NULLIFY(analyticField)
998 ERRORSEXITS("EquationsSet_AnalyticFieldExists",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AnalyticFieldExists

  !
  !================================================================================================================================
  !

  !>Gets the analytic field for an equations set.
  SUBROUTINE EquationsSet_AnalyticFieldGet(equationsSet,analyticField,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to get the analytic field for
    TYPE(FIELD_TYPE), POINTER :: analyticField !<On exit, a pointer to the analytic field in the specified equations set. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_AnalyticFieldGet",err,error,*998)

    IF(ASSOCIATED(analyticField)) CALL FlagError("Analytic field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet%analytic)) THEN
      localError="Analytic information is not associated for equations set number "// &
        & TRIM(NumberToVString(equationsSet%USER_NUMBER,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    analyticField=>equationsSet%analytic%ANALYTIC_FIELD
    IF(.NOT.ASSOCIATED(analyticField)) THEN
      localError="Analytic field is not associated for equations set number "// &
        & TRIM(NumberToVString(equationsSet%USER_NUMBER,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("EquationsSet_AnalyticFieldGet")
    RETURN
999 NULLIFY(analyticField)
998 ERRORSEXITS("EquationsSet_AnalyticFieldGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AnalyticFieldGet

  !
  !================================================================================================================================
  !

  !>Returns the analytic time for an equations set. \see OpenCMISS::cmfe_EquationsSet_AnalyticTimeGet
  SUBROUTINE EquationsSet_AnalyticTimeGet(equationsSet,analyticTime,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to get the time for.
    REAL(DP), INTENT(OUT) :: analyticTime !<On return, the analytic time value.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsSet_AnalyticTimeGet",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)    
    IF(.NOT.ASSOCIATED(equationsSet%analytic)) CALL FlagError("Equations set analytic is not associated.",err,error,*999)

    analyticTime=equationsSet%analytic%ANALYTIC_TIME
      
    EXITS("EquationsSet_AnalyticTimeGet")
    RETURN
999 ERRORSEXITS("EquationsSet_AnalyticTimeGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AnalyticTimeGet

  !
  !================================================================================================================================
  !

  !>Sets the analytic time for an equations set. \see OpenCMISS::cmfe_EquationsSet_AnalyticTimeSet
  SUBROUTINE EquationsSet_AnalyticTimeSet(equationsSet,analyticTime,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to set the time for.
    REAL(DP), INTENT(IN) :: analyticTime !<The analytic time value to set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsSet_AnalyticTimeSet",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)    
    IF(.NOT.ASSOCIATED(equationsSet%analytic)) CALL FlagError("Equations set analytic is not associated.",err,error,*999)

    equationsSet%analytic%ANALYTIC_TIME=analyticTime
      
    EXITS("EquationsSet_AnalyticTimeSet")
    RETURN
999 ERRORSEXITS("EquationsSet_AnalyticTimeSet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AnalyticTimeSet

  !
  !================================================================================================================================
  !

  !>Gets the coordinate system used for an equations set.
  SUBROUTINE EquationsSet_CoordinateSystemGet(equationsSet,coordinateSystem,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to get the coordinate system for
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: coordinateSystem !<On exit, a pointer to the coordinate system for the specified equations set. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(REGION_TYPE), POINTER :: region
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_CoordinateSystemGet",err,error,*998)

    IF(ASSOCIATED(coordinateSystem)) CALL FlagError("Coordinate system is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    region=>equationsSet%region
    IF(.NOT.ASSOCIATED(region)) THEN
      localError="Region is not associated for equations set number "// &
      & TRIM(NumberToVString(equationsSet%USER_NUMBER,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    coordinateSystem=>region%COORDINATE_SYSTEM
    IF(.NOT.ASSOCIATED(coordinateSystem)) THEN
      localError="Coordinate system is not associated for region number "// &
        & TRIM(NumberToVString(region%USER_NUMBER,"*",err,error))//" of equations set number "// &
        & TRIM(NumberToVString(equationsSet%USER_NUMBER,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("EquationsSet_CoordinateSystemGet")
    RETURN
999 NULLIFY(coordinateSystem)
998 ERRORSEXITS("EquationsSet_CoordinateSystemGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_CoordinateSystemGet

  !
  !================================================================================================================================
  !

  !>Asserts that derived has been created for an equations set.
  SUBROUTINE EquationsSet_DerivedCreated(equationsSet,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to check that the derived has been created for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_DerivedCreated",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet%derived)) THEN
      localError="Equations set derived has not been created for equations set number "// &
        & TRIM(NumberToVString(equationsSet%USER_NUMBER,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
     
    EXITS("EquationsSet_DerivedCreated")
    RETURN
999 ERRORSEXITS("EquationsSet_DerivedCreated",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_DerivedCreated

  !
  !================================================================================================================================
  !

  !>Gets the derived field for an equations set if it exists.
  SUBROUTINE EquationsSet_DerivedFieldExists(equationsSet,derivedField,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to get the derived field for
    TYPE(FIELD_TYPE), POINTER :: derivedField !<On exit, a pointer to the derived field in the specified equations set if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsSet_DerivedFieldExists",err,error,*998)

    IF(ASSOCIATED(derivedField)) CALL FlagError("Derived field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    IF(ASSOCIATED(equationsSet%derived)) THEN
      derivedField=>equationsSet%derived%derivedField
    ELSE
      NULLIFY(derivedField)
    ENDIF
      
    EXITS("EquationsSet_DerivedFieldExists")
    RETURN
999 NULLIFY(derivedField)
998 ERRORSEXITS("EquationsSet_DerivedFieldExists",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_DerivedFieldExists

  !
  !================================================================================================================================
  !

  !>Gets the derived field for an equations set.
  SUBROUTINE EquationsSet_DerivedFieldGet(equationsSet,derivedField,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to get the derived field for
    TYPE(FIELD_TYPE), POINTER :: derivedField !<On exit, a pointer to the derived field in the specified equations set. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_DerivedFieldGet",err,error,*998)

    IF(ASSOCIATED(derivedField)) CALL FlagError("Derived field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equationsSet%derived)) CALL FlagError("Equations set derived is not associated.",err,error,*999)
    
    derivedField=>equationsSet%derived%derivedField
    IF(.NOT.ASSOCIATED(derivedField)) THEN
      localError="Derived field is not associated for equations set number "// &
      & TRIM(NumberToVString(equationsSet%USER_NUMBER,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("EquationsSet_DerivedFieldGet")
    RETURN
999 NULLIFY(derivedField)
998 ERRORSEXITS("EquationsSet_DerivedFieldGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_DerivedFieldGet

  !
  !================================================================================================================================
  !

  !>Gets the dependent field for an equations set.
  SUBROUTINE EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to get the dependent field for
    TYPE(FIELD_TYPE), POINTER :: dependentField !<On exit, a pointer to the dependent field in the specified equations set. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_DependentFieldGet",err,error,*998)

    IF(ASSOCIATED(dependentField)) CALL FlagError("Dependent field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    dependentField=>equationsSet%dependent%DEPENDENT_FIELD
    IF(.NOT.ASSOCIATED(dependentField)) THEN
      localError="Dependent field is not associated for equations set number "// &
      & TRIM(NumberToVString(equationsSet%USER_NUMBER,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("EquationsSet_DependentFieldGet")
    RETURN
999 NULLIFY(dependentField)
998 ERRORSEXITS("EquationsSet_DependentFieldGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_DependentFieldGet

  !
  !================================================================================================================================
  !

  !>Gets the equations for an equations set.
  SUBROUTINE EquationsSet_EquationsGet(equationsSet,equations,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to get the equations for
    TYPE(EquationsType), POINTER :: equations !<On exit, a pointer to the equations in the specified equations set. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_EquationsGet",err,error,*998)

    IF(ASSOCIATED(equations)) CALL FlagError("Equations is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.equationsSet%EQUATIONS_SET_FINISHED) CALL FlagError("Equations set has not been finished.",err,error,*999)

    equations=>equationsSet%equations
    IF(.NOT.ASSOCIATED(equations)) THEN
      localError="Equations is not associated for equations set number "// &
      & TRIM(NumberToVString(equationsSet%USER_NUMBER,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("EquationsSet_EquationsGet")
    RETURN
999 NULLIFY(equations)
998 ERRORSEXITS("EquationsSet_EquationsGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_EquationsGet

  !
  !================================================================================================================================
  !

  !>Gets the equations set field for an equations set.
  SUBROUTINE EquationsSet_EquationsSetFieldGet(equationsSet,equationsSetField,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to get the equations set field for
    TYPE(FIELD_TYPE), POINTER :: equationsSetField !<On exit, a pointer to the equations set field in the specified equations set. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_EquationsSetFieldGet",err,error,*998)

    IF(ASSOCIATED(equationsSetField)) CALL FlagError("Equations set field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    equationsSetField=>equationsSet%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
    IF(.NOT.ASSOCIATED(equationsSetField)) THEN
      localError="Equations set field is not associated for equations set number "// &
      & TRIM(NumberToVString(equationsSet%USER_NUMBER,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("EquationsSet_EquationsSetFieldGet")
    RETURN
999 NULLIFY(equationsSetField)
998 ERRORSEXITS("EquationsSet_EquationsSetFieldGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_EquationsSetFieldGet

  !
  !================================================================================================================================
  !

  !>Gets the fibre field for an equations set if it exists.
  SUBROUTINE EquationsSet_FibreFieldExists(equationsSet,fibreField,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to get the fibre field for
    TYPE(FIELD_TYPE), POINTER :: fibreField !<On exit, a pointer to the fibre field in the specified equations set if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsSet_FibreFieldExists",err,error,*998)

    IF(ASSOCIATED(fibreField)) CALL FlagError("Fibre field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    fibreField=>equationsSet%geometry%FIBRE_FIELD
      
    EXITS("EquationsSet_FibreFieldExists")
    RETURN
999 NULLIFY(fibreField)
998 ERRORSEXITS("EquationsSet_FibreFieldExists",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_FibreFieldExists

  !
  !================================================================================================================================
  !

  !>Gets the geometric field for an equations set.
  SUBROUTINE EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to get the geometric field for
    TYPE(FIELD_TYPE), POINTER :: geometricField !<On exit, a pointer to the geometric field in the specified equations set. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_GeometricFieldGet",err,error,*998)

    IF(ASSOCIATED(geometricField)) CALL FlagError("Geometric field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    geometricField=>equationsSet%geometry%GEOMETRIC_FIELD
    IF(.NOT.ASSOCIATED(geometricField)) THEN
      localError="Geometric field is not associated for equations set number "// &
        & TRIM(NumberToVString(equationsSet%USER_NUMBER,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
      
    EXITS("EquationsSet_GeometricFieldGet")
    RETURN
999 NULLIFY(geometricField)
998 ERRORSEXITS("EquationsSet_GeometricFieldGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_GeometricFieldGet

  !
  !================================================================================================================================
  !

  !>Asserts that independent has been created for an equations set.
  SUBROUTINE EquationsSet_IndependentCreated(equationsSet,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to check that the independent has been created for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_IndependentCreated",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet%independent)) THEN
      localError="Equations set independent has not been created for equations set number "// &
        & TRIM(NumberToVString(equationsSet%USER_NUMBER,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
     
    EXITS("EquationsSet_IndependentCreated")
    RETURN
999 ERRORSEXITS("EquationsSet_IndependentCreated",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_IndependentCreated

  !
  !================================================================================================================================
  !

  !>Gets the independent field for an equations set if it exists.
  SUBROUTINE EquationsSet_IndependentFieldExists(equationsSet,independentField,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to get the independent field for
    TYPE(FIELD_TYPE), POINTER :: independentField !<On exit, a pointer to the independent field in the specified equations set if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsSet_IndependentFieldExists",err,error,*998)

    IF(ASSOCIATED(independentField)) CALL FlagError("Independent field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    IF(ASSOCIATED(equationsSet%independent)) THEN
      independentField=>equationsSet%INDEPENDENT%INDEPENDENT_FIELD
    ELSE
      NULLIFY(independentField)
    ENDIF
       
    EXITS("EquationsSet_IndependentFieldExists")
    RETURN
999 NULLIFY(independentField)
998 ERRORSEXITS("EquationsSet_IndependentFieldExists",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_IndependentFieldExists

  !
  !================================================================================================================================
  !

  !>Gets the independent field for an equations set.
  SUBROUTINE EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to get the independent field for
    TYPE(FIELD_TYPE), POINTER :: independentField !<On exit, a pointer to the independent field in the specified equations set. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_IndependentFieldGet",err,error,*998)

    IF(ASSOCIATED(independentField)) CALL FlagError("Independent field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet%independent)) THEN
      localError="Independent information is not associated for equations set number "// &
        & TRIM(NumberToVString(equationsSet%USER_NUMBER,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    independentField=>equationsSet%independent%INDEPENDENT_FIELD
    IF(.NOT.ASSOCIATED(independentField)) THEN
      localError="Independent field is not associated for equations set number "// &
      & TRIM(NumberToVString(equationsSet%USER_NUMBER,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("EquationsSet_IndependentFieldGet")
    RETURN
999 NULLIFY(independentField)
998 ERRORSEXITS("EquationsSet_IndependentFieldGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_IndependentFieldGet

  !
  !================================================================================================================================
  !

  !>Returns the label of an equations set into a character string. \see OpenCMISS::cmfe_EquationsSet_LabelGet
  SUBROUTINE EquationsSet_LabelGetC(equationsSet,label,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to get the label for
    CHARACTER(LEN=*), INTENT(OUT) :: label !<On return, the equations set label.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: cLength,vsLength

    ENTERS("EquationsSet_LabelGetC",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    
    cLength=LEN(label)
    vsLength=LEN_TRIM(equationsSet%label)
    IF(cLength>vsLength) THEN
      label=CHAR(equationsSet%label,vsLength)
    ELSE
      label=CHAR(equationsSet%label,cLength)
    ENDIF
    
    EXITS("EquationsSet_LabelGetC")
    RETURN
999 ERRORSEXITS("EquationsSet_LabelGetC",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_LabelGetC

   !
  !================================================================================================================================
  !

  !>Returns the label of a equations set into a varying string. \see OpenCMISS::cmfe_EquationsSet_LabelGet
  SUBROUTINE EquationsSet_LabelGetVS(equationsSet,label,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to get the label for
    TYPE(VARYING_STRING), INTENT(OUT) :: label !<On return, the equations set label.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsSet_LabelGetVS",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    
    label=VAR_STR(CHAR(equationsSet%label))
          
    EXITS("EquationsSet_LabelGetVS")
    RETURN
999 ERRORSEXITS("EquationsSet_LabelGetVS",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_LabelGetVS

  !
  !================================================================================================================================
  !

  !>Sets the label of an equations set from a character string. \see OpenCMISS::cmfe_EquationsSet_LabelSet
  SUBROUTINE EquationsSet_LabelSetC(equationsSet,label,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to set the label for 
    CHARACTER(LEN=*), INTENT(IN) :: label !<The label to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsSet_LabelSetC",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(equationsSet%EQUATIONS_SET_FINISHED) CALL FlagError("Equations set has been finished.",err,error,*999)
    
    equationsSet%label=label
        
    EXITS("EquationsSet_LabelSetC")
    RETURN
999 ERRORSEXITS("EquationsSet_LabelSetC",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_LabelSetC

  !
  !================================================================================================================================
  !

  !>Sets the label of an equations set from a varying string. \see OpenCMISS::cmfe_EquationsSet_LabelSet
  SUBROUTINE EquationsSet_LabelSetVS(equationsSet,label,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to set the label for 
    TYPE(VARYING_STRING), INTENT(IN) :: label !<The label to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsSet_LabelSetVS",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(equationsSet%EQUATIONS_SET_FINISHED) CALL FlagError("Equations set has been finished.",err,error,*999)
    
    equationsSet%label=label
    
    EXITS("EquationsSet_LabelSetVS")
    RETURN
999 ERRORSEXITS("EquationsSet_LabelSetVS",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_LabelSetVS

  !
  !================================================================================================================================
  !

  !>Asserts that materials has been created for an equations set.
  SUBROUTINE EquationsSet_MaterialsCreated(equationsSet,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to check that the materials has been created for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_MaterialsCreated",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet%materials)) THEN
      localError="Equations set materials has not been created for equations set number "// &
        & TRIM(NumberToVString(equationsSet%USER_NUMBER,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
     
    EXITS("EquationsSet_MaterialsCreated")
    RETURN
999 ERRORSEXITS("EquationsSet_MaterialsCreated",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_MaterialsCreated

  !
  !================================================================================================================================
  !

  !>Gets the materials field for an equations set if it exists.
  SUBROUTINE EquationsSet_MaterialsFieldExists(equationsSet,materialsField,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to get the materials field for
    TYPE(FIELD_TYPE), POINTER :: materialsField !<On exit, a pointer to the materials field in the specified equations set if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsSet_MaterialsFieldExists",err,error,*998)

    IF(ASSOCIATED(materialsField)) CALL FlagError("Materials field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    IF(ASSOCIATED(equationsSet%materials)) THEN
      materialsField=>equationsSet%materials%MATERIALS_FIELD
    ELSE
      NULLIFY(materialsField)
    ENDIF
     
    EXITS("EquationsSet_MaterialsFieldExists")
    RETURN
999 NULLIFY(materialsField)
998 ERRORSEXITS("EquationsSet_MaterialsFieldExists",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_MaterialsFieldExists

  !
  !================================================================================================================================
  !

  !>Gets the materials field for an equations set.
  SUBROUTINE EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to get the materials field for
    TYPE(FIELD_TYPE), POINTER :: materialsField !<On exit, a pointer to the materials field in the specified equations set. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_MaterialsFieldGet",err,error,*998)

    IF(ASSOCIATED(materialsField)) CALL FlagError("Materials field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet%materials)) THEN
      localError="Materials information is not associated for equations set number "// &
        & TRIM(NumberToVString(equationsSet%USER_NUMBER,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    materialsField=>equationsSet%materials%MATERIALS_FIELD
    IF(.NOT.ASSOCIATED(materialsField)) THEN
      localError="Materials field is not associated for equations set number "// &
      & TRIM(NumberToVString(equationsSet%USER_NUMBER,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("EquationsSet_MaterialsFieldGet")
    RETURN
999 NULLIFY(materialsField)
998 ERRORSEXITS("EquationsSet_MaterialsFieldGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_MaterialsFieldGet

  !
  !================================================================================================================================
  !

  !>Gets the region used for an equations set.
  SUBROUTINE EquationsSet_RegionGet(equationsSet,region,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to get the region for
    TYPE(REGION_TYPE), POINTER :: region !<On exit, a pointer to the region for the specified equations set. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_RegionGet",err,error,*998)

    IF(ASSOCIATED(region)) CALL FlagError("Region is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    region=>equationsSet%region
    IF(.NOT.ASSOCIATED(region)) THEN
      localError="Region is not associated for equations set number "// &
      & TRIM(NumberToVString(equationsSet%USER_NUMBER,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("EquationsSet_RegionGet")
    RETURN
999 NULLIFY(region)
998 ERRORSEXITS("EquationsSet_RegionGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_RegionGet

  !
  !================================================================================================================================
  !

  !>Asserts that source has been created for an equations set.
  SUBROUTINE EquationsSet_SourceCreated(equationsSet,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to check that the source has been created for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_SourceCreated",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet%source)) THEN
      localError="Equations set source has not been created for equations set number "// &
        & TRIM(NumberToVString(equationsSet%USER_NUMBER,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
     
    EXITS("EquationsSet_SourceCreated")
    RETURN
999 ERRORSEXITS("EquationsSet_SourceCreated",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_SourceCreated

  !
  !================================================================================================================================
  !

  !>Gets the source field for an equations set if it exists
  SUBROUTINE EquationsSet_SourceFieldExists(equationsSet,sourceField,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to get the source field for
    TYPE(FIELD_TYPE), POINTER :: sourceField !<On exit, a pointer to the source field in the specified equations set if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsSet_SourceFieldExists",err,error,*998)

    IF(ASSOCIATED(sourceField)) CALL FlagError("Source field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    IF(ASSOCIATED(equationsSet%source)) THEN
      sourceField=>equationsSet%source%SOURCE_FIELD
    ELSE
      NULLIFY(sourceField)
    ENDIF
       
    EXITS("EquationsSet_SourceFieldExists")
    RETURN
999 NULLIFY(sourceField)
998 ERRORSEXITS("EquationsSet_SourceFieldExists",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_SourceFieldExists

  !
  !================================================================================================================================
  !

  !>Gets the source  field for an equations set.
  SUBROUTINE EquationsSet_SourceFieldGet(equationsSet,sourceField,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to get the source field for
    TYPE(FIELD_TYPE), POINTER :: sourceField !<On exit, a pointer to the source field in the specified equations set. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_SourceFieldGet",err,error,*998)

    IF(ASSOCIATED(sourceField)) CALL FlagError("Source field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet%source)) THEN
      localError="Source information is not associated for equations set number "// &
        & TRIM(NumberToVString(equationsSet%USER_NUMBER,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    sourceField=>equationsSet%source%SOURCE_FIELD
    IF(.NOT.ASSOCIATED(sourceField)) THEN
      localError="Source field is not associated for equations set number "// &
      & TRIM(NumberToVString(equationsSet%USER_NUMBER,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("EquationsSet_SourceFieldGet")
    RETURN
999 NULLIFY(sourceField)
998 ERRORSEXITS("EquationsSet_SourceFieldGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_SourceFieldGet

  !
  !================================================================================================================================
  !

  !>Finds and returns a pointer to the equations set identified by user number in the given region. If no equations set with that user number exists the equations set is left nullified.
  SUBROUTINE EquationsSet_UserNumberFind(userNumber,region,equationsSet,err,error,*)

    !Argument variables 
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the equation set to find.
    TYPE(REGION_TYPE), POINTER :: region !<The region to find the equations set in
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<On return, a pointer to the equations set if an equations set with the specified user number exists in the given region. If no equation set with the specified number exists a NULL pointer is returned. The pointer must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsSet_UserNumberFind",err,error,*999)

    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated.",err,error,*999)
    IF(ASSOCIATED(equationsSet)) CALL FlagError("Equations set is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(region%EQUATIONS_SETS)) THEN
      localError="The equations sets on region number "//TRIM(NumberToVString(region%USER_NUMBER,"*",err,error))// &
        & " are not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    NULLIFY(equationsSet)
    IF(ASSOCIATED(region%EQUATIONS_SETS%EQUATIONS_SETS)) THEN
      DO equationsSetIdx=1,region%EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS
        IF(ASSOCIATED(region%EQUATIONS_SETS%EQUATIONS_SETS(equationsSetIdx)%ptr)) THEN
          IF(region%EQUATIONS_SETS%EQUATIONS_SETS(equationsSetIdx)%ptr%USER_NUMBER==userNumber) THEN
            equationsSet=>region%EQUATIONS_SETS%EQUATIONS_SETS(equationsSetIdx)%ptr
            EXIT
          ENDIF
        ELSE
          localError="The equations set pointer in region equations sets is not associated for equations set index "// &
            & TRIM(NumberToVString(equationsSetIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDDO !equationsSetIdx
    ENDIF
    
    EXITS("EquationsSet_UserNumberFind")
    RETURN
999 ERRORSEXITS("EquationsSet_UserNumberFind",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_UserNumberFind

  !
  !================================================================================================================================
  !

END MODULE EquationsSetAccessRoutines
