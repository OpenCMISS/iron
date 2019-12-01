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

  PUBLIC EquationsSet_AnalyticFieldExists
  
  PUBLIC EquationsSet_AnalyticFieldGet
  
  PUBLIC EquationsSet_AnalyticFunctionTypeGet

  PUBLIC EquationsSet_AnalyticTimeGet,EquationsSet_AnalyticTimeSet

  PUBLIC EQUATIONS_SET_ANALYTIC_TIME_GET,EQUATIONS_SET_ANALYTIC_TIME_SET

  PUBLIC EquationsSet_AssertIsFinished,EquationsSet_AssertNotFinished
  
  PUBLIC EquationsSet_AssertAnalyticIsCreated,EquationsSet_AssertAnalyticNotCreated
  
  PUBLIC EquationsSet_AssertAnalyticIsFinished,EquationsSet_AssertAnalyticNotFinished
  
  PUBLIC EquationsSet_AssertDependentIsFinished,EquationsSet_AssertDependentNotFinished
  
  PUBLIC EquationsSet_AssertDerivedIsCreated,EquationsSet_AssertDerivedNotCreated
  
  PUBLIC EquationsSet_AssertDerivedIsFinished,EquationsSet_AssertDerivedNotFinished
  
  PUBLIC EquationsSet_AssertEquationsIsCreated,EquationsSet_AssertEquationsNotCreated
  
  PUBLIC EquationsSet_AssertEquationsIsFinished,EquationsSet_AssertEquationsNotFinished
  
  PUBLIC EquationsSet_AssertIndependentIsCreated,EquationsSet_AssertIndependentNotCreated
  
  PUBLIC EquationsSet_AssertIndependentIsFinished,EquationsSet_AssertIndependentNotFinished
  
  PUBLIC EquationsSet_AssertMaterialsIsCreated,EquationsSet_AssertMaterialsNotCreated
  
  PUBLIC EquationsSet_AssertMaterialsIsFinished,EquationsSet_AssertMaterialsNotFinished
  
  PUBLIC EquationsSet_AssertSourceIsCreated,EquationsSet_AssertSourceNotCreated
  
  PUBLIC EquationsSet_AssertSourceIsFinished,EquationsSet_AssertSourceNotFinished
  
  PUBLIC EquationsSet_CoordinateSystemGet

  PUBLIC EquationsSet_DerivedFieldExists
  
  PUBLIC EquationsSet_DerivedFieldGet
  
  PUBLIC EquationsSet_DependentFieldGet

  PUBLIC EquationsSet_EquationsGet

  PUBLIC EQUATIONS_SET_EQUATIONS_GET
  
  PUBLIC EquationsSet_EquationsSetsGet

  PUBLIC EquationsSet_EquationsSetFieldGet
  
  PUBLIC EquationsSet_FibreFieldExists
  
  PUBLIC EquationsSet_FibreFieldGet
  
  PUBLIC EquationsSet_GeometricFieldGet

  PUBLIC EquationsSet_IndependentFieldExists
  
  PUBLIC EquationsSet_IndependentFieldGet
  
  PUBLIC EquationsSet_LabelGet,EquationsSet_LabelSet

  PUBLIC EquationsSet_MaterialsFieldExists

  PUBLIC EquationsSet_MaterialsFieldGet

  PUBLIC EquationsSet_OutputTypeGet

  PUBLIC EquationsSet_RegionGet

  PUBLIC EquationsSet_SolutionMethodGet
  
  PUBLIC EquationsSet_SourceFieldExists
  
  PUBLIC EquationsSet_SourceFieldGet

  PUBLIC EquationsSet_SpecificationGet
  
  PUBLIC EquationsSet_SpecificationSizeGet

  PUBLIC EquationsSet_TimesGet

  PUBLIC EquationsSet_UserNumberFind

  PUBLIC EQUATIONS_SET_USER_NUMBER_FIND

CONTAINS

  !
  !================================================================================================================================
  !

  !>Gets the analytic field for an equations set if it exists.
  SUBROUTINE EquationsSet_AnalyticFieldExists(equationsSet,analyticField,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to get the analytic field for
    TYPE(FieldType), POINTER :: analyticField !<On exit, a pointer to the analytic field in the specified equations set if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsSet_AnalyticFieldExists",err,error,*998)

    IF(ASSOCIATED(analyticField)) CALL FlagError("Analytic field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    IF(ASSOCIATED(equationsSet%analytic)) THEN
      analyticField=>equationsSet%analytic%analyticField
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
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to get the analytic field for
    TYPE(FieldType), POINTER :: analyticField !<On exit, a pointer to the analytic field in the specified equations set. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_AnalyticFieldGet",err,error,*998)

    IF(ASSOCIATED(analyticField)) CALL FlagError("Analytic field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet%analytic)) THEN
      localError="Analytic information is not associated for equations set number "// &
        & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    analyticField=>equationsSet%analytic%analyticField
    IF(.NOT.ASSOCIATED(analyticField)) THEN
      localError="Analytic field is not associated for equations set number "// &
        & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))//"."
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

  !>Gets the analytic function type for an equations set.
  SUBROUTINE EquationsSet_AnalyticFunctionTypeGet(equationsSet,analyticFunctionType,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to get the analytic function type for
    INTEGER(INTG), INTENT(OUT):: analyticFunctionType !<On exit, the analytic function type in the specified equations set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_AnalyticFunctionTypeGet",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet%analytic)) THEN
      localError="Analytic information is not associated for equations set number "// &
        & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    analyticFunctionType=equationsSet%analytic%analyticFunctionType
       
    EXITS("EquationsSet_AnalyticFunctionTypeGet")
    RETURN
999 ERRORSEXITS("EquationsSet_AnalyticFunctionTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AnalyticFunctionTypeGet

  !
  !================================================================================================================================
  !

  !>Returns the analytic time for an equations set. \see OpenCMISS::cmfe_EquationsSet_AnalyticTimeGet
  SUBROUTINE EquationsSet_AnalyticTimeGet(equationsSet,analyticTime,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to get the time for.
    REAL(DP), INTENT(OUT) :: analyticTime !<On return, the analytic time value.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsSet_AnalyticTimeGet",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)    
    IF(.NOT.ASSOCIATED(equationsSet%analytic)) CALL FlagError("Equations set analytic is not associated.",err,error,*999)

    analyticTime=equationsSet%analytic%analyticTime
      
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
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the time for.
    REAL(DP), INTENT(IN) :: analyticTime !<The analytic time value to set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsSet_AnalyticTimeSet",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)    
    IF(.NOT.ASSOCIATED(equationsSet%analytic)) CALL FlagError("Equations set analytic is not associated.",err,error,*999)

    equationsSet%analytic%analyticTime=analyticTime
      
    EXITS("EquationsSet_AnalyticTimeSet")
    RETURN
999 ERRORSEXITS("EquationsSet_AnalyticTimeSet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AnalyticTimeSet

  !
  !=================================================================================================================================
  !

  !>Assert that an equations set has been finished.
  SUBROUTINE EquationsSet_AssertIsFinished(equationsSet,err,error,*)

    !Argument Variables
    TYPE(EquationsSetType), POINTER, INTENT(IN) :: equationsSet !<The equations set to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_AssertIsFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.equationsSet%equationsSetFinished) THEN
      localError="Equations set number "// &
        & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) THEN
        localError=localError//" of region number "//TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//" has not been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("EquationsSet_AssertIsFinished")
    RETURN
999 ERRORSEXITS("EquationsSet_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AssertIsFinished

  !
  !=================================================================================================================================
  !

  !>Assert that an equations set has not been finished.
  SUBROUTINE EquationsSet_AssertNotFinished(equationsSet,err,error,*)

    !Argument Variables
    TYPE(EquationsSetType), POINTER, INTENT(IN) :: equationsSet !<The equations set to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_AssertNotFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(equationsSet%equationsSetFinished) THEN
      localError="Equations set number "// &
        & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) THEN
        localError=localError//" of region number "//TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//" has already been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("EquationsSet_AssertNotFinished")
    RETURN
999 ERRORSEXITS("EquationsSet_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AssertNotFinished

  !
  !=================================================================================================================================
  !

  !>Assert that the analytic setup has been created for an equations set
  SUBROUTINE EquationsSet_AssertAnalyticIsCreated(equationsSet,err,error,*)

    !Argument Variables
    TYPE(EquationsSetType), POINTER, INTENT(IN) :: equationsSet !<The equations set to assert the analytic created status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_AssertAnalyticIsCreated",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equationsSet%analytic)) THEN
      localError="Analytic has not been created for equations set number "// &
        & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) THEN
        localError=localError//" of region number "//TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("EquationsSet_AssertAnalyticIsCreated")
    RETURN
999 ERRORSEXITS("EquationsSet_AssertAnalyticIsCreated",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AssertAnalyticIsCreated

  !
  !=================================================================================================================================
  !

  !>Assert that the analytic setup has been finished for an equations set
  SUBROUTINE EquationsSet_AssertAnalyticIsFinished(equationsSet,err,error,*)

    !Argument Variables
    TYPE(EquationsSetType), POINTER, INTENT(IN) :: equationsSet !<The equations set to assert the analytic finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_AssertAnalyticIsFinished",err,error,*999)

    CALL EquationsSet_AssertAnalyticIsCreated(equationsSet,err,error,*999)

    IF(.NOT.equationsSet%analytic%analyticFinished) THEN
      localError="The analytic setup has not been finished for equations set number "// &
        & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) THEN
        localError=localError//" of region number "//TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("EquationsSet_AssertAnalyticIsFinished")
    RETURN
999 ERRORSEXITS("EquationsSet_AssertAnalyticIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AssertAnalyticIsFinished

  !
  !=================================================================================================================================
  !

  !>Assert that the analytic setup has not been created for an equations set
  SUBROUTINE EquationsSet_AssertAnalyticNotCreated(equationsSet,err,error,*)

    !Argument Variables
    TYPE(EquationsSetType), POINTER, INTENT(IN) :: equationsSet !<The equations set to assert the analytic created status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_AssertAnalyticNotCreated",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(ASSOCIATED(equationsSet%analytic)) THEN
      localError="Analytic has alredy been created for equations set number "// &
        & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) THEN
        localError=localError//" of region number "//TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    EXITS("EquationsSet_AssertAnalyticNotCreated")
    RETURN
999 ERRORSEXITS("EquationsSet_AssertAnalyticNotCreated",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AssertAnalyticNotCreated

  !
  !=================================================================================================================================
  !

  !>Assert that the analytic setup has not been finished for an equations set
  SUBROUTINE EquationsSet_AssertAnalyticNotFinished(equationsSet,err,error,*)

    !Argument Variables
    TYPE(EquationsSetType), POINTER, INTENT(IN) :: equationsSet !<The equations set to assert the analytic finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_AssertAnalyticNotFinished",err,error,*999)

    CALL EquationsSet_AssertAnalyticIsCreated(equationsSet,err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equationsSet%analytic)) THEN
      localError="Analytic is not associated for equations set number "// &
        & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) THEN
        localError=localError//" of region number "//TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    IF(equationsSet%analytic%analyticFinished) THEN
      localError="The analytic setup has already been finished for equations set number "// &
        & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) THEN
        localError=localError//" of region number "//TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("EquationsSet_AssertAnalyticNotFinished")
    RETURN
999 ERRORSEXITS("EquationsSet_AssertAnalyticNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AssertAnalyticNotFinished

  !
  !=================================================================================================================================
  !

  !>Assert that the dependent setup has been finished for an equations set
  SUBROUTINE EquationsSet_AssertDependentIsFinished(equationsSet,err,error,*)

    !Argument Variables
    TYPE(EquationsSetType), POINTER, INTENT(IN) :: equationsSet !<The equations set to assert the dependent finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_AssertDependentIsFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    IF(.NOT.equationsSet%dependent%dependentFinished) THEN
      localError="The dependent setup has not been finished for equations set number "// &
        & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) THEN
        localError=localError//" of region number "//TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("EquationsSet_AssertDependentIsFinished")
    RETURN
999 ERRORSEXITS("EquationsSet_AssertDependentIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AssertDependentIsFinished

  !
  !=================================================================================================================================
  !

  !>Assert that the dependent setup has not been finished for an equations set
  SUBROUTINE EquationsSet_AssertDependentNotFinished(equationsSet,err,error,*)

    !Argument Variables
    TYPE(EquationsSetType), POINTER, INTENT(IN) :: equationsSet !<The equations set to assert the dependent finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_AssertDependentNotFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    IF(equationsSet%dependent%dependentFinished) THEN
      localError="The dependent setup has already been finished for equations set number "// &
        & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) THEN
        localError=localError//" of region number "//TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("EquationsSet_AssertDependentNotFinished")
    RETURN
999 ERRORSEXITS("EquationsSet_AssertDependentNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AssertDependentNotFinished

  !
  !=================================================================================================================================
  !

  !>Assert that the derived setup has been created for an equations set
  SUBROUTINE EquationsSet_AssertDerivedIsCreated(equationsSet,err,error,*)

    !Argument Variables
    TYPE(EquationsSetType), POINTER, INTENT(IN) :: equationsSet !<The equations set to assert the derived created status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_AssertDerivedIsCreated",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equationsSet%derived)) THEN
      localError="Derived has not been created for equations set number "// &
        & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) THEN
        localError=localError//" of region number "//TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    EXITS("EquationsSet_AssertDerivedIsCreated")
    RETURN
999 ERRORSEXITS("EquationsSet_AssertDerivedIsCreated",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AssertDerivedIsCreated

  !
  !=================================================================================================================================
  !

  !>Assert that the derived setup has been finished for an equations set
  SUBROUTINE EquationsSet_AssertDerivedIsFinished(equationsSet,err,error,*)

    !Argument Variables
    TYPE(EquationsSetType), POINTER, INTENT(IN) :: equationsSet !<The equations set to assert the derived finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_AssertDerivedIsFinished",err,error,*999)

    CALL EquationsSet_AssertDerivedIsCreated(equationsSet,err,error,*999)

    IF(.NOT.equationsSet%derived%derivedFinished) THEN
      localError="The derived setup has not been finished for equations set number "// &
        & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) THEN
        localError=localError//" of region number "//TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("EquationsSet_AssertDerivedIsFinished")
    RETURN
999 ERRORSEXITS("EquationsSet_AssertDerivedIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AssertDerivedIsFinished

  !
  !=================================================================================================================================
  !

  !>Assert that the derived setup has not been created for an equations set
  SUBROUTINE EquationsSet_AssertDerivedNotCreated(equationsSet,err,error,*)

    !Argument Variables
    TYPE(EquationsSetType), POINTER, INTENT(IN) :: equationsSet !<The equations set to assert the derived created status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_AssertDerivedNotCreated",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(ASSOCIATED(equationsSet%derived)) THEN
      localError="Derived has already been created for equations set number "// &
        & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) THEN
        localError=localError//" of region number "//TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    EXITS("EquationsSet_AssertDerivedNotCreated")
    RETURN
999 ERRORSEXITS("EquationsSet_AssertDerivedNotCreated",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AssertDerivedNotCreated

  !
  !=================================================================================================================================
  !

  !>Assert that the derived setup has not been finished for an equations set
  SUBROUTINE EquationsSet_AssertDerivedNotFinished(equationsSet,err,error,*)

    !Argument Variables
    TYPE(EquationsSetType), POINTER, INTENT(IN) :: equationsSet !<The equations set to assert the derived finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_AssertDerivedNotFinished",err,error,*999)

    CALL EquationsSet_AssertDerivedIsCreated(equationsSet,err,error,*999)

    IF(equationsSet%derived%derivedFinished) THEN
      localError="The derived setup has already been finished for equations set number "// &
        & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) THEN
        localError=localError//" of region number "//TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("EquationsSet_AssertDerivedNotFinished")
    RETURN
999 ERRORSEXITS("EquationsSet_AssertDerivedNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AssertDerivedNotFinished

  !
  !=================================================================================================================================
  !

  !>Assert that the equations setup has been created for an equations set
  SUBROUTINE EquationsSet_AssertEquationsIsCreated(equationsSet,err,error,*)

    !Argument Variables
    TYPE(EquationsSetType), POINTER, INTENT(IN) :: equationsSet !<The equations set to assert the equations created status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_AssertEquationsIsCreated",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equationsSet%equations)) THEN
      localError="Equations has not been created for equations set number "// &
        & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) THEN
        localError=localError//" of region number "//TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    EXITS("EquationsSet_AssertEquationsIsCreated")
    RETURN
999 ERRORSEXITS("EquationsSet_AssertEquationsIsCreated",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AssertEquationsIsCreated

  !
  !=================================================================================================================================
  !

  !>Assert that the equations setup has been finished for an equations set
  SUBROUTINE EquationsSet_AssertEquationsIsFinished(equationsSet,err,error,*)

    !Argument Variables
    TYPE(EquationsSetType), POINTER, INTENT(IN) :: equationsSet !<The equations set to assert the equations finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_AssertEquationsIsFinished",err,error,*999)

    CALL EquationsSet_AssertEquationsIsCreated(equationsSet,err,error,*999)

    IF(.NOT.equationsSet%equations%equationsFinished) THEN
      localError="The equations setup has not been finished for equations set number "// &
        & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) THEN
        localError=localError//" of region number "//TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("EquationsSet_AssertEquationsIsFinished")
    RETURN
999 ERRORSEXITS("EquationsSet_AssertEquationsIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AssertEquationsIsFinished

  !
  !=================================================================================================================================
  !

  !>Assert that the equations setup has not been created for an equations set
  SUBROUTINE EquationsSet_AssertEquationsNotCreated(equationsSet,err,error,*)

    !Argument Variables
    TYPE(EquationsSetType), POINTER, INTENT(IN) :: equationsSet !<The equations set to assert the equations created status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_AssertEquationsNotCreated",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(ASSOCIATED(equationsSet%equations)) THEN
      localError="Equations has already been created for equations set number "// &
        & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) THEN
        localError=localError//" of region number "//TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    EXITS("EquationsSet_AssertEquationsNotCreated")
    RETURN
999 ERRORSEXITS("EquationsSet_AssertEquationsNotCreated",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AssertEquationsNotCreated

  !
  !=================================================================================================================================
  !

  !>Assert that the equations setup has not been finished for an equations set
  SUBROUTINE EquationsSet_AssertEquationsNotFinished(equationsSet,err,error,*)

    !Argument Variables
    TYPE(EquationsSetType), POINTER, INTENT(IN) :: equationsSet !<The equations set to assert the equations finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_AssertEquationsNotFinished",err,error,*999)

    CALL EquationsSet_AssertEquationsIsCreated(equationsSet,err,error,*999)

    IF(equationsSet%equations%equationsFinished) THEN
      localError="The equations setup has already been finished for equations set number "// &
        & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) THEN
        localError=localError//" of region number "//TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("EquationsSet_AssertEquationsNotFinished")
    RETURN
999 ERRORSEXITS("EquationsSet_AssertEquationsNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AssertEquationsNotFinished

  !
  !=================================================================================================================================
  !

  !>Assert that the independent setup has been created for an equations set
  SUBROUTINE EquationsSet_AssertIndependentIsCreated(equationsSet,err,error,*)

    !Argument Variables
    TYPE(EquationsSetType), POINTER, INTENT(IN) :: equationsSet !<The equations set to assert the independent created status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_AssertIndependentIsCreated",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equationsSet%independent)) THEN
      localError="Independent has not been created for equations set number "// &
        & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) THEN
        localError=localError//" of region number "//TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    EXITS("EquationsSet_AssertIndependentIsCreated")
    RETURN
999 ERRORSEXITS("EquationsSet_AssertIndependentIsCreated",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AssertIndependentIsCreated

  !
  !=================================================================================================================================
  !

  !>Assert that the independent setup has been finished for an equations set
  SUBROUTINE EquationsSet_AssertIndependentIsFinished(equationsSet,err,error,*)

    !Argument Variables
    TYPE(EquationsSetType), POINTER, INTENT(IN) :: equationsSet !<The equations set to assert the independent finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_AssertIndependentIsFinished",err,error,*999)

    CALL EquationsSet_AssertIndependentIsCreated(equationsSet,err,error,*999)

    IF(.NOT.equationsSet%independent%independentFinished) THEN
      localError="The independent setup has not been finished for equations set number "// &
        & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) THEN
        localError=localError//" of region number "//TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("EquationsSet_AssertIndependentIsFinished")
    RETURN
999 ERRORSEXITS("EquationsSet_AssertIndependentIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AssertIndependentIsFinished

  !
  !=================================================================================================================================
  !

  !>Assert that the independent setup has not been created for an equations set
  SUBROUTINE EquationsSet_AssertIndependentNotCreated(equationsSet,err,error,*)

    !Argument Variables
    TYPE(EquationsSetType), POINTER, INTENT(IN) :: equationsSet !<The equations set to assert the independent created status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_AssertIndependentNotCreated",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(ASSOCIATED(equationsSet%independent)) THEN
      localError="Independent has already been created for equations set number "// &
        & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) THEN
        localError=localError//" of region number "//TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    EXITS("EquationsSet_AssertIndependentNotCreated")
    RETURN
999 ERRORSEXITS("EquationsSet_AssertIndependentNotCreated",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AssertIndependentNotCreated

  !
  !=================================================================================================================================
  !

  !>Assert that the independent setup has not been finished for an equations set
  SUBROUTINE EquationsSet_AssertIndependentNotFinished(equationsSet,err,error,*)

    !Argument Variables
    TYPE(EquationsSetType), POINTER, INTENT(IN) :: equationsSet !<The equations set to assert the independent finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_AssertIndependentNotFinished",err,error,*999)

    CALL EquationsSet_AssertIndependentIsCreated(equationsSet,err,error,*999)

    IF(equationsSet%independent%independentFinished) THEN
      localError="The independent setup has already been finished for equations set number "// &
        & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) THEN
        localError=localError//" of region number "//TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("EquationsSet_AssertIndependentNotFinished")
    RETURN
999 ERRORSEXITS("EquationsSet_AssertIndependentNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AssertIndependentNotFinished

  !
  !=================================================================================================================================
  !

  !>Assert that the materials setup has been created for an equations set
  SUBROUTINE EquationsSet_AssertMaterialsIsCreated(equationsSet,err,error,*)

    !Argument Variables
    TYPE(EquationsSetType), POINTER, INTENT(IN) :: equationsSet !<The equations set to assert the materials created status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_AssertMaterialsIsCreated",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equationsSet%materials)) THEN
      localError="Materials has not been created for equations set number "// &
        & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) THEN
        localError=localError//" of region number "//TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("EquationsSet_AssertMaterialsIsCreated")
    RETURN
999 ERRORSEXITS("EquationsSet_AssertMaterialsIsCreated",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AssertMaterialsIsCreated

  !
  !=================================================================================================================================
  !

  !>Assert that the materials setup has been finished for an equations set
  SUBROUTINE EquationsSet_AssertMaterialsIsFinished(equationsSet,err,error,*)

    !Argument Variables
    TYPE(EquationsSetType), POINTER, INTENT(IN) :: equationsSet !<The equations set to assert the materials finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_AssertMaterialsIsFinished",err,error,*999)

    CALL EquationsSet_AssertMaterialsIsCreated(equationsSet,err,error,*999)

    IF(.NOT.equationsSet%materials%materialsFinished) THEN
      localError="The materials setup has not been finished for equations set number "// &
        & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) THEN
        localError=localError//" of region number "//TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("EquationsSet_AssertMaterialsIsFinished")
    RETURN
999 ERRORSEXITS("EquationsSet_AssertMaterialsIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AssertMaterialsIsFinished

  !
  !=================================================================================================================================
  !

  !>Assert that the materials setup has not been created for an equations set
  SUBROUTINE EquationsSet_AssertMaterialsNotCreated(equationsSet,err,error,*)

    !Argument Variables
    TYPE(EquationsSetType), POINTER, INTENT(IN) :: equationsSet !<The equations set to assert the materials created status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_AssertMaterialsNotCreated",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(ASSOCIATED(equationsSet%materials)) THEN
      localError="Materials has already been created for equations set number "// &
        & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) THEN
        localError=localError//" of region number "//TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    EXITS("EquationsSet_AssertMaterialsNotCreated")
    RETURN
999 ERRORSEXITS("EquationsSet_AssertMaterialsNotCreated",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AssertMaterialsNotCreated

  !
  !=================================================================================================================================
  !

  !>Assert that the materials setup has not been finished for an equations set
  SUBROUTINE EquationsSet_AssertMaterialsNotFinished(equationsSet,err,error,*)

    !Argument Variables
    TYPE(EquationsSetType), POINTER, INTENT(IN) :: equationsSet !<The equations set to assert the materials finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_AssertMaterialsNotFinished",err,error,*999)

    CALL EquationsSet_AssertMaterialsIsCreated(equationsSet,err,error,*999)

    IF(equationsSet%materials%materialsFinished) THEN
      localError="The materials setup has already been finished for equations set number "// &
        & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) THEN
        localError=localError//" of region number "//TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("EquationsSet_AssertMaterialsNotFinished")
    RETURN
999 ERRORSEXITS("EquationsSet_AssertMaterialsNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AssertMaterialsNotFinished

  !
  !=================================================================================================================================
  !

  !>Assert that the source setup has been created for an equations set
  SUBROUTINE EquationsSet_AssertSourceIsCreated(equationsSet,err,error,*)

    !Argument Variables
    TYPE(EquationsSetType), POINTER, INTENT(IN) :: equationsSet !<The equations set to assert the source created status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_AssertSourceIsCreated",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equationsSet%source)) THEN
      localError="Source has not been created for equations set number "// &
        & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) THEN
        localError=localError//" of region number "//TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
   
    EXITS("EquationsSet_AssertSourceIsCreated")
    RETURN
999 ERRORSEXITS("EquationsSet_AssertSourceIsCreated",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AssertSourceIsCreated

  !
  !=================================================================================================================================
  !

  !>Assert that the source setup has been finished for an equations set
  SUBROUTINE EquationsSet_AssertSourceIsFinished(equationsSet,err,error,*)

    !Argument Variables
    TYPE(EquationsSetType), POINTER, INTENT(IN) :: equationsSet !<The equations set to assert the source finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_AssertSourceIsFinished",err,error,*999)

    CALL EquationsSet_AssertSourceIsCreated(equationsSet,err,error,*999)

    IF(.NOT.equationsSet%source%sourceFinished) THEN
      localError="The source setup has not been finished for equations set number "// &
        & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) THEN
        localError=localError//" of region number "//TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("EquationsSet_AssertSourceIsFinished")
    RETURN
999 ERRORSEXITS("EquationsSet_AssertSourceIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AssertSourceIsFinished

  !
  !=================================================================================================================================
  !

  !>Assert that the source has not been created for an equations set
  SUBROUTINE EquationsSet_AssertSourceNotCreated(equationsSet,err,error,*)

    !Argument Variables
    TYPE(EquationsSetType), POINTER, INTENT(IN) :: equationsSet !<The equations set to assert the source created status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_AssertSourceNotCreated",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(ASSOCIATED(equationsSet%source)) THEN
      localError="Source has already been created for equations set number "// &
        & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) THEN
        localError=localError//" of region number "//TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    EXITS("EquationsSet_AssertSourceNotCreated")
    RETURN
999 ERRORSEXITS("EquationsSet_AssertSourceNotCreated",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AssertSourceNotCreated

  !
  !=================================================================================================================================
  !

  !>Assert that the source setup has not been finished for an equations set
  SUBROUTINE EquationsSet_AssertSourceNotFinished(equationsSet,err,error,*)

    !Argument Variables
    TYPE(EquationsSetType), POINTER, INTENT(IN) :: equationsSet !<The equations set to assert the source finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_AssertSourceNotFinished",err,error,*999)

    CALL EquationsSet_AssertSourceIsCreated(equationsSet,err,error,*999)

    IF(equationsSet%source%sourceFinished) THEN
      localError="The source setup has already been finished for equations set number "// &
        & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) THEN
        localError=localError//" of region number "//TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      ENDIF
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("EquationsSet_AssertSourceNotFinished")
    RETURN
999 ERRORSEXITS("EquationsSet_AssertSourceNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AssertSourceNotFinished

  !
  !================================================================================================================================
  !

  !>Gets the coordinate system used for an equations set.
  SUBROUTINE EquationsSet_CoordinateSystemGet(equationsSet,coordinateSystem,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to get the coordinate system for
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem !<On exit, a pointer to the coordinate system for the specified equations set. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(RegionType), POINTER :: region
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_CoordinateSystemGet",err,error,*998)

    IF(ASSOCIATED(coordinateSystem)) CALL FlagError("Coordinate system is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    region=>equationsSet%region
    IF(.NOT.ASSOCIATED(region)) THEN
      localError="Region is not associated for equations set number "// &
      & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    coordinateSystem=>region%coordinateSystem
    IF(.NOT.ASSOCIATED(coordinateSystem)) THEN
      localError="Coordinate system is not associated for region number "// &
        & TRIM(NumberToVString(region%userNumber,"*",err,error))//" of equations set number "// &
        & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))//"."
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

  !>Gets the derived field for an equations set if it exists.
  SUBROUTINE EquationsSet_DerivedFieldExists(equationsSet,derivedField,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to get the derived field for
    TYPE(FieldType), POINTER :: derivedField !<On exit, a pointer to the derived field in the specified equations set if it exists. Must not be associated on entry
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
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to get the derived field for
    TYPE(FieldType), POINTER :: derivedField !<On exit, a pointer to the derived field in the specified equations set. Must not be associated on entry
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
      & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) &
        & localError=localError//" of region number "// &
        & TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      localError=localError//"."
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
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to get the dependent field for
    TYPE(FieldType), POINTER :: dependentField !<On exit, a pointer to the dependent field in the specified equations set. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_DependentFieldGet",err,error,*998)

    IF(ASSOCIATED(dependentField)) CALL FlagError("Dependent field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    dependentField=>equationsSet%dependent%dependentField
    IF(.NOT.ASSOCIATED(dependentField)) THEN
      localError="Dependent field is not associated for equations set number "// &
      & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) &
        & localError=localError//" of region number "// &
        & TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      localError=localError//"."
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
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to get the equations for
    TYPE(EquationsType), POINTER :: equations !<On exit, a pointer to the equations in the specified equations set. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_EquationsGet",err,error,*998)

    IF(ASSOCIATED(equations)) CALL FlagError("Equations is already associated.",err,error,*998)
    CALL EquationsSet_AssertIsFinished(equationsSet,err,error,*999)

    equations=>equationsSet%equations
    IF(.NOT.ASSOCIATED(equations)) THEN
      localError="Equations is not associated for equations set number "// &
      & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) &
        & localError=localError//" of region number "// &
        & TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      localError=localError//"."
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

  !>Gets the equation sets for an equations set.
  SUBROUTINE EquationsSet_EquationsSetsGet(equationsSet,equationsSets,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to get the equations sets for
    TYPE(EquationsSetsType), POINTER :: equationsSets !<On exit, a pointer to the equations sets for the specified equations set. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_EquationsSetsGet",err,error,*998)

    IF(ASSOCIATED(equationsSets)) CALL FlagError("Equations sets is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
 
    equationsSets=>equationsSet%equationsSets
    IF(.NOT.ASSOCIATED(equationsSets)) THEN
      localError="Equations sets is not associated for equations set number "// &
        & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) &
        & localError=localError//" of region number "// &
        & TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("EquationsSet_EquationsSetsGet")
    RETURN
999 NULLIFY(equationsSets)
998 ERRORSEXITS("EquationsSet_EquationsSetsGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_EquationsSetsGet

  !
  !================================================================================================================================
  !

  !>Gets the equations set field for an equations set.
  SUBROUTINE EquationsSet_EquationsSetFieldGet(equationsSet,equationsSetField,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to get the equations set field for
    TYPE(FieldType), POINTER :: equationsSetField !<On exit, a pointer to the equations set field in the specified equations set. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_EquationsSetFieldGet",err,error,*998)

    IF(ASSOCIATED(equationsSetField)) CALL FlagError("Equations set field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    equationsSetField=>equationsSet%equationsSetField%equationsSetFieldField
    IF(.NOT.ASSOCIATED(equationsSetField)) THEN
      localError="Equations set field is not associated for equations set number "// &
      & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) &
        & localError=localError//" of region number "// &
        & TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      localError=localError//"."
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
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to get the fibre field for
    TYPE(FieldType), POINTER :: fibreField !<On exit, a pointer to the fibre field in the specified equations set if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsSet_FibreFieldExists",err,error,*998)

    IF(ASSOCIATED(fibreField)) CALL FlagError("Fibre field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    fibreField=>equationsSet%geometry%fibreField
      
    EXITS("EquationsSet_FibreFieldExists")
    RETURN
999 NULLIFY(fibreField)
998 ERRORSEXITS("EquationsSet_FibreFieldExists",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_FibreFieldExists

  !
  !================================================================================================================================
  !

  !>Gets the fibre field for an equations set.
  SUBROUTINE EquationsSet_FibreFieldGet(equationsSet,fibreField,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to get the fibre field for
    TYPE(FieldType), POINTER :: fibreField !<On exit, a pointer to the fibre field in the specified equations set. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_FibreFieldGet",err,error,*998)

    IF(ASSOCIATED(fibreField)) CALL FlagError("Fibre field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    fibreField=>equationsSet%geometry%fibreField
    IF(.NOT.ASSOCIATED(fibreField)) THEN
      localError="Fibre field is not associated for equations set number "// &
        & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) &
        & localError=localError//" of region number "// &
        & TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
      
    EXITS("EquationsSet_FibreFieldGet")
    RETURN
999 NULLIFY(fibreField)
998 ERRORSEXITS("EquationsSet_FibreFieldGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_FibreFieldGet

  !
  !================================================================================================================================
  !

  !>Gets the geometric field for an equations set.
  SUBROUTINE EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to get the geometric field for
    TYPE(FieldType), POINTER :: geometricField !<On exit, a pointer to the geometric field in the specified equations set. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_GeometricFieldGet",err,error,*998)

    IF(ASSOCIATED(geometricField)) CALL FlagError("Geometric field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    geometricField=>equationsSet%geometry%geometricField
    IF(.NOT.ASSOCIATED(geometricField)) THEN
      localError="Geometric field is not associated for equations set number "// &
        & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) &
        & localError=localError//" of region number "// &
        & TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      localError=localError//"."
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

  !>Gets the independent field for an equations set if it exists.
  SUBROUTINE EquationsSet_IndependentFieldExists(equationsSet,independentField,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to get the independent field for
    TYPE(FieldType), POINTER :: independentField !<On exit, a pointer to the independent field in the specified equations set if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsSet_IndependentFieldExists",err,error,*998)

    IF(ASSOCIATED(independentField)) CALL FlagError("Independent field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    IF(ASSOCIATED(equationsSet%independent)) THEN
      independentField=>equationsSet%INDEPENDENT%independentField
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
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to get the independent field for
    TYPE(FieldType), POINTER :: independentField !<On exit, a pointer to the independent field in the specified equations set. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_IndependentFieldGet",err,error,*998)

    IF(ASSOCIATED(independentField)) CALL FlagError("Independent field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet%independent)) THEN
      localError="Independent information is not associated for equations set number "// &
        & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) &
        & localError=localError//" of region number "// &
        & TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    independentField=>equationsSet%independent%independentField
    IF(.NOT.ASSOCIATED(independentField)) THEN
      localError="Independent field is not associated for equations set number "// &
      & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) &
        & localError=localError//" of region number "// &
        & TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      localError=localError//"."
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
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to get the label for
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
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to get the label for
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
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the label for 
    CHARACTER(LEN=*), INTENT(IN) :: label !<The label to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsSet_LabelSetC",err,error,*999)

    CALL EquationsSet_AssertNotFinished(equationsSet,err,error,*999)
   
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
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the label for 
    TYPE(VARYING_STRING), INTENT(IN) :: label !<The label to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsSet_LabelSetVS",err,error,*999)

    CALL EquationsSet_AssertNotFinished(equationsSet,err,error,*999)
    
    equationsSet%label=label
    
    EXITS("EquationsSet_LabelSetVS")
    RETURN
999 ERRORSEXITS("EquationsSet_LabelSetVS",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_LabelSetVS

  !
  !================================================================================================================================
  !

  !>Gets the materials field for an equations set if it exists.
  SUBROUTINE EquationsSet_MaterialsFieldExists(equationsSet,materialsField,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to get the materials field for
    TYPE(FieldType), POINTER :: materialsField !<On exit, a pointer to the materials field in the specified equations set if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsSet_MaterialsFieldExists",err,error,*998)

    IF(ASSOCIATED(materialsField)) CALL FlagError("Materials field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    IF(ASSOCIATED(equationsSet%materials)) THEN
      materialsField=>equationsSet%materials%materialsField
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
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to get the materials field for
    TYPE(FieldType), POINTER :: materialsField !<On exit, a pointer to the materials field in the specified equations set. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_MaterialsFieldGet",err,error,*998)

    IF(ASSOCIATED(materialsField)) CALL FlagError("Materials field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet%materials)) THEN
      localError="Materials information is not associated for equations set number "// &
        & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) &
        & localError=localError//" of region number "// &
        & TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    materialsField=>equationsSet%materials%materialsField
    IF(.NOT.ASSOCIATED(materialsField)) THEN
      localError="Materials field is not associated for equations set number "// &
      & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) &
        & localError=localError//" of region number "// &
        & TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      localError=localError//"."
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

  !>Gets the output type for an equations set.
  SUBROUTINE EquationsSet_OutputTypeGet(equationsSet,outputType,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to get the output type for
    INTEGER(INTG), INTENT(OUT) :: outputType !<On exit, the output type of the equations set \see EquationsSetConstants_OutputTypes,EquationsSetConstants
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsSet_OutputTypeGet",err,error,*999)

    CALL EquationsSet_AssertIsFinished(equationsSet,err,error,*999)
   
    outputType=equationsSet%outputType
       
    EXITS("EquationsSet_OutputTypeGet")
    RETURN
999 ERRORSEXITS("EquationsSet_OutputTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_OutputTypeGet
  
  !
  !================================================================================================================================
  !

  !>Gets the region used for an equations set.
  SUBROUTINE EquationsSet_RegionGet(equationsSet,region,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to get the region for
    TYPE(RegionType), POINTER :: region !<On exit, a pointer to the region for the specified equations set. Must not be associated on entry
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
      & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))//"."
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

  !>Returns the solution method for an equations set. \see OpenCMISS::Iron::cmfe_EquationsSet_SolutionMethodGet
  SUBROUTINE EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to get the solution method for
    INTEGER(INTG), INTENT(OUT) :: solutionMethod !<On return, the equations set solution method \see EquationsSetConstants_SolutionMethods,EquationsSetConstants
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsSet_SolutionMethodGet",err,error,*999)

    CALL EquationsSet_AssertIsFinished(equationsSet,err,error,*999)
    
    solutionMethod=equationsSet%solutionMethod
     
    EXITS("EquationsSet_SolutionMethodGet")
    RETURN
999 ERRORSEXITS("EquationsSet_SolutionMethodGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_SolutionMethodGet
  
  !
  !================================================================================================================================
  !

  !>Gets the source field for an equations set if it exists
  SUBROUTINE EquationsSet_SourceFieldExists(equationsSet,sourceField,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to get the source field for
    TYPE(FieldType), POINTER :: sourceField !<On exit, a pointer to the source field in the specified equations set if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsSet_SourceFieldExists",err,error,*998)

    IF(ASSOCIATED(sourceField)) CALL FlagError("Source field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    IF(ASSOCIATED(equationsSet%source)) THEN
      sourceField=>equationsSet%source%sourceField
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
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to get the source field for
    TYPE(FieldType), POINTER :: sourceField !<On exit, a pointer to the source field in the specified equations set. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_SourceFieldGet",err,error,*998)

    IF(ASSOCIATED(sourceField)) CALL FlagError("Source field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet%source)) THEN
      localError="Source information is not associated for equations set number "// &
        & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) &
        & localError=localError//" of region number "// &
        & TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      localError=localError//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    sourceField=>equationsSet%source%sourceField
    IF(.NOT.ASSOCIATED(sourceField)) THEN
      localError="Source field is not associated for equations set number "// &
      & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))
      IF(ASSOCIATED(equationsSet%region)) &
        & localError=localError//" of region number "// &
        & TRIM(NumberToVString(equationsSet%region%userNumber,"*",err,error))
      localError=localError//"."
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

  !>Returns the equations set specification i.e., equations set class, type and subtype for an equations set. \see OpenCMISS::Iron::cmfe_EquationsSet_SpecificationGet
  SUBROUTINE EquationsSet_SpecificationGet(equationsSet,equationsSetSpecification,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to get the specification for
    INTEGER(INTG), INTENT(INOUT) :: equationsSetSpecification(:) !<On return, The equations set specifcation array. Must be allocated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: specificationLength,specificationIdx
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsSet_SpecificationGet",err,error,*999)
   
    CALL EquationsSet_AssertIsFinished(equationsSet,err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) &
      & CALL FlagError("Equations set specification has not been allocated.",err,error,*999)
    specificationLength=0
    DO specificationIdx=1,SIZE(equationsSet%specification,1)
      IF(equationsSet%specification(specificationIdx)>0) THEN
        specificationLength=specificationIdx
      ENDIF
    ENDDO !specificationIdx
    IF(SIZE(equationsSetSpecification,1)<specificationLength) THEN
      localError="The equations set specification array size is "//TRIM(NumberToVstring(specificationLength,"*",err,error))// &
        & " and it needs to be >= "//TRIM(NumberToVstring(SIZE(equationsSetSpecification,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    equationsSetSpecification(1:specificationLength)=equationsSet%specification(1:specificationLength)

    EXITS("EquationsSet_SpecificationGet")
    RETURN
999 ERRORSEXITS("EquationsSet_SpecificationGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_SpecificationGet

  !
  !================================================================================================================================
  !

  !>Gets the size of the equations set specification array for a problem identified by a pointer. \see OpenCMISS::Iron::cmfe_EquationsSet_SpecificationSizeGet
  SUBROUTINE EquationsSet_SpecificationSizeGet(equationsSet,specificationSize,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations Set to get the specification for.
    INTEGER(INTG), INTENT(OUT) :: specificationSize !<On return, the size of the problem specifcation array.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsSet_SpecificationSizeGet",err,error,*999)

    specificationSize=0
    CALL EquationsSet_AssertIsFinished(equationsSet,err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    
    specificationSize=SIZE(equationsSet%specification,1)

    EXITS("EquationsSet_SpecificationSizeGet")
    RETURN
999 ERRORSEXITS("EquationsSet_SpecificationSizeGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_SpecificationSizeGet

  !
  !================================================================================================================================
  !

  !>Gets the current times for an equations set. \see OpenCMISS::Iron::cmfe_EquationsSet_TimesGet
  SUBROUTINE EquationsSet_TimesGet(equationsSet,currentTime,deltaTime,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to get the times for
    REAL(DP), INTENT(OUT) :: currentTime !<The current time for the equations set to get.
    REAL(DP), INTENT(OUT) :: deltaTime !<The current time incremenet for the equations set to get
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsSet_TimesGet",err,error,*999) 

    CALL EquationsSet_AssertIsFinished(equationsSet,err,error,*999)

    currentTime=equationsSet%currentTime
    deltaTime=equationsSet%deltaTime
      
    EXITS("EquationsSet_TimesGet")
    RETURN
999 ERRORSEXITS("EquationsSet_TimesGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_TimesGet
  
  !
  !================================================================================================================================
  !

  !>Finds and returns a pointer to the equations set identified by user number in the given region. If no equations set with that user number exists the equations set is left nullified.
  SUBROUTINE EquationsSet_UserNumberFind(userNumber,region,equationsSet,err,error,*)

    !Argument variables 
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the equation set to find.
    TYPE(RegionType), POINTER :: region !<The region to find the equations set in
    TYPE(EquationsSetType), POINTER :: equationsSet !<On return, a pointer to the equations set if an equations set with the specified user number exists in the given region. If no equation set with the specified number exists a NULL pointer is returned. The pointer must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsSet_UserNumberFind",err,error,*999)

    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated.",err,error,*999)
    IF(ASSOCIATED(equationsSet)) CALL FlagError("Equations set is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(region%equationsSets)) THEN
      localError="The equations sets on region number "//TRIM(NumberToVString(region%userNumber,"*",err,error))// &
        & " are not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    NULLIFY(equationsSet)
    IF(ASSOCIATED(region%equationsSets%equationsSets)) THEN
      DO equationsSetIdx=1,region%equationsSets%numberOfEquationsSets
        IF(ASSOCIATED(region%equationsSets%equationsSets(equationsSetIdx)%ptr)) THEN
          IF(region%equationsSets%equationsSets(equationsSetIdx)%ptr%userNumber==userNumber) THEN
            equationsSet=>region%equationsSets%equationsSets(equationsSetIdx)%ptr
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
