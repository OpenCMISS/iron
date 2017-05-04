!> \file
!> \author Chris Bradley
!> \brief This module contains all equations access method routines.
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

!> This module contains all equations access method routines.
MODULE EquationsAccessRoutines
  
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

  PUBLIC Equations_EquationsSetGet
  
  PUBLIC EquationsScalar_EquationsGet,EquationsVector_EquationsGet
  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Gets the equations set for an equations.
  SUBROUTINE Equations_EquationsSetGet(equations,equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER :: equations !<A pointer to the equations to get the equations set for
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<On exit, a pointer to the equations set for the specified equations. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Equations_EquationsSetGet",err,error,*998)

    IF(ASSOCIATED(equationsSet)) CALL FlagError("Equations set is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)

    equationsSet=>equations%equationsSet
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated for the equations.",err,error,*999)
       
    EXITS("Equations_EquationsSetGet")
    RETURN
999 NULLIFY(equationsSet)
998 ERRORSEXITS("Equations_EquationsSetGet",err,error)
    RETURN 1
    
  END SUBROUTINE Equations_EquationsSetGet

  !
  !================================================================================================================================
  !

  !>Gets the equations for a scalar equations.
  SUBROUTINE EquationsScalar_EquationsGet(scalarEquations,equations,err,error,*)

    !Argument variables
    TYPE(EquationsScalarType), POINTER :: scalarEquations !<A pointer to the scalar equations to get the equations for
    TYPE(EquationsType), POINTER :: equations !<On exit, a pointer to the equations in the specified scalar equations. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsScalar_EquationsGet",err,error,*998)

    IF(ASSOCIATED(equations)) CALL FlagError("Equations is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(scalarEquations)) CALL FlagError("Scalar equations is not associated.",err,error,*999)

    equations=>scalarEquations%equations
    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated for the scalar equations.",err,error,*999)
       
    EXITS("EquationsScalar_EquationsGet")
    RETURN
999 NULLIFY(equations)
998 ERRORSEXITS("EquationsScalar_EquationsGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsScalar_EquationsGet
  
  !
  !================================================================================================================================
  !

  !>Gets the equations for a vector equations.
  SUBROUTINE EquationsVector_EquationsGet(vectorEquations,equations,err,error,*)

    !Argument variables
    TYPE(EquationsVectorType), POINTER :: vectorEquations !<A pointer to the vector equations to get the equations for
    TYPE(EquationsType), POINTER :: equations !<On exit, a pointer to the equations in the specified vector equations. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsVector_EquationsGet",err,error,*998)

    IF(ASSOCIATED(equations)) CALL FlagError("Equations is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorEquations)) CALL FlagError("Vector equations is not associated.",err,error,*999)

    equations=>vectorEquations%equations
    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated for the vector equations.",err,error,*999)
       
    EXITS("EquationsVector_EquationsGet")
    RETURN
999 NULLIFY(equations)
998 ERRORSEXITS("EquationsVector_EquationsGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsVector_EquationsGet

  !
  !================================================================================================================================
  !

END MODULE EquationsAccessRoutines
