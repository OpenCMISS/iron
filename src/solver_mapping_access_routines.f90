!> \file
!> \author Chris Bradley
!> \brief This module contains all solver mapping access method routines.
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

!> This module contains all solver mapping access method routines.
MODULE SolverMappingAccessRoutines
  
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

  PUBLIC SolverMapping_EquationsSetGet

  PUBLIC SolverMapping_InterfaceConditionGet
  
CONTAINS

  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the equations set for solver mapping.
  SUBROUTINE SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*)

    !Argument variables
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping !<A pointer to the solver mapping to get the equations set for
    INTEGER(INTG), INTENT(IN) :: equationsSetIdx !<The equations set index in the solver mapping to get the equations set for
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<On exit, a pointer to the specified equations set. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("SolverMapping_EquationsSetGet",err,error,*998)

    IF(ASSOCIATED(equationsSet)) CALL FlagError("Equations set is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMapping)) CALL FlagError("Solver mapping is not associated.",err,error,*999)
    IF(equationsSetIdx<0.OR.equationsSetIdx>solverMapping%NUMBER_OF_EQUATIONS_SETS) THEN
      localError="The specified equations set index of "//TRIM(NumberToVString(equationsSetIdx,"*",err,error))// &
        & " is invalid. The index must be > 0 and <= "// &
          & TRIM(NumberToVString(SolverMapping%NUMBER_OF_EQUATIONS_SETS,"*",err,error))//"."      
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(solverMapping%EQUATIONS_SETS)) &
      & CALL FlagError("Solver mapping equations sets is not allocated.",err,error,*999)

    equationsSet=>solverMapping%EQUATIONS_SETS(equationsSetIdx)%ptr
    IF(.NOT.ASSOCIATED(equationsSet)) THEN
      localError="The equations set for the specified equations set index of "// &
        & TRIM(NumberToVString(equationsSetIdx,"*",err,error))//" is not associated."      
      CALL FlagError(localError,err,error,*999)
    ENDIF
      
    EXITS("SolverMapping_EquationsSetGet")
    RETURN
999 NULLIFY(equationsSet)
998 ERRORSEXITS("SolverMapping_EquationsSetGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_EquationsSetGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the interface condition for solver mapping.
  SUBROUTINE SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*)

    !Argument variables
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping !<A pointer to the solver mapping to get the interface condition for
    INTEGER(INTG), INTENT(IN) :: interfaceConditionIdx !<The interface condition index in the solver mapping to get the interface condition for
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: interfaceCondition !<On exit, a pointer to the specified interface condition. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("SolverMapping_InterfaceConditionGet",err,error,*998)

    IF(ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMapping)) CALL FlagError("Solver mapping is not associated.",err,error,*999)
    IF(interfaceConditionIdx<1.OR.interfaceConditionIdx>solverMapping%NUMBER_OF_INTERFACE_CONDITIONS) THEN
      localError="The specified interface condition index of "//TRIM(NumberToVString(interfaceConditionIdx,"*",err,error))// &
        & " is invalid. The index must be > 1 and <= "// &
        & TRIM(NumberToVString(SolverMapping%NUMBER_OF_INTERFACE_CONDITIONS,"*",err,error))//"."      
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(solverMapping%INTERFACE_CONDITIONS)) &
      & CALL FlagError("Solver mapping interface conditions is not allocated.",err,error,*999)

    interfaceCondition=>solverMapping%INTERFACE_CONDITIONS(interfaceConditionIdx)%ptr
    IF(.NOT.ASSOCIATED(interfaceCondition)) THEN
      localError="The interface condition for the specified interface condition index of "// &
        & TRIM(NumberToVString(interfaceConditionIdx,"*",err,error))//" is not associated."      
      CALL FlagError(localError,err,error,*999)
    ENDIF
      
    EXITS("SolverMapping_InterfaceConditionGet")
    RETURN
999 NULLIFY(interfaceCondition)
998 ERRORSEXITS("SolverMapping_InterfaceConditionGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_InterfaceConditionGet
  
  !
  !================================================================================================================================
  !

END MODULE SolverMappingAccessRoutines
