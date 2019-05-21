!> \file
!> \author Chris Bradley
!> \brief This module contains all interface mapping access method routines.
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
MODULE InterfaceMappingAccessRoutines
  
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

  PUBLIC InterfaceMapping_AssertIsFinished,InterfaceMapping_AssertNotFinished

  PUBLIC InterfaceMapping_InterfaceEquationsGet

  PUBLIC InterfaceMapping_LagrangeVariableGet

  PUBLIC InterfaceMapping_MatrixVariableGet

  PUBLIC InterfaceMapping_RHSMappingGet

  CONTAINS

  !
  !================================================================================================================================
  !

  !>Assert that an interface mapping has been finished
  SUBROUTINE InterfaceMapping_AssertIsFinished(interfaceMapping,err,error,*)

    !Argument Variables
    TYPE(InterfaceMappingType), POINTER, INTENT(IN) :: interfaceMapping !<The interface mapping to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMapping_AssertIsFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceMapping)) CALL FlagError("Interface mapping is not associated.",err,error,*999)

    IF(.NOT.interfaceMapping%interfaceMappingFinished) CALL FlagError("Interface mapping has not been finished.",err,error,*999)
    
    EXITS("InterfaceMapping_AssertIsFinished")
    RETURN
999 ERRORSEXITS("InterfaceMapping_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_AssertIsFinished

  !
  !================================================================================================================================
  !

  !>Assert that an interface mapping has not been finished
  SUBROUTINE InterfaceMapping_AssertNotFinished(interfaceMapping,err,error,*)

    !Argument Variables
    TYPE(InterfaceMappingType), POINTER, INTENT(IN) :: interfaceMapping !<The interface mapping to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMapping_AssertNotFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceMapping)) CALL FlagError("Interface mapping is not associated.",err,error,*999)

    IF(interfaceMapping%interfaceMappingFinished) CALL FlagError("Interface mapping has already been finished.",err,error,*999)
    
    EXITS("InterfaceMapping_AssertNotFinished")
    RETURN
999 ERRORSEXITS("InterfaceMapping_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_AssertNotFinished

  !
  !================================================================================================================================
  !

  !>Gets the interface equations for an interface mapping.
  SUBROUTINE InterfaceMapping_InterfaceEquationsGet(interfaceMapping,interfaceEquations,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping !<A pointer to the interface mapping to get the interface equations for
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: interfaceEquations !<On exit, a pointer to the interface equations in the specified interface mapping. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMapping_InterfaceEquationsGet",err,error,*998)

    IF(ASSOCIATED(interfaceEquations)) CALL FlagError("Interface equations is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMapping)) CALL FlagError("Interface mapping is not associated.",err,error,*999)

    interfaceEquations=>interfaceMapping%interfaceEquations
    IF(.NOT.ASSOCIATED(interfaceEquations)) &
      & CALL FlagError("Interface mapping interface equations is not associated.",err,error,*999)
       
    EXITS("InterfaceMapping_InterfaceEquationsGet")
    RETURN
999 NULLIFY(interfaceEquations)
998 ERRORSEXITS("InterfaceMapping_InterfaceEquationsGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_InterfaceEquationsGet

  !
  !================================================================================================================================
  !

  !>Gets the lagrange variable for an interface mapping.
  SUBROUTINE InterfaceMapping_LagrangeVariableGet(interfaceMapping,lagrangeVariable,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping !<A pointer to the interface mapping to get the lagrange variable for
    TYPE(FieldVariableType), POINTER :: lagrangeVariable !<On exit, a pointer to the lagrange variable in the specified interface mapping. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMapping_LagrangeVariableGet",err,error,*998)

    IF(ASSOCIATED(lagrangeVariable)) CALL FlagError("Lagrange variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMapping)) CALL FlagError("Interface mapping is not associated.",err,error,*999)

    lagrangeVariable=>interfaceMapping%lagrangeVariable
    IF(.NOT.ASSOCIATED(lagrangeVariable)) &
      & CALL FlagError("Interface mapping Lagrange variable is not associated.",err,error,*999)
       
    EXITS("InterfaceMapping_LagrangeVariableGet")
    RETURN
999 NULLIFY(lagrangeVariable)
998 ERRORSEXITS("InterfaceMapping_LagrangeVariableGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_LagrangeVariableGet

  !
  !================================================================================================================================
  !

  !>Gets the matrix (row) variable for an interface matrix in an interface mapping.
  SUBROUTINE InterfaceMapping_MatrixVariableGet(interfaceMapping,matrixNumber,matrixVariable,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping !<A pointer to the interface mapping to get the matrix variable for
    INTEGER(INTG), INTENT(IN) :: matrixNumber !<The number of the matrix to get the variable for
    TYPE(FieldVariableType), POINTER :: matrixVariable !<On exit, a pointer to the field variable for the specified interface matrix. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("InterfaceMapping_MatrixVariableGet",err,error,*998)

    IF(ASSOCIATED(matrixVariable)) CALL FlagError("Matrix variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMapping)) CALL FlagError("Interface mapping is not associated.",err,error,*999)
    IF(matrixNumber<=0.OR.matrixNumber>interfaceMapping%numberOfInterfaceMatrices) THEN
      localError="The specified matrix number of "//TRIM(NumberToVString(matrixNumber,"*",err,error))// &
        & " is invalid. The matrix number must be >= 1 and <= "// &
        & TRIM(NumberToVString(interfaceMapping%numberOfInterfaceMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(interfaceMapping%interfaceMatrixRowsToVarMaps)) &
      & CALL FlagError("The interface matrix row to variable maps is not allocated for the interface mapping.",err,error,*999)

    matrixVariable=>interfaceMapping%interfaceMatrixRowsToVarMaps(matrixNumber)%variable
    IF(.NOT.ASSOCIATED(matrixVariable)) THEN
      localError="The interface mapping matrix variable is not associted for matrix number "// &
        & TRIM(NumberToVString(matrixNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("InterfaceMapping_MatrixVariableGet")
    RETURN
999 NULLIFY(matrixVariable)
998 ERRORSEXITS("InterfaceMapping_MatrixVariableGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_MatrixVariableGet

  !
  !================================================================================================================================
  !

  !>Gets the RHS mapping for an interface mapping.
  SUBROUTINE InterfaceMapping_RHSMappingGet(interfaceMapping,rhsMapping,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping !<A pointer to the interface mapping to get the RHS mapping for
    TYPE(InterfaceMappingRHSType), POINTER :: rhsMapping !<On exit, a pointer to the RHS mapping in the specified interface mapping. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMapping_RHSMappingGet",err,error,*998)

    IF(ASSOCIATED(rhsMapping)) CALL FlagError("RHS mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMapping)) CALL FlagError("Interface mapping is not associated.",err,error,*999)

    rhsMapping=>interfaceMapping%rhsMapping
    IF(.NOT.ASSOCIATED(rhsMapping)) &
      & CALL FlagError("Interface mapping RHS mapping is not associated.",err,error,*999)
       
    EXITS("InterfaceMapping_RHSMappingGet")
    RETURN
999 NULLIFY(rhsMapping)
998 ERRORSEXITS("InterfaceMapping_RHSMappingGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_RHSMappingGet

  !
  !================================================================================================================================
  !
  
END MODULE InterfaceMappingAccessRoutines
