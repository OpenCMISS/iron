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

  PUBLIC InterfaceMapping_ColumnDOFsMappingGet

  PUBLIC InterfaceMapping_CreateValuesCacheGet

  PUBLIC InterfaceMapping_InterfaceEquationsGet

  PUBLIC InterfaceMapping_InterfaceMatrixRowsToVarMapGet

  PUBLIC InterfaceMapping_LagrangeVariableGet

  PUBLIC InterfaceMapping_MatrixEquationsSetGet

  PUBLIC InterfaceMapping_MatrixVariableGet

  PUBLIC InterfaceMapping_RHSMappingExists

  PUBLIC InterfaceMapping_RHSMappingGet

  PUBLIC InterfaceMappingCVC_HasTransposeGet

  PUBLIC InterfaceMappingCVC_MatrixCoefficientGet

  PUBLIC InterfaceMappingCVS_MatrixColVariableIndexGet

  PUBLIC InterfaceMappingCVS_MatrixRowVariableIndexGet

  PUBLIC InterfaceMappingCVC_TransposeMatrixCoefficientGet

  PUBLIC InterfaceMappingIMToVM_InterfaceMatrixGet

  PUBLIC InterfaceMappingIMToVM_MeshIndexGet

  PUBLIC InterfaceMappingIMToVM_RowDOFsMappingGet

  PUBLIC InterfaceMappingIMToVM_VariableGet

  PUBLIC InterfaceMappingIMToVM_VariableDOFToRowMapGet

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

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(interfaceMapping)) CALL FlagError("Interface mapping is not associated.",err,error,*999)
#endif    

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

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(interfaceMapping)) CALL FlagError("Interface mapping is not associated.",err,error,*999)
#endif    

    IF(interfaceMapping%interfaceMappingFinished) CALL FlagError("Interface mapping has already been finished.",err,error,*999)
    
    EXITS("InterfaceMapping_AssertNotFinished")
    RETURN
999 ERRORSEXITS("InterfaceMapping_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_AssertNotFinished

  !
  !================================================================================================================================
  !

  !>Gets the column DOFs mapping for an interface mapping.
  SUBROUTINE InterfaceMapping_ColumnDOFsMappingGet(interfaceMapping,columnDOFsMapping,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping !<A pointer to the interface mapping to get the column DOFs mapping for
    TYPE(DomainMappingType), POINTER :: columnDOFsMapping !<On exit, a pointer to the column DOFs mapping in the interface mapping. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMapping_ColumnDOFsMappingGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(columnDOFsMapping)) CALL FlagError("Column DOFs mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMapping)) CALL FlagError("Interface mapping is not associated.",err,error,*999)
#endif    

   columnDOFsMapping=>interfaceMapping%columnDOFSMapping

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(columnDOFsMapping)) &
      & CALL FlagError("Interface mapping column DOFs mapping is not associated.",err,error,*999)
#endif    
       
    EXITS("InterfaceMapping_ColumnDOFsMappingGet")
    RETURN
999 NULLIFY(columnDOFsMapping)
998 ERRORSEXITS("InterfaceMapping_ColumnDOFsMappingGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_ColumnDOFsMappingGet

  !
  !================================================================================================================================
  !

  !>Gets the create values cache for an interface mapping.
  SUBROUTINE InterfaceMapping_CreateValuesCacheGet(interfaceMapping,createValuesCache,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping !<A pointer to the interface mapping to get the create values cache for
    TYPE(InterfaceMappingCreateValuesCacheType), POINTER :: createValuesCache !<On exit, a pointer to the create values cache in the specified interface mapping. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMapping_CreateValuesCacheGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(createValuesCache)) CALL FlagError("Create values cache is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMapping)) CALL FlagError("Interface mapping is not associated.",err,error,*999)
#endif    

    createValuesCache=>interfaceMapping%createValuesCache

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(createValuesCache)) &
      & CALL FlagError("Interface mapping create values cache is not associated.",err,error,*999)
#endif    
       
    EXITS("InterfaceMapping_CreateValuesCacheGet")
    RETURN
999 NULLIFY(createValuesCache)
998 ERRORSEXITS("InterfaceMapping_CreateValuesCacheGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_CreateValuesCacheGet

  !
  !================================================================================================================================
  !

  !>Gets the interface equations for an interface mapping.
  SUBROUTINE InterfaceMapping_InterfaceEquationsGet(interfaceMapping,interfaceEquations,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping !<A pointer to the interface mapping to get the interface equations for
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations !<On exit, a pointer to the interface equations in the specified interface mapping. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMapping_InterfaceEquationsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(interfaceEquations)) CALL FlagError("Interface equations is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMapping)) CALL FlagError("Interface mapping is not associated.",err,error,*999)
#endif    

    interfaceEquations=>interfaceMapping%interfaceEquations

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(interfaceEquations)) &
      & CALL FlagError("Interface mapping interface equations is not associated.",err,error,*999)
#endif    
       
    EXITS("InterfaceMapping_InterfaceEquationsGet")
    RETURN
999 NULLIFY(interfaceEquations)
998 ERRORSEXITS("InterfaceMapping_InterfaceEquationsGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_InterfaceEquationsGet

  !
  !================================================================================================================================
  !

  !>Gets a interface matrix rows to variable map for an interface mapping.
  SUBROUTINE InterfaceMapping_InterfaceMatrixRowsToVarMapGet(interfaceMapping,matrixIdx,interfaceMatrixRowsToVarMap,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping !<A pointer to the interface mapping to get the interface matrix rows to variable map for
    INTEGER(INTG), INTENT(IN) :: matrixIdx !<The interface matrix index to get the interface matrix rows to variable map for
    TYPE(InterfaceMatrixToVarMapType), POINTER :: interfaceMatrixRowsToVarMap !<On exit, a pointer to the specified interface matrix rows to variable map in the interface mapping. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif
 
    ENTERS("InterfaceMapping_InterfaceMatrixRowsToVarMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(interfaceMatrixRowsToVarMap)) &
      & CALL FlagError("Interface matrix rows to variable map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMapping)) CALL FlagError("Interface mapping is not associated.",err,error,*999)
    IF(matrixIdx<1.OR.matrixIdx>interfaceMapping%numberOfInterfaceMatrices) THEN
      localError="The specified matrix index of "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
        & " is invalid. The matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(interfaceMapping%numberOfInterfaceMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(interfaceMapping%interfaceMatrixRowsToVarMaps)) &
      & CALL FlagError("The interface matrix rows to variable maps is not allocated for the interface mapping.",err,error,*999)
#endif    

    interfaceMatrixRowsToVarMap=>interfaceMapping%interfaceMatrixRowsToVarMaps(matrixIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(interfaceMatrixRowsToVarMap)) THEN
      localError="The interface matrix rows to variable map is not associated for matrix index "// &
        & TRIM(NumberToVString(matrixIdx,"*",err,error))//" in the interface mapping."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
    EXITS("InterfaceMapping_InterfaceMatrixRowsToVarMapGet")
    RETURN
999 NULLIFY(interfaceMatrixRowsToVarMap)
998 ERRORSEXITS("InterfaceMapping_InterfaceMatrixRowsToVarMapGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_InterfaceMatrixRowsToVarMapGet

  !
  !================================================================================================================================
  !

  !>Gets the lagrange DOF to column map for an interface mapping.
  SUBROUTINE InterfaceMapping_LagrangeDOFToColumnMapGet(interfaceMapping,lagrangeDOFToColumnMap,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping !<A pointer to the interface mapping to get the lagrange DOF to column map for
    INTEGER(INTG), POINTER :: lagrangeDOFToColumnMap(:) !<On exit, a pointer to the lagrange DOF to column map in the specified interface mapping. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMapping_LagrangeDOFToColumnMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(lagrangeDOFToColumnMap)) CALL FlagError("Lagrange DOF to column map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMapping)) CALL FlagError("Interface mapping is not associated.",err,error,*999)
#endif    

    lagrangeDOFToColumnMap=>interfaceMapping%lagrangeDOFToColumnMap

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(lagrangeDOFToColumnMap)) &
      & CALL FlagError("Interface mapping Lagrange DOF to column map is not associated.",err,error,*999)
#endif    
       
    EXITS("InterfaceMapping_LagrangeDOFToColumnMapGet")
    RETURN
999 NULLIFY(lagrangeDOFToColumnMap)
998 ERRORSEXITS("InterfaceMapping_LagrangeDOFToColumnMapGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_LagrangeDOFToColumnMapGet

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

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(lagrangeVariable)) CALL FlagError("Lagrange variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMapping)) CALL FlagError("Interface mapping is not associated.",err,error,*999)
#endif    

    lagrangeVariable=>interfaceMapping%lagrangeVariable

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(lagrangeVariable)) &
      & CALL FlagError("Interface mapping Lagrange variable is not associated.",err,error,*999)
#endif    
       
    EXITS("InterfaceMapping_LagrangeVariableGet")
    RETURN
999 NULLIFY(lagrangeVariable)
998 ERRORSEXITS("InterfaceMapping_LagrangeVariableGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_LagrangeVariableGet

  !
  !================================================================================================================================
  !

  !>Gets the matrix (row) variable equations set for an interface matrix in an interface mapping.
  SUBROUTINE InterfaceMapping_MatrixEquationsSetGet(interfaceMapping,matrixNumber,equationsSet,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping !<A pointer to the interface mapping to get the matrix variable equations set for
    INTEGER(INTG), INTENT(IN) :: matrixNumber !<The number of the matrix to get the variable equations set for
    TYPE(EquationsSetType), POINTER :: equationsSet !<On exit, a pointer to the equations set for the specified interface matrix field variable. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("InterfaceMapping_MatrixEquationsSetGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(equationsSet)) CALL FlagError("Matrix equations set is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMapping)) CALL FlagError("Interface mapping is not associated.",err,error,*999)
    IF(matrixNumber<=0.OR.matrixNumber>interfaceMapping%numberOfInterfaceMatrices) THEN
      localError="The specified matrix number of "//TRIM(NumberToVString(matrixNumber,"*",err,error))// &
        & " is invalid. The matrix number must be >= 1 and <= "// &
        & TRIM(NumberToVString(interfaceMapping%numberOfInterfaceMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(interfaceMapping%interfaceMatrixRowsToVarMaps)) &
      & CALL FlagError("The interface matrix row to variable maps is not allocated for the interface mapping.",err,error,*999)
#endif    

    equationsSet=>interfaceMapping%interfaceMatrixRowsToVarMaps(matrixNumber)%equationsSet

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(equationsSet)) THEN
      localError="The interface mapping matrix variable equations set is not associted for matrix number "// &
        & TRIM(NumberToVString(matrixNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
    EXITS("InterfaceMapping_MatrixEquationsSetGet")
    RETURN
999 NULLIFY(equationsSet)
998 ERRORSEXITS("InterfaceMapping_MatrixEquationsSetGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_MatrixEquationsSetGet

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
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("InterfaceMapping_MatrixVariableGet",err,error,*998)

#ifdef WITH_PRECHECKS    
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
#endif    

    matrixVariable=>interfaceMapping%interfaceMatrixRowsToVarMaps(matrixNumber)%variable

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(matrixVariable)) THEN
      localError="The interface mapping matrix variable is not associted for matrix number "// &
        & TRIM(NumberToVString(matrixNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
    EXITS("InterfaceMapping_MatrixVariableGet")
    RETURN
999 NULLIFY(matrixVariable)
998 ERRORSEXITS("InterfaceMapping_MatrixVariableGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_MatrixVariableGet

  !
  !================================================================================================================================
  !

  !>Checks if the RHS mapping exists for an interface mapping.
  SUBROUTINE InterfaceMapping_RHSMappingExists(interfaceMapping,rhsMapping,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping !<A pointer to the interface mapping to check the RHS mapping for
    TYPE(InterfaceMappingRHSType), POINTER :: rhsMapping !<On exit, a pointer to the RHS mapping in the specified interface mapping if it exits. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMapping_RHSMappingExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(rhsMapping)) CALL FlagError("RHS mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMapping)) CALL FlagError("Interface mapping is not associated.",err,error,*999)
#endif    

    rhsMapping=>interfaceMapping%rhsMapping

       
    EXITS("InterfaceMapping_RHSMappingExists")
    RETURN
999 NULLIFY(rhsMapping)
998 ERRORSEXITS("InterfaceMapping_RHSMappingExists",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_RHSMappingExists

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

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(rhsMapping)) CALL FlagError("RHS mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMapping)) CALL FlagError("Interface mapping is not associated.",err,error,*999)
#endif    

    rhsMapping=>interfaceMapping%rhsMapping

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(rhsMapping)) &
      & CALL FlagError("Interface mapping RHS mapping is not associated.",err,error,*999)
#endif    
       
    EXITS("InterfaceMapping_RHSMappingGet")
    RETURN
999 NULLIFY(rhsMapping)
998 ERRORSEXITS("InterfaceMapping_RHSMappingGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_RHSMappingGet

  !
  !================================================================================================================================
  !

  !>Gets the has transpose for an interface matrix for an interface mapping create values cache.
  SUBROUTINE InterfaceMappingCVC_HasTransposeGet(createValuesCache,matrixIdx,hasTranspose,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingCreateValuesCacheType), POINTER :: createValuesCache !<A pointer to the interface mapping create value cache to get the has transpose for
    INTEGER(INTG), INTENT(IN) :: matrixIdx !<The interface matrix index to get the has transpose for
    LOGICAL, INTENT(OUT) :: hasTranspose !<On exit, the has transpose for the specified interface matrix.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif
 
    ENTERS("InterfaceMappingCVC_HasTransposeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(createValuesCache)) CALL FlagError("Create values cache is not associated.",err,error,*999)
    IF(matrixIdx<1.OR.matrixIdx>createValuesCache%numberOfInterfaceMatrices) THEN
      localError="The specified matrix index of "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
        & " is invalid. The matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(createValuesCache%numberOfInterfaceMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(createValuesCachce%hasTranspose)) &
      & CALL FlagError("The has transpose array is not allocated for the create values cache.",err,error,*999)
#endif    

    hasTranspose=createValuesCache%hasTranpose(matrixIdx)
       
    EXITS("InterfaceMappingCVC_HasTransposeGet")
    RETURN
999 ERRORSEXITS("InterfaceMappingCVC_HasTransposeGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMappingCVC_HasTransposeGet

  !
  !================================================================================================================================
  !

  !>Gets the matrix coefficient for an interface matrix for an interface mapping create values cache.
  SUBROUTINE InterfaceMappingCVC_MatrixCoefficientGet(createValuesCache,matrixIdx,matrixCoefficient,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingCreateValuesCacheType), POINTER :: createValuesCache !<A pointer to the interface mapping create value cache to get the matrix coefficient for
    INTEGER(INTG), INTENT(IN) :: matrixIdx !<The interface matrix index to get the matrix coefficient for
    REAL(DP), INTENT(OUT) :: matrixCoefficient !<On exit, the matrix coefficient for the specified interface matrix.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif
 
    ENTERS("InterfaceMappingCVC_MatrixCoefficientGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(createValuesCache)) CALL FlagError("Create values cache is not associated.",err,error,*999)
    IF(matrixIdx<1.OR.matrixIdx>createValuesCache%numberOfInterfaceMatrices) THEN
      localError="The specified matrix index of "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
        & " is invalid. The matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(createValuesCache%numberOfInterfaceMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(createValuesCachce%matrixCoefficients)) &
      & CALL FlagError("The matrix coefficients array is not allocated for the create values cache.",err,error,*999)
#endif    

    matrixCoefficient=createValuesCache%matrixCoefficients(matrixIdx)
       
    EXITS("InterfaceMappingCVC_MatrixCoefficientGet")
    RETURN
999 ERRORSEXITS("InterfaceMappingCVC_MatrixCoefficientGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMappingCVC_MatrixCoefficientGet

  !
  !================================================================================================================================
  !

  !>Gets the column variable index for an interface matrix for an interface mapping create values cache.
  SUBROUTINE InterfaceMappingCVC_MatrixColVariableIndexGet(createValuesCache,matrixIdx,columnVariableIndex,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingCreateValuesCacheType), POINTER :: createValuesCache !<A pointer to the interface mapping create value cache to get the column variable index for
    INTEGER(INTG), INTENT(IN) :: matrixIdx !<The interface matrix index to get the column variable index for
    INTEGER(INTG), INTENT(OUT) :: columnVariableIndex !<On exit, the column variable index for the specified interface matrix.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif
 
    ENTERS("InterfaceMappingCVC_MatrixColVariableIndexGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(createValuesCache)) CALL FlagError("Create values cache is not associated.",err,error,*999)
    IF(matrixIdx<1.OR.matrixIdx>createValuesCache%numberOfInterfaceMatrices) THEN
      localError="The specified matrix index of "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
        & " is invalid. The matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(createValuesCache%numberOfInterfaceMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(createValuesCachce%matrixColFieldVariableIndices)) &
      & CALL FlagError("The matrix column variable indices array is not allocated for the create values cache.",err,error,*999)
#endif    

    columnVariableIndex=createValuesCache%matrixColFieldVariableIndices(matrixIdx)
       
    EXITS("InterfaceMappingCVC_MatrixColVariableIndexGet")
    RETURN
999 ERRORSEXITS("InterfaceMappingCVC_MatrixColVariableIndexGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMappingCVC_MatrixColVariableIndexGet

  !
  !================================================================================================================================
  !

  !>Gets the row variable index for an interface matrix for an interface mapping create values cache.
  SUBROUTINE InterfaceMappingCVC_MatrixColVariableIndexGet(createValuesCache,matrixIdx,rowVariableIndex,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingCreateValuesCacheType), POINTER :: createValuesCache !<A pointer to the interface mapping create value cache to get the row variable index for
    INTEGER(INTG), INTENT(IN) :: matrixIdx !<The interface matrix index to get the row variable index for
    INTEGER(INTG), INTENT(OUT) :: rowVariableIndex !<On exit, the row variable index for the specified interface matrix.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif
 
    ENTERS("InterfaceMappingCVC_MatrixColVariableIndexGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(createValuesCache)) CALL FlagError("Create values cache is not associated.",err,error,*999)
    IF(matrixIdx<1.OR.matrixIdx>createValuesCache%numberOfInterfaceMatrices) THEN
      localError="The specified matrix index of "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
        & " is invalid. The matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(createValuesCache%numberOfInterfaceMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(createValuesCachce%matrixRowFieldVariableIndices)) &
      & CALL FlagError("The matrix row variable indices array is not allocated for the create values cache.",err,error,*999)
#endif    

    rowVariableIndex=createValuesCache%matrixRowFieldVariableIndices(matrixIdx)
       
    EXITS("InterfaceMappingCVC_MatrixColVariableIndexGet")
    RETURN
999 ERRORSEXITS("InterfaceMappingCVC_MatrixColVariableIndexGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMappingCVC_MatrixColVariableIndexGet

  !
  !================================================================================================================================
  !

  !>Gets the transpose matrix coefficient for an interface matrix for an interface mapping create values cache.
  SUBROUTINE InterfaceMappingCVC_TransposeMatrixCoefficientGet(createValuesCache,matrixIdx,transposeMatrixCoefficient,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingCreateValuesCacheType), POINTER :: createValuesCache !<A pointer to the interface mapping create value cache to get the transpose matrix coefficient for
    INTEGER(INTG), INTENT(IN) :: matrixIdx !<The interface matrix index to get the transpose matrix coefficient for
    REAL(DP), INTENT(OUT) :: transposeMatrixCoefficient !<On exit, the transpose matrix coefficient for the specified interface matrix.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif
 
    ENTERS("InterfaceMappingCVC_TransposeMatrixCoefficientGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(createValuesCache)) CALL FlagError("Create values cache is not associated.",err,error,*999)
    IF(matrixIdx<1.OR.matrixIdx>createValuesCache%numberOfInterfaceMatrices) THEN
      localError="The specified matrix index of "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
        & " is invalid. The matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(createValuesCache%numberOfInterfaceMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(createValuesCachce%transposeMatrixCoefficients)) &
      & CALL FlagError("The transpose matrix coefficients array is not allocated for the create values cache.",err,error,*999)
#endif    

    transposeMatrixCoefficient=createValuesCache%transposeMatrixCoefficients(matrixIdx)
       
    EXITS("InterfaceMappingCVC_TransposeMatrixCoefficientGet")
    RETURN
999 ERRORSEXITS("InterfaceMappingCVC_TransposeMatrixCoefficientGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMappingCVC_TransposeMatrixCoefficientGet

  !
  !================================================================================================================================
  !

  !>Gets the interface matrix for an interface matrix to variable map.
  SUBROUTINE InterfaceMappingIMToVM_InterfaceMatrixGet(interfaceMatrixToVarMap,interfaceMatrix,err,error,*)

    !Argument variables
    TYPE(InterfaceMatrixToVarMapType), POINTER :: interfaceMatrixToVarMap !<A pointer to the interface matrix to variable map to get the interface matrix for
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix !<On exit, the interface matrix for the interface matrix to variable map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMappingIMToVM_InterfaceMatrixGet",err,error,*998)

#ifdef WITH_PRECHECKS
    IF(ASSOCIATED(interfaceMatrix)) CALL FlagError("Interface matrix is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMatrixToVarMap)) &
      & CALL FlagError("Interface matrix to variable map is not associated.",err,error,*999)
#endif    

    interfaceMatrix=>interfaceMatrixToVarMap%interfaceMatrix

#ifdef WITH_POSTCHECKS
    IF(.NOT.ASSOCIATED(interfaceMatrix)) &
      & CALL FlagError("Interface matrix is not associated for the interface matrix to variable map.",err,error,*999)
#endif        
       
    EXITS("InterfaceMappingIMToVM_InterfaceMatrixGet")
    RETURN
999 NULLIFY(interfaceMatrix)
998 ERRORSEXITS("InterfaceMappingIMToVM_InterfaceMatrixGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMappingIMToVM_InterfaceMatrixGet

  !
  !================================================================================================================================
  !

  !>Gets the mesh index for an interface matrix to variable map.
  SUBROUTINE InterfaceMappingIMToVM_MeshIndexGet(interfaceMatrixToVarMap,meshIndex,err,error,*)

    !Argument variables
    TYPE(InterfaceMatrixToVarMapType), POINTER :: interfaceMatrixToVarMap !<A pointer to the interface matrix to variable map to get the mesh index for
    INTEGER(INTG), INTENT(OUT) :: meshIndex !<On exit, the mesh index for the interface matrix to variable map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMappingIMToVM_MeshIndexGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(interfaceMatrixToVarMap)) &
      & CALL FlagError("Interface matrix to variable map is not associated.",err,error,*999)
#endif    

    meshIndex=interfaceMatrixToVarMap%meshIndex
       
    EXITS("InterfaceMappingIMToVM_MeshIndexGet")
    RETURN
999 ERRORSEXITS("InterfaceMappingIMToVM_MeshIndexGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMappingIMToVM_MeshIndexGet

  !
  !================================================================================================================================
  !

  !>Gets the row DOFs mapping for an interface matrix to variable map.
  SUBROUTINE InterfaceMappingIMToVM_RowDOFsMappingGet(interfaceMatrixToVarMap,rowDOFsMapping,err,error,*)

    !Argument variables
    TYPE(InterfaceMatrixToVarMapType), POINTER :: interfaceMatrixToVarMap !<A pointer to the interface matrix to variable map to get the row DOFs mapping for
    TYPE(DomainMappingType), POINTER :: rowDOFsMapping !<On exit, the row DOFs mapping for the interface matrix to variable map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMappingIMToVM_RowDOFsMappingGet",err,error,*998)

#ifdef WITH_PRECHECKS
    IF(ASSOCIATED(rowDOFsMapping)) CALL FlagError("Row DOFs mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMatrixToVarMap)) &
      & CALL FlagError("Interface matrix to variable map is not associated.",err,error,*999)
#endif    

    rowDOFsMapping=>interfaceMatrixToVarMap%rowDOFsMapping

#ifdef WITH_POSTCHECKS
    IF(.NOT.ASSOCIATED(rowDOFsMapping)) &
      & CALL FlagError("Row DOFs mapping is not associated for the interface matrix to variable map.",err,error,*999)
#endif        
       
    EXITS("InterfaceMappingIMToVM_RowDOFsMappingGet")
    RETURN
999 NULLIFY(rowDOFsMapping)
998 ERRORSEXITS("InterfaceMappingIMToVM_RowDOFsMappingGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMappingIMToVM_RowDOFsMappingGet

  !
  !================================================================================================================================
  !

  !>Gets the dependent variable for an interface matrix to variable map.
  SUBROUTINE InterfaceMappingIMToVM_VariableGet(interfaceMatrixToVarMap,fieldVariable,err,error,*)

    !Argument variables
    TYPE(InterfaceMatrixToVarMapType), POINTER :: interfaceMatrixToVarMap !<A pointer to the interface matrix to variable map to get the dependent variable for
    TYPE(FieldVariableType), POINTER :: fieldVariable !<On exit, the dependent field variable for the interface matrix to variable map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMappingIMToVM_VariableGet",err,error,*998)

#ifdef WITH_PRECHECKS
    IF(ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMatrixToVarMap)) &
      & CALL FlagError("Interface matrix to variable map is not associated.",err,error,*999)
#endif    

    fieldVariable=>interfaceMatrixToVarMap%variable

#ifdef WITH_POSTCHECKS
    IF(.NOT.ASSOCIATED(fieldVariable)) &
      & CALL FlagError("Field variable is not associated for the interface matrix to variable map.",err,error,*999)
#endif        
       
    EXITS("InterfaceMappingIMToVM_VariableGet")
    RETURN
999 NULLIFY(fieldVariable)
998 ERRORSEXITS("InterfaceMappingIMToVM_VariableGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMappingIMToVM_VariableGet

  !
  !================================================================================================================================
  !

  !>Gets the variable DOF to row map for an interface matrix to variable map.
  SUBROUTINE InterfaceMappingIMToVM_VariableDOFToRowMapGet(interfaceMatrixToVarMap,variableDOFToRowMap,err,error,*)

    !Argument variables
    TYPE(InterfaceMatrixToVarMapType), POINTER :: interfaceMatrixToVarMap !<A pointer to the interface matrix to variable map to get the variable DOF to row map for
    INTEGER(INTG), POINTER :: variableDOFToRowMap(:) !<On exit, the variable DOF to row map for the interface matrix to variable map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMappingIMToVM_VariableDOFToRowMapGet",err,error,*998)

#ifdef WITH_PRECHECKS
    IF(ASSOCIATED(variableDOFToRowMap)) CALL FlagError("Variable DOF to row map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMatrixToVarMap)) &
      & CALL FlagError("Interface matrix to variable map is not associated.",err,error,*999)
#endif    

    variableDOFToRowMap=>interfaceMatrixToVarMap%variableDOFToRowMap

#ifdef WITH_POSTCHECKS
    IF(.NOT.ASSOCIATED(variableDOFToRowMap)) &
      & CALL FlagError("The variable DOF to row map is not associated for the interface matrix to variable map.",err,error,*999)
#endif        
       
    EXITS("InterfaceMappingIMToVM_VariableDOFToRowMapGet")
    RETURN
999 NULLIFY(variableDOFToRowMap)
998 ERRORSEXITS("InterfaceMappingIMToVM_VariableDOFToRowMapGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMappingIMToVM_VariableDOFToRowMapGet

  !
  !================================================================================================================================
  !
  
END MODULE InterfaceMappingAccessRoutines
