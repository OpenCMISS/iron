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

  PUBLIC InterfaceMapping_InterfaceMatrixToVarMapGet

  PUBLIC InterfaceMapping_LagrangeVariableGet

  PUBLIC InterfaceMapping_MatrixEquationsSetGet

  PUBLIC InterfaceMapping_MatrixVariableGet

  PUBLIC InterfaceMapping_NumberOfColumnsGet

  PUBLIC InterfaceMapping_NumberOfGlobalColumnsGet

  PUBLIC InterfaceMapping_NumberOfInterfaceMatricesGet

  PUBLIC InterfaceMapping_RHSMappingExists

  PUBLIC InterfaceMapping_RHSMappingGet

  PUBLIC InterfaceMapping_TotalNumberOfColumnsGet

  PUBLIC InterfaceMappingCVC_HasTransposeGet

  PUBLIC InterfaceMappingCVC_MatrixCoefficientGet

  PUBLIC InterfaceMappingCVC_MatrixColVariableIndexGet

  PUBLIC InterfaceMappingCVC_MatrixRowVariableIndexGet

  PUBLIC InterfaceMappingCVC_TransposeMatrixCoefficientGet

  PUBLIC InterfaceMappingIMToVMap_HasTransposeGet

  PUBLIC InterfaceMappingIMToVMap_InterfaceMatrixGet

  PUBLIC InterfaceMappingIMToVMap_MatrixCoefficientGet

  PUBLIC InterfaceMappingIMToVMap_MeshIndexGet

  PUBLIC InterfaceMappingIMToVMap_NumberOfRowsGet

  PUBLIC InterfaceMappingIMToVMap_NumberOfGlobalRowsGet

  PUBLIC InterfaceMappingIMToVMap_RowDOFsMappingGet

  PUBLIC InterfaceMappingIMToVMap_TotalNumberOfRowsGet

  PUBLIC InterfaceMappingIMToVMap_TransposeMatrixCoefficientGet

  PUBLIC InterfaceMappingIMToVMap_VariableGet

  PUBLIC InterfaceMappingIMToVMap_VariableDOFToRowMapGet

  PUBLIC InterfaceMappingRHS_InterfaceMappingGet

  PUBLIC InterfaceMappingRHS_InterfaceRowToRHSDOFMapGet

  PUBLIC InterfaceMappingRHS_RHSDOFToInterfaceRowMapGet

  PUBLIC InterfaceMappingRHS_RHSVariableGet
  
  PUBLIC InterfaceMappingRHS_RHSVariableTypeGet

  PUBLIC InterfaceMappingRHS_VectorCoefficientGet
  
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
  SUBROUTINE InterfaceMapping_InterfaceMatrixToVarMapGet(interfaceMapping,matrixIdx,interfaceMatrixToVarMap,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping !<A pointer to the interface mapping to get the interface matrix rows to variable map for
    INTEGER(INTG), INTENT(IN) :: matrixIdx !<The interface matrix index to get the interface matrix rows to variable map for
    TYPE(InterfaceMatrixToVarMapType), POINTER :: interfaceMatrixToVarMap !<On exit, a pointer to the specified interface matrix to variable map in the interface mapping. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif
 
    ENTERS("InterfaceMapping_InterfaceMatrixToVarMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(interfaceMatrixToVarMap)) &
      & CALL FlagError("Interface matrix to variable map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMapping)) CALL FlagError("Interface mapping is not associated.",err,error,*999)
    IF(matrixIdx<1.OR.matrixIdx>interfaceMapping%numberOfInterfaceMatrices) THEN
      localError="The specified matrix index of "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
        & " is invalid. The matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(interfaceMapping%numberOfInterfaceMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(interfaceMapping%interfaceMatrixToVarMaps)) &
      & CALL FlagError("The interface matrix to variable maps is not allocated for the interface mapping.",err,error,*999)
#endif    

    interfaceMatrixToVarMap=>interfaceMapping%interfaceMatrixToVarMaps(matrixIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(interfaceMatrixToVarMap)) THEN
      localError="The interface matrix to variable map is not associated for matrix index "// &
        & TRIM(NumberToVString(matrixIdx,"*",err,error))//" in the interface mapping."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
    EXITS("InterfaceMapping_InterfaceMatrixToVarMapGet")
    RETURN
999 NULLIFY(interfaceMatrixToVarMap)
998 ERRORSEXITS("InterfaceMapping_InterfaceMatrixToVarMapGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_InterfaceMatrixToVarMapGet

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
    IF(.NOT.ALLOCATED(interfaceMapping%interfaceMatrixToVarMaps)) &
      & CALL FlagError("The interface matrix rows to variable maps is not allocated for the interface mapping.",err,error,*999)
    IF(.NOT.ASSOCIATED(interfaceMapping%interfaceMatrixToVarMaps(matrixNumber)%ptr)) THEN
      localError="The interface matrix rows to variable maps is not associated for matrix number "// &
        & TRIM(NumberToVString(matrixNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    equationsSet=>interfaceMapping%interfaceMatrixToVarMaps(matrixNumber)%ptr%equationsSet

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
    IF(.NOT.ALLOCATED(interfaceMapping%interfaceMatrixToVarMaps)) &
      & CALL FlagError("The interface matrix rows to variable maps is not allocated for the interface mapping.",err,error,*999)
    IF(.NOT.ASSOCIATED(interfaceMapping%interfaceMatrixToVarMaps(matrixNumber)%ptr)) THEN
      localError="The interface matrix rows to variable maps is not associated for matrix number "// &
        & TRIM(NumberToVString(matrixNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    matrixVariable=>interfaceMapping%interfaceMatrixToVarMaps(matrixNumber)%ptr%variable

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

  !>Gets the number of columns of matrices in an interface mapping.
  SUBROUTINE InterfaceMapping_NumberOfColumnsGet(interfaceMapping,numberOfColumns,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping !<A pointer to the interface mapping to get the number of columns for
    INTEGER(INTG), INTENT(OUT) :: numberOfColumns !<On exit, the number of columns of matrices in the interface mapping.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMapping_NumberOfColumnsGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(interfaceMapping)) CALL FlagError("Interface mapping is not associated.",err,error,*999)
#endif    

    numberOfColumns=interfaceMapping%numberOfColumns
       
    EXITS("InterfaceMapping_NumberOfColumnsGet")
    RETURN
999 ERRORSEXITS("InterfaceMapping_NumberOfColumnsGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_NumberOfColumnsGet

  !
  !================================================================================================================================
  !

  !>Gets the number of global columns of matrices in an interface mapping.
  SUBROUTINE InterfaceMapping_NumberOfGlobalColumnsGet(interfaceMapping,numberOfGlobalColumns,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping !<A pointer to the interface mapping to get the number of global columns for
    INTEGER(INTG), INTENT(OUT) :: numberOfGlobalColumns !<On exit, the number of global columns of matrices in the interface mapping.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMapping_NumberOfGlobalColumnsGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(interfaceMapping)) CALL FlagError("Interface mapping is not associated.",err,error,*999)
#endif    

    numberOfGlobalColumns=interfaceMapping%numberOfGlobalColumns
       
    EXITS("InterfaceMapping_NumberOfGlobalColumnsGet")
    RETURN
999 ERRORSEXITS("InterfaceMapping_NumberOfGlobalColumnsGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_NumberOfGlobalColumnsGet

  !
  !================================================================================================================================
  !

  !>Gets the number of interface matrices in an interface mapping.
  SUBROUTINE InterfaceMapping_NumberOfInterfaceMatricesGet(interfaceMapping,numberOfInterfaceMatrices,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping !<A pointer to the interface mapping to get the number of interface matrices for
    INTEGER(INTG), INTENT(OUT) :: numberOfInterfaceMatrices !<On exit, the number of the interface matrices in the interface mapping.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMapping_NumberOfInterfaceMatricesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(interfaceMapping)) CALL FlagError("Interface mapping is not associated.",err,error,*999)
#endif    

    numberOfInterfaceMatrices=interfaceMapping%numberOfInterfaceMatrices
       
    EXITS("InterfaceMapping_NumberOfInterfaceMatricesGet")
    RETURN
999 ERRORSEXITS("InterfaceMapping_NumberOfInterfaceMatricesGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_NumberOfInterfaceMatricesGet

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

  !>Gets the total number of columns of matrices in an interface mapping.
  SUBROUTINE InterfaceMapping_TotalNumberOfColumnsGet(interfaceMapping,totalNumberOfColumns,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping !<A pointer to the interface mapping to get the total number of columns for
    INTEGER(INTG), INTENT(OUT) :: totalNumberOfColumns !<On exit, the total number of columns of matrices in the interface mapping.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMapping_TotalNumberOfColumnsGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(interfaceMapping)) CALL FlagError("Interface mapping is not associated.",err,error,*999)
#endif    

    totalNumberOfColumns=interfaceMapping%totalNumberOfColumns
       
    EXITS("InterfaceMapping_TotalNumberOfColumnsGet")
    RETURN
999 ERRORSEXITS("InterfaceMapping_TotalNumberOfColumnsGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMapping_TotalNumberOfColumnsGet

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
    IF(.NOT.ALLOCATED(createValuesCache%hasTranspose)) &
      & CALL FlagError("The has transpose array is not allocated for the create values cache.",err,error,*999)
#endif    

    hasTranspose=createValuesCache%hasTranspose(matrixIdx)
       
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
    IF(.NOT.ALLOCATED(createValuesCache%matrixCoefficients)) &
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
    IF(.NOT.ALLOCATED(createValuesCache%matrixColFieldVariableIndices)) &
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
  SUBROUTINE InterfaceMappingCVC_MatrixRowVariableIndexGet(createValuesCache,matrixIdx,rowVariableIndex,err,error,*)

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
 
    ENTERS("InterfaceMappingCVC_MatrixRowVariableIndexGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(createValuesCache)) CALL FlagError("Create values cache is not associated.",err,error,*999)
    IF(matrixIdx<1.OR.matrixIdx>createValuesCache%numberOfInterfaceMatrices) THEN
      localError="The specified matrix index of "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
        & " is invalid. The matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(createValuesCache%numberOfInterfaceMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(createValuesCache%matrixRowFieldVariableIndices)) &
      & CALL FlagError("The matrix row variable indices array is not allocated for the create values cache.",err,error,*999)
#endif    

    rowVariableIndex=createValuesCache%matrixRowFieldVariableIndices(matrixIdx)
       
    EXITS("InterfaceMappingCVC_MatrixRowVariableIndexGet")
    RETURN
999 ERRORSEXITS("InterfaceMappingCVC_MatrixRowVariableIndexGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMappingCVC_MatrixRowVariableIndexGet

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
    IF(.NOT.ALLOCATED(createValuesCache%transposeMatrixCoefficients)) &
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

  !>Gets the has transpose flag for an interface matrix to variable map.
  SUBROUTINE InterfaceMappingIMToVMap_HasTransposeGet(interfaceMatrixToVarMap,hasTranspose,err,error,*)

    !Argument variables
    TYPE(InterfaceMatrixToVarMapType), POINTER :: interfaceMatrixToVarMap !<A pointer to the interface matrix to variable map to get the has transpose flag for
    LOGICAL, INTENT(OUT) :: hasTranspose !<On exit, the has transpose flag for the interface matrix to variable map.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMappingIMToVMap_HasTransposeGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(interfaceMatrixToVarMap)) &
      & CALL FlagError("Interface matrix to variable map is not associated.",err,error,*999)
#endif    

    hasTranspose=interfaceMatrixToVarMap%hasTranspose
       
    EXITS("InterfaceMappingIMToVMap_HasTransposeGet")
    RETURN
999 ERRORSEXITS("InterfaceMappingIMToVMap_HasTransposeGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMappingIMToVMap_HasTransposeGet

  !
  !================================================================================================================================
  !

  !>Gets the interface matrix for an interface matrix to variable map.
  SUBROUTINE InterfaceMappingIMToVMap_InterfaceMatrixGet(interfaceMatrixToVarMap,interfaceMatrix,err,error,*)

    !Argument variables
    TYPE(InterfaceMatrixToVarMapType), POINTER :: interfaceMatrixToVarMap !<A pointer to the interface matrix to variable map to get the interface matrix for
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix !<On exit, the interface matrix for the interface matrix to variable map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMappingIMToVMap_InterfaceMatrixGet",err,error,*998)

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
       
    EXITS("InterfaceMappingIMToVMap_InterfaceMatrixGet")
    RETURN
999 NULLIFY(interfaceMatrix)
998 ERRORSEXITS("InterfaceMappingIMToVMap_InterfaceMatrixGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMappingIMToVMap_InterfaceMatrixGet

  !
  !================================================================================================================================
  !

  !>Gets the matrix coefficient for an interface matrix to variable map.
  SUBROUTINE InterfaceMappingIMToVMap_MatrixCoefficientGet(interfaceMatrixToVarMap,matrixCoefficient,err,error,*)

    !Argument variables
    TYPE(InterfaceMatrixToVarMapType), POINTER :: interfaceMatrixToVarMap !<A pointer to the interface matrix to variable map to get the matrix coefficient for
    REAL(DP), INTENT(OUT) :: matrixCoefficient !<On exit, the matrix coefficient for the interface matrix to variable map.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMappingIMToVMap_MatrixCoefficientGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(interfaceMatrixToVarMap)) &
      & CALL FlagError("Interface matrix to variable map is not associated.",err,error,*999)
#endif    

    matrixCoefficient=interfaceMatrixToVarMap%matrixCoefficient
       
    EXITS("InterfaceMappingIMToVMap_MatrixCoefficientGet")
    RETURN
999 ERRORSEXITS("InterfaceMappingIMToVMap_MatrixCoefficientGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMappingIMToVMap_MatrixCoefficientGet

  !
  !================================================================================================================================
  !

  !>Gets the mesh index for an interface matrix to variable map.
  SUBROUTINE InterfaceMappingIMToVMap_MeshIndexGet(interfaceMatrixToVarMap,meshIndex,err,error,*)

    !Argument variables
    TYPE(InterfaceMatrixToVarMapType), POINTER :: interfaceMatrixToVarMap !<A pointer to the interface matrix to variable map to get the mesh index for
    INTEGER(INTG), INTENT(OUT) :: meshIndex !<On exit, the mesh index for the interface matrix to variable map.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMappingIMToVMap_MeshIndexGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(interfaceMatrixToVarMap)) &
      & CALL FlagError("Interface matrix to variable map is not associated.",err,error,*999)
#endif    

    meshIndex=interfaceMatrixToVarMap%meshIndex
       
    EXITS("InterfaceMappingIMToVMap_MeshIndexGet")
    RETURN
999 ERRORSEXITS("InterfaceMappingIMToVMap_MeshIndexGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMappingIMToVMap_MeshIndexGet

  !
  !================================================================================================================================
  !

  !>Gets the number of rows for an interface matrix to variable map.
  SUBROUTINE InterfaceMappingIMToVMap_NumberOfRowsGet(interfaceMatrixToVarMap,numberOfRows,err,error,*)

    !Argument variables
    TYPE(InterfaceMatrixToVarMapType), POINTER :: interfaceMatrixToVarMap !<A pointer to the interface matrix to variable map to get the number of rows for
    INTEGER(INTG), INTENT(OUT) :: numberOfRows !<On exit, the number of rows for the interface matrix to variable map.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMappingIMToVMap_NumberOfRowsGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(interfaceMatrixToVarMap)) &
      & CALL FlagError("Interface matrix to variable map is not associated.",err,error,*999)
#endif    

    numberOfRows=interfaceMatrixToVarMap%numberOfRows
       
    EXITS("InterfaceMappingIMToVMap_NumberOfRowsGet")
    RETURN
999 ERRORSEXITS("InterfaceMappingIMToVMap_NumberOfRowsGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMappingIMToVMap_NumberOfRowsGet

  !
  !================================================================================================================================
  !

  !>Gets the number of global rows for an interface matrix to variable map.
  SUBROUTINE InterfaceMappingIMToVMap_NumberOfGlobalRowsGet(interfaceMatrixToVarMap,numberOfGlobalRows,err,error,*)

    !Argument variables
    TYPE(InterfaceMatrixToVarMapType), POINTER :: interfaceMatrixToVarMap !<A pointer to the interface matrix to variable map to get the number of global rows for
    INTEGER(INTG), INTENT(OUT) :: numberOfGlobalRows !<On exit, the number of global rows for the interface matrix to variable map.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMappingIMToVMap_NumberOfGlobalRowsGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(interfaceMatrixToVarMap)) &
      & CALL FlagError("Interface matrix to variable map is not associated.",err,error,*999)
#endif    

    numberOfGlobalRows=interfaceMatrixToVarMap%numberOfGlobalRows
       
    EXITS("InterfaceMappingIMToVMap_NumberOfGlobalRowsGet")
    RETURN
999 ERRORSEXITS("InterfaceMappingIMToVMap_NumberOfGlobalRowsGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMappingIMToVMap_NumberOfGlobalRowsGet

  !
  !================================================================================================================================
  !

  !>Gets the row DOFs mapping for an interface matrix to variable map.
  SUBROUTINE InterfaceMappingIMToVMap_RowDOFsMappingGet(interfaceMatrixToVarMap,rowDOFsMapping,err,error,*)

    !Argument variables
    TYPE(InterfaceMatrixToVarMapType), POINTER :: interfaceMatrixToVarMap !<A pointer to the interface matrix to variable map to get the row DOFs mapping for
    TYPE(DomainMappingType), POINTER :: rowDOFsMapping !<On exit, the row DOFs mapping for the interface matrix to variable map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMappingIMToVMap_RowDOFsMappingGet",err,error,*998)

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
       
    EXITS("InterfaceMappingIMToVMap_RowDOFsMappingGet")
    RETURN
999 NULLIFY(rowDOFsMapping)
998 ERRORSEXITS("InterfaceMappingIMToVMap_RowDOFsMappingGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMappingIMToVMap_RowDOFsMappingGet

  !
  !================================================================================================================================
  !

  !>Gets the total number of rows for an interface matrix to variable map.
  SUBROUTINE InterfaceMappingIMToVMap_TotalNumberOfRowsGet(interfaceMatrixToVarMap,totalNumberOfRows,err,error,*)

    !Argument variables
    TYPE(InterfaceMatrixToVarMapType), POINTER :: interfaceMatrixToVarMap !<A pointer to the interface matrix to variable map to get the total number of rows for
    INTEGER(INTG), INTENT(OUT) :: totalNumberOfRows !<On exit, the total number of rows for the interface matrix to variable map.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMappingIMToVMap_TotalNumberOfRowsGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(interfaceMatrixToVarMap)) &
      & CALL FlagError("Interface matrix to variable map is not associated.",err,error,*999)
#endif    

    totalNumberOfRows=interfaceMatrixToVarMap%totalNumberOfRows
       
    EXITS("InterfaceMappingIMToVMap_TotalNumberOfRowsGet")
    RETURN
999 ERRORSEXITS("InterfaceMappingIMToVMap_TotalNumberOfRowsGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMappingIMToVMap_TotalNumberOfRowsGet

  !
  !================================================================================================================================
  !

  !>Gets the transpose matrix coefficient for an interface matrix to variable map.
  SUBROUTINE InterfaceMappingIMToVMap_TransposeMatrixCoefficientGet(interfaceMatrixToVarMap,transposeMatrixCoefficient,err,error,*)

    !Argument variables
    TYPE(InterfaceMatrixToVarMapType), POINTER :: interfaceMatrixToVarMap !<A pointer to the interface matrix to variable map to get the transpose matrix coefficient for
    REAL(DP), INTENT(OUT) :: transposeMatrixCoefficient !<On exit, the transpose matrix coefficient for the interface matrix to variable map.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMappingIMToVMap_TransposeMatrixCoefficientGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(interfaceMatrixToVarMap)) &
      & CALL FlagError("Interface matrix to variable map is not associated.",err,error,*999)
#endif    

    IF(interfaceMatrixToVarMap%hasTranspose) THEN
      transposeMatrixCoefficient=interfaceMatrixToVarMap%transposeMatrixCoefficient
    ELSE
      transposeMatrixCoefficient=0.0_DP
    ENDIF
       
    EXITS("InterfaceMappingIMToVMap_TransposeMatrixCoefficientGet")
    RETURN
999 ERRORSEXITS("InterfaceMappingIMToVMap_TransposeMatrixCoefficientGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMappingIMToVMap_TransposeMatrixCoefficientGet

  !
  !================================================================================================================================
  !

  !>Gets the dependent variable for an interface matrix to variable map.
  SUBROUTINE InterfaceMappingIMToVMap_VariableGet(interfaceMatrixToVarMap,fieldVariable,err,error,*)

    !Argument variables
    TYPE(InterfaceMatrixToVarMapType), POINTER :: interfaceMatrixToVarMap !<A pointer to the interface matrix to variable map to get the dependent variable for
    TYPE(FieldVariableType), POINTER :: fieldVariable !<On exit, the dependent field variable for the interface matrix to variable map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMappingIMToVMap_VariableGet",err,error,*998)

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
       
    EXITS("InterfaceMappingIMToVMap_VariableGet")
    RETURN
999 NULLIFY(fieldVariable)
998 ERRORSEXITS("InterfaceMappingIMToVMap_VariableGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMappingIMToVMap_VariableGet

  !
  !================================================================================================================================
  !

  !>Gets the variable DOF to row map for an interface matrix to variable map.
  SUBROUTINE InterfaceMappingIMToVMap_VariableDOFToRowMapGet(interfaceMatrixToVarMap,variableDOFToRowMap,err,error,*)

    !Argument variables
    TYPE(InterfaceMatrixToVarMapType), POINTER :: interfaceMatrixToVarMap !<A pointer to the interface matrix to variable map to get the variable DOF to row map for
    INTEGER(INTG), POINTER :: variableDOFToRowMap(:) !<On exit, the variable DOF to row map for the interface matrix to variable map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMappingIMToVMap_VariableDOFToRowMapGet",err,error,*998)

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
       
    EXITS("InterfaceMappingIMToVMap_VariableDOFToRowMapGet")
    RETURN
999 NULLIFY(variableDOFToRowMap)
998 ERRORSEXITS("InterfaceMappingIMToVMap_VariableDOFToRowMapGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceMappingIMToVMap_VariableDOFToRowMapGet

  !
  !================================================================================================================================
  !

  !>Gets the interface mapping for a RHS mapping.
  SUBROUTINE InterfaceMappingRHS_InterfaceMappingGet(rhsMapping,interfaceMapping,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingRHSType), POINTER :: rhsMapping !<A pointer to the rhs mapping to get the interface mapping for
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping !<On exit, a pointer to the interface mapping for the RHS mapping. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMappingRHS_InterfaceMappingGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(interfaceMapping)) CALL FlagError("Interface mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(rhsMapping)) CALL FlagError("Interface RHS mapping is not associated.",err,error,*999)
#endif    
    
    interfaceMapping=>rhsMapping%interfaceMapping

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(interfaceMapping)) &
      & CALL FlagError("The interface mapping is not associated for the interface RHS mapping.",err,error,*999)
#endif    
    
    EXITS("InterfaceMappingRHS_InterfaceMappingGet")
    RETURN
999 NULLIFY(interfaceMapping)
998 ERRORS("InterfaceMappingRHS_InterfaceMappingGet",err,error)
    EXITS("InterfaceMappingRHS_InterfaceMappingGet")
    RETURN 1
    
  END SUBROUTINE InterfaceMappingRHS_InterfaceMappingGet

  !
  !================================================================================================================================
  !

  !>Gets the interface row to RHS DOF map for a RHS mapping.
  SUBROUTINE InterfaceMappingRHS_InterfaceRowToRHSDOFMapGet(rhsMapping,interfaceRowToRHSDOFMap,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingRHSType), POINTER :: rhsMapping !<A pointer to the rhs mapping to get the interface row to RHS DOF map for
    INTEGER(INTG), POINTER :: interfaceRowToRHSDOFMap(:) !<On exit, a pointer to the interface row to RHS DOF map for the RHS mapping. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMappingRHS_InterfaceRowToRHSDOFMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(interfaceRowToRHSDOFMap)) CALL FlagError("Interface row to RHS DOF map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(rhsMapping)) CALL FlagError("RHS mapping is not associated.",err,error,*999)
#endif    
    
    interfaceRowToRHSDOFMap=>rhsMapping%interfaceRowToRHSDOFMap

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(interfaceRowToRHSDOFMap)) &
      & CALL FlagError("The interface row to RHS DOF map is not associated for the RHS mapping.",err,error,*999)
#endif    
    
    EXITS("InterfaceMappingRHS_InterfaceRowToRHSDOFMapGet")
    RETURN
999 NULLIFY(interfaceRowToRHSDOFMap)
998 ERRORS("InterfaceMappingRHS_InterfaceRowToRHSDOFMapGet",err,error)
    EXITS("InterfaceMappingRHS_InterfaceRowToRHSDOFMapGet")
    RETURN 1
    
  END SUBROUTINE InterfaceMappingRHS_InterfaceRowToRHSDOFMapGet

  !
  !================================================================================================================================
  !

  !>Gets the RHS DOF to interface row  map for a RHS mapping.
  SUBROUTINE InterfaceMappingRHS_RHSDOFToInterfaceRowMapGet(rhsMapping,rhsDOFToInterfaceRowMap,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingRHSType), POINTER :: rhsMapping !<A pointer to the rhs mapping to get the RHS DOF to interface row map for
    INTEGER(INTG), POINTER :: rhsDOFToInterfaceRowMap(:) !<On exit, a pointer to the RHS DOF to interface row map for the RHS mapping. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMappingRHS_RHSDOFToInterfaceRowMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(rhsDOFToInterfaceRowMap)) CALL FlagError("RHS DOF to interface row map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(rhsMapping)) CALL FlagError("RHS mapping is not associated.",err,error,*999)
#endif    
    
    rhsDOFToInterfaceRowMap=>rhsMapping%rhsDOFToInterfaceRowMap

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(rhsDOFToInterfaceRowMap)) &
      & CALL FlagError("The RHS DOF to interface row map is not associated for the RHS mapping.",err,error,*999)
#endif    
    
    EXITS("InterfaceMappingRHS_RHSDOFToInterfaceRowMapGet")
    RETURN
999 NULLIFY(rhsDOFToInterfaceRowMap)
998 ERRORS("InterfaceMappingRHS_RHSDOFToInterfaceRowMapGet",err,error)
    EXITS("InterfaceMappingRHS_RHSDOFToInterfaceRowMapGet")
    RETURN 1
    
  END SUBROUTINE InterfaceMappingRHS_RHSDOFToInterfaceRowMapGet

  !
  !================================================================================================================================
  !

  !>Gets the specified RHS variable for a rhs mapping.
  SUBROUTINE InterfaceMappingRHS_RHSVariableGet(rhsMapping,rhsVariable,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingRHSType), POINTER :: rhsMapping !<A pointer to the RHS mapping to get the RHS variable for
    TYPE(FieldVariableType), POINTER :: rhsVariable !<On exit, a pointer to the requested RHS field variable. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMappingRHS_RHSVariableGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(rhsVariable)) CALL FlagError("RHS Field variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(rhsMapping)) CALL FlagError("RHS mapping is not associated.",err,error,*999)
#endif    
    
    rhsVariable=>rhsMapping%rhsVariable

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(rhsVariable)) CALL FlagError("The RHS field variable is not associated in the RHS mapping.",err,error,*999)
#endif    
    
    EXITS("InterfaceMappingRHS_RHSVariableGet")
    RETURN
999 NULLIFY(rhsVariable)
998 ERRORS("InterfaceMappingRHS_RHSVariableGet",err,error)
    EXITS("InterfaceMappingRHS_RHSVariableGet")
    RETURN 1
    
  END SUBROUTINE InterfaceMappingRHS_RHSVariableGet

  !
  !================================================================================================================================
  !

  !>Gets the specified RHS variable type for a rhs mapping.
  SUBROUTINE InterfaceMappingRHS_RHSVariableTypeGet(rhsMapping,rhsVariableType,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingRHSType), POINTER :: rhsMapping !<A pointer to the RHS mapping to get the RHS variable for
    INTEGER(INTG), INTENT(OUT) :: rhsVariableType !<On exit, the requested RHS variable type.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceMappingRHS_RHSVariableTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(rhsMapping)) CALL FlagError("RHS mapping is not associated.",err,error,*999)
#endif    
    
    rhsVariableType=rhsMapping%rhsVariableType
    
    EXITS("InterfaceMappingRHS_RHSVariableTypeGet")
    RETURN
999 ERRORS("InterfaceMappingRHS_RHSVariableTypeGet",err,error)
    EXITS("InterfaceMappingRHS_RHSVariableTypeGet")
    RETURN 1
    
  END SUBROUTINE InterfaceMappingRHS_RHSVariableTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the RHS vector coefficient for a RHS mapping.
  SUBROUTINE InterfaceMappingRHS_VectorCoefficientGet(rhsMapping,vectorCoefficient,err,error,*)

    !Argument variables
    TYPE(InterfaceMappingRHSType), POINTER :: rhsMapping !<A pointer to the RHS mapping to get RHS vector coefficient for
    REAL(DP), INTENT(OUT) :: vectorCoefficient !<On exit, the vector coefficient for a RHS vector of the RHS mapping.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("InterfaceMappingRHS_VectorCoefficientGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(rhsMapping)) CALL FlagError("RHS mapping is not associated.",err,error,*999)
#endif    
      
    vectorCoefficient=rhsMapping%rhsCoefficient
   
    EXITS("InterfaceMappingRHS_VectorCoefficientGet")
    RETURN
999 ERRORS("InterfaceMappingRHS_VectorCoefficientGet",err,error)
    EXITS("InterfaceMappingRHS_VectorCoefficientGet")
    RETURN 1
    
  END SUBROUTINE InterfaceMappingRHS_VectorCoefficientGet

  !
  !================================================================================================================================
  !
  
END MODULE InterfaceMappingAccessRoutines
