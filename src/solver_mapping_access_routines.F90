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
  USE ISO_VARYING_STRING
  USE Kinds
  USE Strings
  USE Types
  
#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters


  !> \addtogroup SolverMappingRoutines_EquationsMatrixTypes SolverMappingRoutines::EquationsMatrixTypes
  !> \brief Equations matrix types
  !> \see SolverMappingRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_MAPPING_EQUATIONS_DYNAMIC_MATRIX=1 !<The equations matrix in the solver mapping is a dynamic equations matrix \see SolverMappingRoutines_EquationsMatrixTypes,SolverMappingRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_MAPPING_EQUATIONS_LINEAR_MATRIX=2 !<The equations matrix in the solver mapping is a linear equations matrix \see SolverMappingRoutines_EquationsMatrixTypes,SolverMappingRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_MAPPING_EQUATIONS_NONLINEAR_MATRIX=3 !<The equations matrix in the solver mapping is a nonlinear equations (Jacobian) matrix \see SolverMappingRoutines_EquationsMatrixTypes,SolverMappingRoutines
  !>@}
 
  !> \addtogroup SolverMappingRoutines_EquationsTypes SolverMappingRoutines::EquationsTypes
  !> \brief Equations Matrix types
  !> \see SolverMappingRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_MAPPING_EQUATIONS_EQUATIONS_SET=1 !<The equations in the solver mapping is from an equations set \see SolverMappingRoutines_EquationsTypes,SolverMappingRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_MAPPING_EQUATIONS_INTERFACE_CONDITION=2 !<The equations in the solver mapping is from an interface condition \see SolverMappingRoutines_EquationsTypes,SolverMappingRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_MAPPING_EQUATIONS_INTERFACE_TRANSPOSE=3 !<The equations in the solver mapping is from a transposed interface condition \see SolverMappingRoutines_EquationsTypes,SolverMappingRoutines
  !>@}
 
  !Module types

  !Module variables 

  !Interfaces

  PUBLIC SOLVER_MAPPING_EQUATIONS_DYNAMIC_MATRIX,SOLVER_MAPPING_EQUATIONS_LINEAR_MATRIX

  PUBLIC SOLVER_MAPPING_EQUATIONS_EQUATIONS_SET,SOLVER_MAPPING_EQUATIONS_INTERFACE_CONDITION

  PUBLIC SolverMapping_AssertIsFinished,SolverMapping_AssertNotFinished

  PUBLIC SolverMapping_ColumnDOFSMappingGet

  PUBLIC SolverMapping_CreateValuesCacheGet
  
  PUBLIC SolverMapping_DynamicEquationsMatrixToSolverMatrixGet
  
  PUBLIC SolverMapping_EMSToSMMapNumberOfDynamicMatricesGet

  PUBLIC SolverMapping_EMSToSMMapNumberOfJacobianMatricesGet

  PUBLIC SolverMapping_EMSToSMMapNumberOfLinearMatricesGet

  PUBLIC SolverMapping_EMSToSMMapNumberOfVariablesGet

  PUBLIC SolverMapping_EMSToSMMapVariableGet

  PUBLIC SolverMapping_EMSToSMMapVariableToSolverColMapGet

  PUBLIC SolverMapping_EquationsSetGet

  PUBLIC SolverMapping_EquationsRowToSolverRowsMapGet

  PUBLIC SolverMapping_EquationsMatrixToSolverMatrixMapGet
  
  PUBLIC SolverMapping_EquationsSetToSolverMatricesMapGet

  PUBLIC SolverMapping_IMSToSMMapLagrangeVariableToSolverColMapGet

  PUBLIC SolverMapping_IMSToSMMapDependentVariableToSolverColMapGet

  PUBLIC SolverMapping_InterfaceColToSolverColsMapGet

  PUBLIC SolverMapping_InterfaceColToSolverRowsMapGet

  PUBLIC SolverMapping_InterfaceConditionGet

  PUBLIC SolverMapping_InterfaceConditionToSolverMatricesMapGet

  PUBLIC SolverMapping_InterfaceEquationsMatrixToSolverMatrixGet
  
  PUBLIC SolverMapping_InterfaceMatrixToSolverMatrixMapGet 

  PUBLIC SolverMapping_InterfaceRowToSolverRowsMapGet

  PUBLIC SolverMapping_JacobianEquationsMatrixToSolverMatrixGet
  
  PUBLIC SolverMapping_JacobianMatrixToSolverMatrixMapGet

  PUBLIC SolverMapping_LinearEquationsMatrixToSolverMatrixGet
  
  PUBLIC SolverMapping_RowDOFSMappingGet

  PUBLIC SolverMapping_SolverColToEquationsColMapGet

  PUBLIC SolverMapping_SolverMatrixToEquationsMapGet

  PUBLIC SolverMappingCreateValuesCache_EquationsVariableListGet
  
  PUBLIC SolverMappingCreateValuesCache_InterfaceVariableListGet
  
  PUBLIC SolverMappingCreateValuesCache_RHSVariableListGet
  
  PUBLIC SolverMappingEMSToSMMap_DynamicEquationsMatrixToSolverMatrixMapGet

  PUBLIC SolverMappingEMSToSMMap_JacobianEquationsMatrixToSolverMatrixMapGet

  PUBLIC SolverMappingEMSToSMMap_LinearEquationsMatrixToSolverMatrixMapGet

  PUBLIC SolverMappingEMSToSMMap_NumberOfDynamicMatricesGet

  PUBLIC SolverMappingEMSToSMMap_NumberOfJacobianMatricesGet

  PUBLIC SolverMappingEMSToSMMap_NumberOfLinearMatricesGet

  PUBLIC SolverMappingEMToSMMap_EquationsColToSolverColsMapGet

  PUBLIC SolverMappingEMToSMMap_EquationsMatrixGet

  PUBLIC SolverMappingEMToSMMap_SolverMatrixGet

  PUBLIC SolverMappingEMToSMMap_VariableGet
  
  PUBLIC SolverMappingESToSMSMap_EquationsMatricesToSolverMatrixMapGet

  PUBLIC SolverMappingESToSMSMap_EquationsMatrixToSolverMatricesMapGet

  PUBLIC SolverMappingESToSMSMap_JacobianMatrixToSolverMatrixMapGet

  PUBLIC SolverMappingESToSMSMap_NumberOfEquationsMatricesGet
  
  PUBLIC SolverMappingESToSMSMap_NumberOfJacobianMatricesGet

  PUBLIC SolverMappingESToSMSMap_NumberOfSolverMatricesGet

  PUBLIC SolverMappingICToSMSMap_InterfaceMatricesToSolverMatrixMapGet

  PUBLIC SolverMappingIMSToSMMap_InterfaceMatrixToSolverMatrixMapGet
  
  PUBLIC SolverMappingIMSToSMMap_NumberOfInterfaceMatricesGet
  
  PUBLIC SolverMappingIMToSMMap_InterfaceMatrixGet

  PUBLIC SolverMappingIMToSMMap_InterfaceRowToSolverColsMapGet

  PUBLIC SolverMappingIMToSMMap_SolverMatrixGet

  PUBLIC SolverMappingJMToSMMap_JacobianColToSolverColsMapGet

  PUBLIC SolverMappingJMToSMMap_JacobianMatrixGet
  
  PUBLIC SolverMappingJMToSMMap_SolverMatrixGet

CONTAINS

  !
  !================================================================================================================================
  !

  !>Assert that an solver mapping has been finished
  SUBROUTINE SolverMapping_AssertIsFinished(solverMapping,err,error,*)

    !Argument Variables
    TYPE(SolverMappingType), POINTER, INTENT(IN) :: solverMapping !<The solver mapping to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMapping_AssertIsFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solverMapping)) CALL FlagError("Solver mapping is not associated.",err,error,*999)
#endif    

    IF(.NOT.solverMapping%solverMappingFinished) CALL FlagError("Solver mapping has not been finished.",err,error,*999)
    
    EXITS("SolverMapping_AssertIsFinished")
    RETURN
999 ERRORSEXITS("SolverMapping_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_AssertIsFinished

  !
  !================================================================================================================================
  !

  !>Assert that an solver mapping has not been finished
  SUBROUTINE SolverMapping_AssertNotFinished(solverMapping,err,error,*)

    !Argument Variables
    TYPE(SolverMappingType), POINTER, INTENT(IN) :: solverMapping !<The solver mapping to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMapping_AssertNotFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solverMapping)) CALL FlagError("Solver mapping is not associated.",err,error,*999)
#endif    

    IF(solverMapping%solverMappingFinished) CALL FlagError("Solver mapping has already been finished.",err,error,*999)
    
    EXITS("SolverMapping_AssertNotFinished")
    RETURN
999 ERRORSEXITS("SolverMapping_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_AssertNotFinished

  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the column DOFs mapping for a solver matrix index for solver mapping.
  SUBROUTINE SolverMapping_ColumnDOFSMappingGet(solverMapping,matrixIndex,columnDOFSMapping,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the column DOFs mapping for
    INTEGER(INTG), INTENT(IN) :: matrixIndex !<The matrix index in the solver mapping to get the column DOFs mapping for
    TYPE(DomainMappingType), POINTER :: columnDOFSMapping !<On exit, a pointer to the specified column DOFs mapping. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMapping_ColumnDOFSMappingGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(columnDOFSMapping)) CALL FlagError("Column DOFs mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMapping)) CALL FlagError("Solver mapping is not associated.",err,error,*999)
    IF(matrixIndex<1.OR.matrixIndex>solverMapping%numberOfSolverMatrices) THEN
      localError="The specified matrix index of "//TRIM(NumberToVString(matrixIndex,"*",err,error))// &
        & " is invalid. The index must be >= 1 and <= "// &
        & TRIM(NumberToVString(solverMapping%numberOfSolverMatrices,"*",err,error))//"."      
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(ALLOCATED(solverMapping%solverColToEquationsColsMap)) &
      & CALL FlagError("The solver column to equations columns map is not allocated for the solver mapping.",err,error,*999)
#endif    

    columnDOFSMapping=>solverMapping%solverColToEquationsColsMap(matrixIndex)%columnDOFSMapping

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(columnDOFSMapping)) THEN
      localError="The column DOFs mapping for the specified matrix index of "// &
        & TRIM(NumberToVString(matrixIndex,"*",err,error))//" is not associated."      
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("SolverMapping_ColumnDOFSMappingGet")
    RETURN
999 NULLIFY(columnDOFSMapping)
998 ERRORSEXITS("SolverMapping_ColumnDOFSMappingGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_ColumnDOFSMappingGet
  
  !
  !================================================================================================================================
  !

  !>Gets the create values cache from a solver mapping.
  SUBROUTINE SolverMapping_CreateValuesCacheGet(solverMapping,createValuesCache,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<The solver mapping to get the create values cache for.
    TYPE(SolverMappingCreateValuesCacheType), POINTER :: createValuesCache !<On exit, a pointer to the create values cache for the solver mapping. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMapping_CreateValuesCacheGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(createValuesCache)) CALL FlagError("Create values cache is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMapping)) CALL FlagError("Solver mapping is not associated.",err,error,*999)
#endif

    createValuesCache=>solverMapping%createValuesCache

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(createValuesCache)) &
      & CALL FlagError("Create values cache is not associated for solver mapping.",err,error,*999)
#endif    
    
    EXITS("SolverMapping_CreateValuesCacheGet")
    RETURN
999 NULLIFY(createValuesCache)
998 ERRORSEXITS("SolverMapping_CreateValuesCacheGet",err,error)
    RETURN 1    
    
  END SUBROUTINE SolverMapping_CreateValuesCacheGet

  !
  !================================================================================================================================
  !
  
  !>Returns a dynamic equations matrix in an equations set mapped to a solver matrix in a solver mapping.
  SUBROUTINE SolverMapping_DynamicEquationsMatrixToSolverMatrixGet(solverMapping,solverMatrixIndex,equationsSetIdx, &
    & dynamicMatrixIdx,dynamicEquationsMatrix,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the dynamic matrix for
    INTEGER(INTG), INTENT(IN) :: solverMatrixIndex !<The solver matrix index in the solver mapping to get the number of dynamic matrices for
    INTEGER(INTG), INTENT(IN) :: equationsSetIdx !<The equations set index containing the dynamic matrices mapped to the solver matrix
    INTEGER(INTG), INTENT(IN) :: dynamicMatrixIdx !<The dynamic matrix index of the equations matrix to the solver matrix
    TYPE(EquationsMatrixIdx), POINTER :: dynamicMatrix !<On return the specified dynamic matrix. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(DynamicEquationsMatrixToSolverMatrixMapType), POINTER :: dynamicEquationsMatrixToSolverMatrixMap
    TYPE(EquationsMatricesToSolverMatrixMapType), POINTER :: equationsMatricesToSolverMatrixMap
    TYPE(EquationsSetToSolverMatricesMapType), POINTER :: equationsSetToSolverMatricesMap
#endif
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMapping_DynamicEquationsMatrixToSolverMatrixGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dynamicMatrix)) CALL FlagError("Dynamic equations matrix is already associated.",err,eror,*998)
    NULLIFY(equationsSetToSolverMatricesMap)
    CALL SolverMapping_EquationsSetToSolverMatricesMapGet(solverMapping,equationsSetIdx,equationsSetToSolverMatricesMap, &
      & err,error,*999)
    NULLIFY(equationsMatricesToSolverMatrixMap)
    CALL SolverMappingESToSMSMap_EquationsMatricesToSolverMatrixMapGet(equationsSetToSolverMatricesMap,solverMatrixIdx, &
      & equationsMatricesToSolverMatrixMap,err,error,*999)
    NULLIFY(dynamicEquationsMatrixToSolverMatrixMap)
    CALL SolverMappingEMSToSMMap_DynamicEquationsMatrixToSolverMatrixMapGet(equationsMatricesToSolverMatrixMap, &
      & dynamicMatrixIdx,dynamicEquationsMatrixToSolverMatrixMap,err,error,*999)
#endif    

    dynamicEquationsMatrix=>solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
      & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
      & dynamicEquationsMatrixToSolverMatrixMaps(dynamicMatrixIdx)%ptr%equationsMatrix

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(dynamicEquationsMatrix)) THEN
      localError="The dynamic equations matrix is not associated for the dynamic matrix index of "// &
        & TRIM(NumberToVString(dynamicMatrixIdx,"*",err,error))// &
        & " of the dynamic equations matrix to solver matrix maps for the solver matrix index of "// &
        & TRIM(NumberToVString(solverMatrixIndex,"*",err,error))// &
        & " of the equations matrices to solver matrix maps for the equations set index "// &
        & TRIM(NumberToVString(equationsSetIdx,"*",err,error))// &
        & " of the equations set to solver matrices map of the solver mapping."
      CALL FlagError(localError,err,error,*999)    
    ENDIF
#endif    
      
    EXITS("SolverMapping_DynamicEquationsMatrixToSolverMatrixGet")
    RETURN
999 NULLIFY(dynamicEquationsMatrix)
998 ERRORSEXITS("SolverMapping_DynamicEquationsMatrixToSolverMatrixGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_DynamicEquationsMatrixToSolverMatrixGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of dynamic matrices in an equations set mapped to a solver matrix in a solver mapping.
  SUBROUTINE SolverMapping_EMSToSMMapNumberOfDynamicMatricesGet(solverMapping,solverMatrixIndex,equationsSetIdx, &
    & numberOfDynamicMatrices,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the number of dynamic matrices for
    INTEGER(INTG), INTENT(IN) :: solverMatrixIndex !<The solver matrix index in the solver mapping to get the number of dynamic matrices for
    INTEGER(INTG), INTENT(IN) :: equationsSetIdx !<The equations set index containing the dynamic matrices mapped to the solver matrix
    INTEGER(INTG), INTENT(OUT) :: numberOfDynamicMatrices !<On exit, the number of dynamic matrices mapped to the solver matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(EquationsMatricesToSolverMatrixMapType), POINTER :: equationsMatricesToSolverMatrixMap
    TYPE(EquationsSetToSolverMatricesMapType), POINTER :: equationsSetToSolverMatricesMap
#endif
 
    ENTERS("SolverMapping_EMSToSMMapNumberOfDynamicMatricesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    NULLIFY(equationsSetToSolverMatricesMap)
    CALL SolverMapping_EquationsSetToSolverMatricesMapGet(solverMapping,equationsSetIdx,equationsSetToSolverMatricesMap, &
      & err,error,*999)
    NULLIFY(equationsMatricesToSolverMatrixMap)
    CALL SolverMappingESToSMSMap_EquationsMatricesToSolverMatrixMapGet(equationsSetToSolverMatricesMap,solverMatrixIdx, &
      & equationsMatricesToSolverMatrixMap,err,error,*999)
#endif    
    
    numberOfDynamicMatrices=solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
      & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%numberOfDynamicEquationsMatrices
      
    EXITS("SolverMapping_EMSToSMMapNumberOfDynamicMatricesGet")
    RETURN
999 ERRORSEXITS("SolverMapping_EMSToSMMapNumberOfDynamicMatricesGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_EMSToSMMapNumberOfDynamicMatricesGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of interface matrices in an interface condition mapped to a solver matrix in a solver mapping.
  SUBROUTINE SolverMapping_EMSToSMMapNumberOfInterfaceMatricesGet(solverMapping,solverMatrixIndex,interfaceConditionIdx, &
    & numberOfInterfaceMatrices,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the number of interface matrices for
    INTEGER(INTG), INTENT(IN) :: solverMatrixIndex !<The solver matrix index in the solver mapping to get the number of interface matrices for
    INTEGER(INTG), INTENT(IN) :: interfaceConditionIdx !<The interface condition index containing the interface matrices mapped to the solver matrix
    INTEGER(INTG), INTENT(OUT) :: numberOfInterfaceMatrices !<On exit, the number of interface matrices mapped to the solver matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(EquationsMatricesToSolverMatrixMapType), POINTER :: equationsMatricesToSolverMatrixMap
    TYPE(EquationsSetToSolverMatricesMapType), POINTER :: equationsSetToSolverMatricesMap
#endif
 
    ENTERS("SolverMapping_EMSToSMMapNumberOfInterfaceMatricesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    NULLIFY(equationsSetToSolverMatricesMap)
    CALL SolverMapping_EquationsSetToSolverMatricesMapGet(solverMapping,equationsSetIdx,equationsSetToSolverMatricesMap, &
      & err,error,*999)
    NULLIFY(equationsMatricesToSolverMatrixMap)
    CALL SolverMappingESToSMSMap_EquationsMatricesToSolverMatrixMapGet(equationsSetToSolverMatricesMap,solverMatrixIdx, &
      & equationsMatricesToSolverMatrixMap,err,error,*999)
#endif    
    
    numberOfInterfaceMatrices=solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
      & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%numberOfInterfaceMatrices
      
    EXITS("SolverMapping_EMSToSMMapNumberOfInterfaceMatricesGet")
    RETURN
999 ERRORSEXITS("SolverMapping_EMSToSMMapNumberOfInterfaceMatricesGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_EMSToSMMapNumberOfInterfaceMatricesGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of Jacobian matrices in an equations set mapped to a solver matrix in a solver mapping.
  SUBROUTINE SolverMapping_EMSToSMMapNumberOfJacobianMatricesGet(solverMapping,solverMatrixIndex,equationsSetIdx, &
    & numberOfJacobianMatrices,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the number of Jacobian matrices for
    INTEGER(INTG), INTENT(IN) :: solverMatrixIndex !<The solver matrix index in the solver mapping to get the number of Jacobian matrices for
    INTEGER(INTG), INTENT(IN) :: equationsSetIdx !<The equations set index containing the Jacobian matrices mapped to the solver matrix
    INTEGER(INTG), INTENT(OUT) :: numberOfJacobianMatrices !<On exit, the number of Jacobian matrices mapped to the solver matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(EquationsMatricesToSolverMatrixMapType), POINTER :: equationsMatricesToSolverMatrixMap
    TYPE(EquationsSetToSolverMatricesMapType), POINTER :: equationsSetToSolverMatricesMap
#endif
 
    ENTERS("SolverMapping_EMSToSMMapNumberOfJacobianMatricesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    NULLIFY(equationsSetToSolverMatricesMap)
    CALL SolverMapping_EquationsSetToSolverMatricesMapGet(solverMapping,equationsSetIdx,equationsSetToSolverMatricesMap, &
      & err,error,*999)
    NULLIFY(equationsMatricesToSolverMatrixMap)
    CALL SolverMappingESToSMSMap_EquationsMatricesToSolverMatrixMapGet(equationsSetToSolverMatricesMap,solverMatrixIdx, &
      & equationsMatricesToSolverMatrixMap,err,error,*999)
#endif    
    
    numberOfJacobianMatrices=solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
      & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%numberOfEquationsJacobians
      
    EXITS("SolverMapping_EMSToSMMapNumberOfJacobianMatricesGet")
    RETURN
999 ERRORSEXITS("SolverMapping_EMSToSMMapNumberOfJacobianMatricesGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_EMSToSMMapNumberOfJacobianMatricesGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of linear matrices in an equations set mapped to a solver matrix in a solver mapping.
  SUBROUTINE SolverMapping_EMSToSMMapNumberOfLinearMatricesGet(solverMapping,solverMatrixIndex,equationsSetIdx, &
    & numberOfLinearMatrices,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the number of linear matrices for
    INTEGER(INTG), INTENT(IN) :: solverMatrixIndex !<The solver matrix index in the solver mapping to get the number of linear matrices for
    INTEGER(INTG), INTENT(IN) :: equationsSetIdx !<The equations set index containing the linear matrices mapped to the solver matrix
    INTEGER(INTG), INTENT(OUT) :: numberOfLinearMatrices !<On exit, the number of linear matrices mapped to the solver matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(EquationsMatricesToSolverMatrixMapType), POINTER :: equationsMatricesToSolverMatrixMap
    TYPE(EquationsSetToSolverMatricesMapType), POINTER :: equationsSetToSolverMatricesMap
#endif
 
    ENTERS("SolverMapping_EMSToSMMapNumberOfLinearMatricesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    NULLIFY(equationsSetToSolverMatricesMap)
    CALL SolverMapping_EquationsSetToSolverMatricesMapGet(solverMapping,equationsSetIdx,equationsSetToSolverMatricesMap, &
      & err,error,*999)
    NULLIFY(equationsMatricesToSolverMatrixMap)
    CALL SolverMappingESToSMSMap_EquationsMatricesToSolverMatrixMapGet(equationsSetToSolverMatricesMap,solverMatrixIdx, &
      & equationsMatricesToSolverMatrixMap,err,error,*999)
#endif    
    
    numberOfLinearMatrices=solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
      & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%numberOfLinearEquationsMatrices
      
    EXITS("SolverMapping_EMSToSMMapNumberOfLinearMatricesGet")
    RETURN
999 ERRORSEXITS("SolverMapping_EMSToSMMapNumberOfLinearMatricesGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_EMSToSMMapNumberOfLinearMatricesGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of variables in an equations set mapped to a solver matrix in a solver mapping.
  SUBROUTINE SolverMapping_EMSToSMMapNumberOfVariablesGet(solverMapping,solverMatrixIndex,equationsSetIdx,numberOfVariables, &
    & err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the number of variables for
    INTEGER(INTG), INTENT(IN) :: solverMatrixIndex !<The solver matrix index in the solver mapping to get the number of variables for
    INTEGER(INTG), INTENT(IN) :: equationsSetIdx !<The equations set index containing the variables mapped to the solver matrix
    INTEGER(INTG), INTENT(OUT) :: numberOfVariables !<On exit, the number of variables mapped to the solver matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(EquationsMatricesToSolverMatrixMapType), POINTER :: equationsMatricesToSolverMatrixMap
    TYPE(EquationsSetToSolverMatricesMapType), POINTER :: equationsSetToSolverMatricesMap
#endif
 
    ENTERS("SolverMapping_EMSToSMMapNumberOfVariablesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    NULLIFY(equationsSetToSolverMatricesMap)
    CALL SolverMapping_EquationsSetToSolverMatricesMapGet(solverMapping,equationsSetIdx,equationsSetToSolverMatricesMap, &
      & err,error,*999)
    NULLIFY(equationsMatricesToSolverMatrixMap)
    CALL SolverMappingESToSMSMap_EquationsMatricesToSolverMatrixMapGet(equationsSetToSolverMatricesMap,solverMatrixIdx, &
      & equationsMatricesToSolverMatrixMap,err,error,*999)
#endif    
    
    numberOfVariables=solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
      & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%numberOfVariables
      
    EXITS("SolverMapping_EMSToSMMapNumberOfVariablesGet")
    RETURN
999 ERRORSEXITS("SolverMapping_EMSToSMMapNumberOfVariablesGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_EMSToSMMapNumberOfVariablesGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the variable in an equations set mapped to a solver matrix in a solver mapping.
  SUBROUTINE SolverMapping_EMSToSMMapVariableGet(solverMapping,solverMatrixIndex,equationsSetIdx,variableIdx,fieldVariable, &
    & err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the variable for
    INTEGER(INTG), INTENT(IN) :: solverMatrixIndex !<The solver matrix index in the solver mapping to get the variable for
    INTEGER(INTG), INTENT(IN) :: equationsSetIdx !<The equations set index containing the variables mapped to the solver matrix
    INTEGER(INTG), INTENT(IN) :: variableIdx !<The variable index in the equations set mapped to the solver matrix
    TYPE(FieldVariableType), POINTER :: fieldVariable !<On exit, a pointer to the specified field variable mapped to the solver matrix. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(EquationsMatricesToSolverMatrixMapType), POINTER :: equationsMatricesToSolverMatrixMap
    TYPE(EquationsSetToSolverMatricesMapType), POINTER :: equationsSetToSolverMatricesMap
#endif
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMapping_EMSToSMMapVariableGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is already associated.",err,error,*998)
    NULLIFY(equationsSetToSolverMatricesMap)
    CALL SolverMapping_EquationsSetToSolverMatricesMapGet(solverMapping,equationsSetIdx,equationsSetToSolverMatricesMap, &
      & err,error,*999)
    NULLIFY(equationsMatricesToSolverMatrixMap)
    CALL SolverMappingESToSMSMap_EquationsMatricesToSolverMatrixMapGet(equationsSetToSolverMatricesMap,solverMatrixIdx, &
      & equationsMatricesToSolverMatrixMap,err,error,*999)
    IF(variableIdx<1.OR.variableIdx>equationsMatricesToSolverMatrixMap%numberOfVariables) THEN
      localError="The specified variable index of "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
        & " is invalid for the solver matrix index of "//TRIM(NumberToVString(solverMatrixIndex,"*",err,error))// &
        & " of the equations matrices to solver matrix maps for the equations set index "// &
        & TRIM(NumberToVString(equationsSetIdx,"*",err,error))// &
        & " of the equations set to solver matrices map of the solver mapping. The variable index should be >= 1 and <= "// &
        & TRIM(NumberToVString(equationsMatricesToSolverMatrixMap%numberOfVariables,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(equationsMatricesToSolverMatrixMap%variables)) THEN
      localError="The variables array is not allocated for the solver matrix index of "// &
        & TRIM(NumberToVString(solverMatrixIndex,"*",err,error))// &
        & " of the equations matrices to solver matrix maps sm the equations set index "// &
        & TRIM(NumberToVString(equationsSetIdx,"*",err,error))// &
        & " of the equations set to solver matrices map of the solver mapping."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
     
    fieldVariable=>solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
      & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%variables(variableIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(fieldVariable)) THEN
      localError="The field variable is not associated for the specified variable index of "// &
        & TRIM(NumberToVString(variableIdx,"*",err,error))//" for the solver matrix index of "// &
        & TRIM(NumberToVString(solverMatrixIndex,"*",err,error))// &
        & " of the equations to solver matrix maps sm for the equations set index "// &
        & TRIM(NumberToVString(equationsSetIdx,"*",err,error))// &
        & " of the equations set to solver map of the solver mapping."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("SolverMapping_EMSToSMMapVariableGet")
    RETURN
999 NULLIFY(fieldVariable)
998 ERRORSEXITS("SolverMapping_EMSToSMMapVariableGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_EMSToSMMapVariableGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the variable to solver col map in an equations set mapped to a solver matrix in a solver mapping.
  SUBROUTINE SolverMapping_EMSToSMMapVariableToSolverColMapGet(solverMapping,solverMatrixIndex,equationsSetIdx,variableIdx, &
    & variableToSolverColMap,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the variable to solver col map for
    INTEGER(INTG), INTENT(IN) :: solverMatrixIndex !<The solver matrix index in the solver mapping to get the variable to solver col map for
    INTEGER(INTG), INTENT(IN) :: equationsSetIdx !<The equations set index containing the variables mapped to the solver matrix
    INTEGER(INTG), INTENT(IN) :: variableIdx !<The variable index in the equations set mapped to the solver matrix
    TYPE(VariableToSolverColMapType), POINTER :: variableToSolverColMap !<On exit, a pointer to the specified variable to solver col map for the solver matrix. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(EquationsMatricesToSolverMatrixMapType), POINTER :: equationsMatricesToSolverMatrixMap
    TYPE(EquationsSetToSolverMatricesMapType), POINTER :: equationsSetToSolverMatricesMap
#endif
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMapping_EMSToSMMapVariableToSolverColMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(variableToSolverColMap)) CALL FlagError("Variable to solver col map is already associated.",err,error,*998)
    NULLIFY(equationsSetToSolverMatricesMap)
    CALL SolverMapping_EquationsSetToSolverMatricesMapGet(solverMapping,equationsSetIdx,equationsSetToSolverMatricesMap, &
      & err,error,*999)
    NULLIFY(equationsMatricesToSolverMatrixMap)
    CALL SolverMappingESToSMSMap_EquationsMatricesToSolverMatrixMapGet(equationsSetToSolverMatricesMap,solverMatrixIdx, &
      & equationsMatricesToSolverMatrixMap,err,error,*999)
    IF(variableIdx<1.OR.variableIdx>equationsMatricesToSolverMatrixMap%numberOfVariables) THEN
      localError="The specified variable index of "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
        & " is invalid for the solver matrix index of "//TRIM(NumberToVString(solverMatrixIndex,"*",err,error))// &
        & " of the equations matrices to solver matrix maps for the equations set index "// &
        & TRIM(NumberToVString(equationsSetIdx,"*",err,error))// &
        & " of the equations set to solver matrices map of the solver mapping. The variable index should be >= 1 and <= "// &
        & TRIM(NumberToVString(equationsMatricesToSolverMatrixMap%numberOfVariables,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(equationsMatricesToSolverMatrixMap%variableToSolverColMaps)) THEN
      localError="The variable to solver col maps array is not allocated for the solver matrix index of "// &
        & TRIM(NumberToVString(solverMatrixIndex,"*",err,error))// &
        & " of the equations matrices to solver matrix maps sm the equations set index "// &
        & TRIM(NumberToVString(equationsSetIdx,"*",err,error))// &
        & " of the equations set to solver matrices map of the solver mapping."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
     
   variableToSolverColMap=>solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
     & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
     & variableToSolverColMaps(variableIdx)%ptr

#ifdef WITH_POSTCHECKS   
    IF(.NOT.ASSOCIATED(variableToSolverColMap)) THEN
      localError="The variable to solver col map is not associated for the specified variable index of "// &
        & TRIM(NumberToVString(variableIdx,"*",err,error))//" for the solver matrix index of "// &
        & TRIM(NumberToVString(solverMatrixIndex,"*",err,error))// &
        & " of the equations matrices to solver matrix maps for the equations set index "// &
        & TRIM(NumberToVString(equationsSetIdx,"*",err,error))// &
        & " of the equations set to solver matrices map of the solver mapping."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("SolverMapping_EMSToSMMapVariableToSolverColMapGet")
    RETURN
999 NULLIFY(variableToSolverColMap)
998 ERRORSEXITS("SolverMapping_EMSToSMMapVariableToSolverColMapGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_EMSToSMMapVariableToSolverColMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the equations set for solver mapping.
  SUBROUTINE SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the equations set for
    INTEGER(INTG), INTENT(IN) :: equationsSetIdx !<The equations set index in the solver mapping to get the equations set for
    TYPE(EquationsSetType), POINTER :: equationsSet !<On exit, a pointer to the specified equations set. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMapping_EquationsSetGet",err,error,*998)

#ifdef WITH_PRECHECKS
    IF(ASSOCIATED(equationsSet)) CALL FlagError("Equations set is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMapping)) CALL FlagError("Solver mapping is not associated.",err,error,*999)
    IF(equationsSetIdx<1.OR.equationsSetIdx>solverMapping%numberOfEquationsSets) THEN
      localError="The specified equations set index of "//TRIM(NumberToVString(equationsSetIdx,"*",err,error))// &
        & " is invalid. The index must be >= 1 and <= "// &
        & TRIM(NumberToVString(SolverMapping%numberOfEquationsSets,"*",err,error))//"."      
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(solverMapping%equationsSets)) &
      & CALL FlagError("Solver mapping equations sets is not allocated.",err,error,*999)
#endif    

    equationsSet=>solverMapping%equationsSets(equationsSetIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(equationsSet)) THEN
      localError="The equations set for the specified equations set index of "// &
        & TRIM(NumberToVString(equationsSetIdx,"*",err,error))//" is not associated."      
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("SolverMapping_EquationsSetGet")
    RETURN
999 NULLIFY(equationsSet)
998 ERRORSEXITS("SolverMapping_EquationsSetGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_EquationsSetGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the equations row to solver rows map for solver mapping.
  SUBROUTINE SolverMapping_EquationsRowToSolverRowsMapGet(solverMapping,equationsSetIdx,equationsRowToSolverRowsMap,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the equations row to solver rows map for
    INTEGER(INTG), INTENT(IN) :: equationsSetIdx !<The equations set index in the solver mapping to get the equations row to solver rows map for
    TYPE(MatrixRowColCouplingType), POINTER :: equationsRowToSolverRowsMap(:) !<On exit, a pointer to the specified equations row to solver rows map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(EquationsSetToSolverMatricesMapType), POINTER :: equationsSetToSolverMatricesMap    
#endif        
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMapping_EquationsRowToSolverRowsMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(equationsRowToSolverRowsMap)) &
      & CALL FlagError("Equations row to solver rows map is already associated.",err,error,*998)
    NULLIFY(equationsSetToSolverMatricesMap)
    CALL SolverMapping_EquationsSetToSolverMatricesMapGet(solverMapping,equationsSetIdx,equationsSetToSolverMatricesMap, &
      & err,error,*999)
#endif
    
    equationsRowToSolverRowsMap=>solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr%equationsRowToSolverRowsMap

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(equationsRowToSolverRowsMap)) THEN
      localError="The equations row to solver rows map for the specified equations set index of "// &
        & TRIM(NumberToVString(equationsSetIdx,"*",err,error))//" is not associated for the solver mapping."      
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("SolverMapping_EquationsRowToSolverRowsMapGet")
    RETURN
999 NULLIFY(equationsRowToSolverRowsMap)
998 ERRORS("SolverMapping_EquationsRowToSolverRowsMapGet",err,error)
    EXITS("SolverMapping_EquationsRowToSolverRowsMapGet")
    RETURN 1
    
  END SUBROUTINE SolverMapping_EquationsRowToSolverRowsMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the equations matrix to solver matrix map for a solver mapping.
  SUBROUTINE SolverMapping_EquationsMatrixToSolverMatrixMapGet(solverMapping,equationsSetIndex,equationsMatrixIndex, &
    & solverMatrixIndex,equationsMatrixToSolverMatrixMap,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the equations matrix to solver matrix map for
    INTEGER(INTG), INTENT(IN) :: equationsSetIndex !<The equations set index in the solver mapping to get the equations matrix to solver matrix map for
    INTEGER(INTG), INTENT(IN) :: equationsMatrixIdx !<The equations matrix index for the equations matrix to solver matrix map
    INTEGER(INTG), INTENT(IN) :: solverMatrixIdx !<The solver matrix index for the equations matrix to solver matrix map
    TYPE(EquationsMatrixToSolverMatrixMapType), POINTER :: equationsMatrixToSolverMatrixMap !<On exit, a pointer to the specified equations matrix to solver matrix map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(EquationsMatrixToSolverMatricesMapType), POINTER :: equationsMatrixToSolverMatricesMap
   TYPE(EquationsSetToSolverMatricesMapType), POINTER :: equationsSetToSolverMatricesMap    
#endif
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMapping_EquationsMatrixToSolverMatrixMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(equationsMatrixToSolverMatrixMap)) &
      & CALL FlagError("Equations matrix to solver matrix map is already associated.",err,error,*998)
    NULLIFY(equationsSetToSolverMatricesMap)
    CALL SolverMapping_EquationsSetToSolverMatricesMapGet(solverMapping,equationsSet,equationsSetToSolverMatricesMap, &
      & err,error,*999)
    NULLIFY(equationsMatrixToSolverMatricesMap)
    CALL SolverMappingESToSMSMap_EquationsMatrixToSolverMatricesMapGet(equationsSetToSolverMatricesMap,equationsMatrixIndex, &
      & equationsMatrixToSolverMatricesMap,err,error,*999)
    IF(solverMatrixIdx<1.OR.solverMatrixIdx>equationsMatrixToSolverMatricesMap%numberOfSolverMatrices) THEN
      localError="The specified solver matrix index of "//TRIM(NumberToVString(solverMatrixIdx,"*",err,error))//&
        & " is invalid for the equations matrix index of "//TRIM(NumberToVString(equationsMatrixIdx,"*",err,error))//&
        & " of the equations matrix to solver matrices maps for the equations set index of "// &
        & TRIM(NumberToVString(equationsSetIdx,"*",err,error))//" of the equations set to solver matrices maps "// &
        & "of the solver mapping. The solver matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(equationsMatrixToSolverMatricesMap%numberOfSolverMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)      
    ENDIF
    IF(.NOT.ALLOCATED(equationsMatrixToSolverMatriceMap%equationsMatrixToSolverMatrixMaps)) THEN
      localError="The equations matrix to solver matrix maps is not allocated for the equations matrix index of "// &
        & TRIM(NumberToVString(equationsMatrixIdx,"*",err,error))//&
        & " of the equations matrix to solver matrices maps for the equations set index of "// &
        & TRIM(NumberToVString(equationsSetIdx,"*",err,error))//" of the equations set to solver matrices maps "// &
        & "of the solver mapping."
    ENDIF
 #endif    
    
    equationsMatrixToSolverMatrixMap=>solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
      & equationsMatrixToSolverMatricesMaps(equationsMatrixIdx)%ptr% &
      & equationsMatrixToSolverMatrixMaps(solverMatrixIndex)%ptr
    
#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(equationsMatrixToSolverMatrixMap)) THEN
      localError="The equations matrix to solver matrix map is not associated for the solver matrix index of "// &
        & TRIM(NumberToVString(solverMatrixIndex,"*",err,error))// &
        & " of the equations matrix to solver matrix maps of the equations matrix index of "// &
        & TRIM(NumberToVString(equationsMatrixIndex,"*",err,error))// &
        & " of the equations matrix to solver matrices maps for the equations set index of "// &
        & TRIM(NumberToVString(equationsSetIndex,"*",err,error))// &
        & " of the equations set to solver matrices maps of the solver mapping."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    EXITS("SolverMapping_EquationsMatrixToSolverMatrixMapGet")
    RETURN
999 NULLIFY(equationsMatrixToSolverMatrixMap)
998 ERRORSEXITS("SolverMapping_EquationsMatrixToSolverMatrixMapGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_EquationsMatrixToSolverMatrixMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the equations set to solver matrices map for a solver mapping.
  SUBROUTINE SolverMapping_EquationsSetToSolverMatricesMapGet(solverMapping,equationsSetIdx,equationsSetToSolverMatricesMap, &
    & err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the equations set to solver matrices map for
    INTEGER(INTG), INTENT(IN) :: equationsSetIdx !<The equations set index in the solver mapping to get the equations set to solver matrices map for
    TYPE(EquationsSetToSolverMatricesMapType), POINTER :: equationsSetToSolverMatricesMap !<On exit, a pointer to the specified equations set to solver matrices map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMapping_EquationsSetToSolverMatricesMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(equationsSetToSolverMatricesMap)) &
      & CALL FlagError("Equations set to solver matrices map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMapping)) CALL FlagError("Solver mapping is not associated.",err,error,*999)
    IF(equationsSetIdx1.OR.equationsSetIdx>solverMapping%numberOfEquationsSets) THEN
      localError="The specified equations set index of "//TRIM(NumberToVString(equationsSetIdx,"*",err,error))// &
        & " is invalid for the solver mapping. The index must be >= 1 and <= "// &
        & TRIM(NumberToVString(solverMapping%numberOfEquationsSets,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(solverMapping%equationsSetToSolverMatricesMaps)) &
      & CALL FlagError("Solver mapping equations sets to solver matrices maps is not allocated.",err,error,*999)
#endif    
    
    equationsSetToSolverMatricesMap=>solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(equationsSetToSolverMatricesMap)) THEN
      localError="The equations set to solver matrices map is not associated for the equations set index of "// &
        & TRIM(NumberToVString(equationsSetIdx,"*",err,error))//" of the solver mapping."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
      
    EXITS("SolverMapping_EquationsSetToSolverMatricesMapGet")
    RETURN
999 NULLIFY(equationsSetToSolverMatricesMap)
998 ERRORSEXITS("SolverMapping_EquationsSetToSolverMatricesMapGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_EquationsSetToSolverMatricesMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the Lagrange variable to solver col map in an equations set mapped to a solver matrix in a solver mapping.
  SUBROUTINE SolverMapping_IMSToSMMapLagrangeVariableToSolverColMapGet(solverMapping,solverMatrixIndex,interfaceConditionIdx, &
    & lagrangeVariableToSolverColMap,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the Lagrange variable to solver col map for
    INTEGER(INTG), INTENT(IN) :: solverMatrixIndex !<The solver matrix index in the solver mapping to get the Lagrange variable to solver col map for
    INTEGER(INTG), INTENT(IN) :: interfaceConditionIdx !<The interface condition index containing the Lagrange variable mapped to the solver matrix
    TYPE(VariableToSolverColMapType), POINTER :: lagrangeVariableToSolverColMap !<On exit, a pointer to the specified Lagrange variable to solver col map for the solver matrix. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(InterfaceMatricesToSolverMatrixMapType), POINTER :: interfaceMatricesToSolverMatrixMap
    TYPE(InterfaceConditionToSolverMatricesMapType), POINTER :: interfaceConditionToSolverMatricesMap
#endif
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMapping_IMSToSMMapLagrangeVariableToSolverColMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(lagrangeVariableToSolverColMap)) &
      & CALL FlagError("Lagrange variable to solver col map is already associated.",err,error,*998)
    NULLIFY(interfaceConditionToSolverMatricesMap)
    CALL SolverMapping_InterfaceConditionToSolverMatricesMapGet(solverMapping,interfaceConditionIdx, &
      & interfaceConditionToSolverMatricesMap,err,error,*999)
    NULLIFY(interfaceMatricesToSolverMatrixMap)
    CALL SolverMappingICToSMSMap_InterfaceMatricesToSolverMatrixMapGet(interfaceConditionToSolverMatricesMap,solverMatrixIdx, &
      & interfaceMatricesToSolverMatrixMap,err,error,*999)
#endif    
     
   lagrangeVariableToSolverColMap=>solverMapping%InterfaceConditionToSolverMatricesMaps(equationsSetIdx)%ptr% &
     & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
     & lagrangeVariableToSolverColMap

#ifdef WITH_POSTCHECKS   
    IF(.NOT.ASSOCIATED(lagrangeVariableToSolverColMap)) THEN
      localError="The Lagrange variable to solver col map is not associated for the solver matrix index of "// &
        & TRIM(NumberToVString(solverMatrixIndex,"*",err,error))// &
        & " of the interface matrices to solver matrix maps for the interface condition index "// &
        & TRIM(NumberToVString(interfaceConditionIdx,"*",err,error))// &
        & " of the interface condition to solver matrices map of the solver mapping."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("SolverMapping_IMSToSMMapLagrangeVariableToSolverColMapGet")
    RETURN
999 NULLIFY(lagrangeVariableToSolverColMap)
998 ERRORS("SolverMapping_IMSToSMMapLagrangeVariableToSolverColMapGet",err,error)
    EXITS("SolverMapping_IMSToSMMapLagrangeVariableToSolverColMapGet")
    RETURN 1
    
  END SUBROUTINE SolverMapping_IMSToSMMapLagrangeVariableToSolverColMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the dependent variable to solver col map in an equations set mapped to a solver matrix in a solver mapping.
  SUBROUTINE SolverMapping_IMSToSMMapDependentVariableToSolverColMapGet(solverMapping,solverMatrixIndex,interfaceConditionIdx, &
    & variableIdx,dependentVariableToSolverColMap,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the dependent variable to solver col map for
    INTEGER(INTG), INTENT(IN) :: solverMatrixIndex !<The solver matrix index in the solver mapping to get the dependent variable to solver col map for
    INTEGER(INTG), INTENT(IN) :: interfaceConditionIdx !<The interface condition index containing the dependent variable mapped to the solver matrix
    INTEGER(INTG), INTENT(IN) :: variableIdx !<The variable index of the dependent variable mapped to the solver matrix
    TYPE(VariableToSolverColMapType), POINTER :: dependentVariableToSolverColMap !<On exit, a pointer to the specified dependent variable to solver col map for the solver matrix. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(InterfaceMatricesToSolverMatrixMapType), POINTER :: interfaceMatricesToSolverMatrixMap
    TYPE(InterfaceConditionToSolverMatricesMapType), POINTER :: interfaceConditionToSolverMatricesMap
#endif
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMapping_IMSToSMMapDependentVariableToSolverColMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dependentVariableToSolverColMap)) &
      & CALL FlagError("Dependent variable to solver col map is already associated.",err,error,*998)
    NULLIFY(interfaceConditionToSolverMatricesMap)
    CALL SolverMapping_InterfaceConditionToSolverMatricesMapGet(solverMapping,interfaceConditionIdx, &
      & interfaceConditionToSolverMatricesMap,err,error,*999)
    NULLIFY(interfaceMatricesToSolverMatrixMap)
    CALL SolverMappingICToSMSMap_InterfaceMatricesToSolverMatrixMapGet(interfaceConditionToSolverMatricesMap,solverMatrixIdx, &
      & interfaceMatricesToSolverMatrixMap,err,error,*999)
    IF(variableIdx<1.OR.variableIdx>interfaceMatricesToSolverMatrixMap%numberOfDependentVariables) THEN
      localError="The specified dependent variable index of "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
        & " is invalid for the solver matrix index of "//TRIM(NumberToVString(solverMatrixIndex,"*",err,error))// &
        & " of the interface matrices to solver matrix maps for the interface condition index "// &
        & TRIM(NumberToVString(interfaceConditionIdx,"*",err,error))//" of the interface condition to solver matrices map of "// &
        & "the solver mapping. The dependent variable index should be >= 1 and <= "// &
        & TRIM(NumberToVString(interfaceMatricesToSolverMatrixMap%numberOfDependentVariables,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(interfaceMatricesToSolverMatrixMap%dependentVariableToSolverColMaps)) THEN
      localError="The dependent variable to solver col maps array is not allocated for the solver matrix index of "// &
        & TRIM(NumberToVString(solverMatrixIndex,"*",err,error))// &
        & " of the interface matrices to solver matrix maps of the interface condition index "// &
        & TRIM(NumberToVString(interfaceConditionIdx,"*",err,error))// &
        & " of the interface condition to solver matrices map of the solver mapping."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
     
   dependentVariableToSolverColMap=>solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
     & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
     & dependentVariableToSolverColMaps(variableIdx)%ptr

#ifdef WITH_POSTCHECKS   
    IF(.NOT.ASSOCIATED(dependentVariableToSolverColMap)) THEN
      localError="The dependent variable to solver col map is not associated for the specified variable index of "// &
        & TRIM(NumberToVString(variableIdx,"*",err,error))//" for the solver matrix index of "// &
        & TRIM(NumberToVString(solverMatrixIndex,"*",err,error))// &
        & " of the interface matrices to solver matrix maps for the interface condition index "// &
        & TRIM(NumberToVString(interfaceConditionIdx,"*",err,error))// &
        & " of the interface condition to solver matrices map of the solver mapping."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("SolverMapping_IMSToSMMapDependentVariableToSolverColMapGet")
    RETURN
999 NULLIFY(dependentVariableToSolverColMap)
998 ERRORS("SolverMapping_IMSToSMMapDependentVariableToSolverColMapGet",err,error)
    EXITS("SolverMapping_IMSToSMMapDependentVariableToSolverColMapGet")
    RETURN 1
    
  END SUBROUTINE SolverMapping_IMSToSMMapDependentVariableToSolverColMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the interface col to solver cols map for solver mapping.
  SUBROUTINE SolverMapping_InterfaceColToSolverColsMapGet(solverMapping,interfaceConditionIdx,solverMatrixIndex, &
    & interfaceColToSolverColsMap,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the interface col to solver cols map for
    INTEGER(INTG), INTENT(IN) :: interfaceConditionIdx !<The interface condition index in the solver mapping to get the interface col to solver cols map for
    INTEGER(INTG), INTENT(IN) :: solverMatrixIndex !<The solver matrix index in the solver mapping to get the interface col to solver cols map for
    TYPE(MatrixRowColCouplingType), POINTER :: interfaceColToSolverColsMap(:) !<On exit, a pointer to the specified interface col to solver cols map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(InterfaceConditionToSolverMatricesMapType), POINTER :: interfaceConditionToSolverMatricesMap
    TYPE(InterfaceMatricesToSolverMatrixMapType), POINTER :: interfaceMatrixToSolverMatricesMap
#endif
#ifdef WITH_POSTCHECKS
    TYPE(VARYING_STRING) :: localError
#endif
 
    ENTERS("SolverMapping_InterfaceColToSolverColsMapGet",err,error,*998)

#ifdef WITH_PRECHECKS 
    IF(ASSOCIATED(interfaceColToSolverColsMap)) &
      & CALL FlagError("Interface col to solver cols map is already associated.",err,error,*998)
    NULLIFY(interfaceConditionToSolverMatricesMap)
    CALL SolverMapping_InterfaceConditionToSolverMatricesMapGet(solverMapping,interfaceConditionIdx, &
      & interfaceConditionsToSolverMatricesMap,err,error,*999)
    NULLIFY(interfaceMatricesToSolverMatrixMap)
    CALL SolverMappingICToSMSMap_InterfaceMatricesToSolverMatrixMapGet(interfaceMatricesToSolverMatrixMap,solverMatrixIdx, &
      & interfaceMatricesToSolverMatrixMap,err,error,*999)
#endif

    interfaceColToSolverColsMap=>solverMapping%interfaceConditionToSolverMatricesMap(interfaceConditionIdx)%ptr% &
      & interfaceMatricesToSolverMatrixMaps(solverMatrixIndex)%ptr%interfaceColToSolverColsMap

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(interfaceColToSolverColsMap)) THEN
      localError="The interface col to solver cols map for the specified solver matrix index of "// &
        & TRIM(NumberToVString(solverMatrixIdx,"*",err,error))//" and interface condition index of "// &
        & TRIM(NumberToVString(interfaceConditionIdx,"*",err,error))//" is not associated."      
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    EXITS("SolverMapping_InterfaceColToSolverColsMapGet")
    RETURN
999 NULLIFY(interfaceColToSolverColsMap)
998 ERRORS("SolverMapping_InterfaceColToSolverColsMapGet",err,error)
    EXITS("SolverMapping_InterfaceColToSolverColsMapGet")
    RETURN 1
    
  END SUBROUTINE SolverMapping_InterfaceColToSolverColsMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the interface col to solver rows map for solver mapping.
  SUBROUTINE SolverMapping_InterfaceColToSolverRowsMapGet(solverMapping,interfaceConditionIdx,interfaceColToSolverRowsMap, &
    & err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the interface row to solver rows map for
    INTEGER(INTG), INTENT(IN) :: interfaceConditionIdx !<The interface condition index in the solver mapping to get the interface col to solver rows map for
    TYPE(MatrixRowColCouplingType), POINTER :: interfaceColToSolverRowsMap(:) !<On exit, a pointer to the specified interface col to solver rows map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(InterfaceConditionToSolverMatricesMapType), POINTER :: interfaceConditionToSolverMatricesMap
#endif
#ifdef WITH_POSTCHECKS
    TYPE(VARYING_STRING) :: localError
#endif
 
    ENTERS("SolverMapping_InterfaceColToSolverRowsMapGet",err,error,*998)

#ifdef WITH_PRECHEKS    
    IF(ASSOCIATED(interfaceColToSolverRowsMap)) &
      & CALL FlagError("Interface col to solver rows map is already associated.",err,error,*998)
    NULLIFY(interfaceConditionToSolverMatricesMap)
    CALL SolverMapping_InterfaceConditionToSolverMatricesMapGet(solverMapping,interfaceConditionIdx, &
      & interfaceConditionToSolverMatricesMap,err,error,*999)
#endif

    interfaceColToSolverRowsMap=>solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
      & interfaceColToSolverRowsMap

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(interfaceColToSolverRowsMap)) THEN
      localError="The interface col to solver rows map for the specified interface condition index of "// &
        & TRIM(NumberToVString(interfaceConditionIdx,"*",err,error))//" is not associated."      
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("SolverMapping_InterfaceColToSolverRowsMapGet")
    RETURN
999 NULLIFY(interfaceColToSolverRowsMap)
998 ERRORS("SolverMapping_InterfaceColToSolverRowsMapGet",err,error)
    EXITS("SolverMapping_InterfaceColToSolverRowsMapGet")
    RETURN 1
    
  END SUBROUTINE SolverMapping_InterfaceColToSolverRowsMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the interface condition for solver mapping.
  SUBROUTINE SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the interface condition for
    INTEGER(INTG), INTENT(IN) :: interfaceConditionIdx !<The interface condition index in the solver mapping to get the interface condition for
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<On exit, a pointer to the specified interface condition. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif
 
    ENTERS("SolverMapping_InterfaceConditionGet",err,error,*998)

 #ifdef WITH_PRECHECKS    
   IF(ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMapping)) CALL FlagError("Solver mapping is not associated.",err,error,*999)
    IF(interfaceConditionIdx<1.OR.interfaceConditionIdx>solverMapping%numberOfInterfaceConditions) THEN
      localError="The specified interface condition index of "//TRIM(NumberToVString(interfaceConditionIdx,"*",err,error))// &
        & " is invalid. The index must be >= 1 and <= "// &
        & TRIM(NumberToVString(SolverMapping%numberOfInterfaceConditions,"*",err,error))//"."      
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(solverMapping%interfaceConditions)) &
      & CALL FlagError("Solver mapping interface conditions is not allocated.",err,error,*999)
#endif

    interfaceCondition=>solverMapping%interfaceConditions(interfaceConditionIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(interfaceCondition)) THEN
      localError="The interface condition for the specified interface condition index of "// &
        & TRIM(NumberToVString(interfaceConditionIdx,"*",err,error))//" is not associated in the solver mapping."      
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("SolverMapping_InterfaceConditionGet")
    RETURN
999 NULLIFY(interfaceCondition)
998 ERRORSEXITS("SolverMapping_InterfaceConditionGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_InterfaceConditionGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the interface condition to solver matrices map for a solver mapping.
  SUBROUTINE SolverMapping_InterfaceConditionToSolverMatricesMapGet(solverMapping,interfaceConditionIdx, &
    & interfaceConditionToSolverMatricesMap,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the equations set to solver matrices map for
    INTEGER(INTG), INTENT(IN) :: interfaceConditionIdx !<The interface condition index in the solver mapping to get the interface condition to solver matrices map for
    TYPE(InterfaceConditionToSolverMatricesMapType), POINTER :: interfaceConditionToSolverMatricesMap !<On exit, a pointer to the specified interface condition to solver matrices map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMapping_InterfaceConditionToSolverMatricesMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(interfaceConditionToSolverMatricesMap)) &
      & CALL FlagError("Interface conditionto solver matrices map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMapping)) CALL FlagError("Solver mapping is not associated.",err,error,*999)
    IF(interfaceConditionIdx<1.OR.interfaceConditionIdx>solverMapping%numberOfInterfaceConditions) THEN
      localError="The specified interface condition index of "//TRIM(NumberToVString(interfaceConditionIdx,"*",err,error))// &
        & " is invalid for the solver mapping. The index must be >= 1 and <= "// &
        & TRIM(NumberToVString(solverMapping%numberOfInterfaceConditions,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(solverMapping%interfaceConditionToSolverMatricesMaps)) &
      & CALL FlagError("Solver mapping interface condition to solver matrices maps is not allocated.",err,error,*999)
#endif    
    
    interfaceConditionToSolverMatricesMap=>solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(interfaceConditionToSolverMatricesMap)) THEN
      localError="The interface condition to solver matrices map is not associated for the interface condition index of "// &
        & TRIM(NumberToVString(interfaceConditionIdx,"*",err,error))//" of the solver mapping."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
      
    EXITS("SolverMapping_InterfaceConditionToSolverMatricesMapGet")
    RETURN
999 NULLIFY(interfaceConditionToSolverMatricesMap)
998 ERRORSEXITS("SolverMapping_InterfaceConditionToSolverMatricesMapGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_InterfaceConditionToSolverMatricesMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a interface equations matrix in an interface condition mapped to a solver matrix in a solver mapping.
  SUBROUTINE SolverMapping_InterfaceEquationsMatrixToSolverMatrixGet(solverMapping,solverMatrixIndex,interfaceConditionIdx, &
    & interfaceMatrixIdx,interfaceEquationsMatrix,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the interface equations matrix for
    INTEGER(INTG), INTENT(IN) :: solverMatrixIndex !<The solver matrix index in the solver mapping to get the interface equations matrix for
    INTEGER(INTG), INTENT(IN) :: interfaceConditionIdx !<The interface condition index containing the interface matrices mapped to the solver matrix
    INTEGER(INTG), INTENT(IN) :: interfaceMatrixIdx !<The interface matrix index of the interface equations matrix to the solver matrix
    TYPE(EquationsMatrixIdx), POINTER :: interfaceEquationsMatrix !<On return the specified interface equations matrix. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(InterfaceConditionToSolverMatricesMapType), POINTER :: interfaceConditionToSolverMatricesMap
    TYPE(InterfaceMatricesToSolverMatrixMapType), POINTER :: interfaceMatricesToSolverMatrixMaps
    TYPE(InterfaceMatrixToSolverMatrixMapType), POINTER :: interfaceMatrixToSolverMatrixMap
#endif
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMapping_InterfaceEquationsMatrixToSolverMatrixGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(interfaceEquationsMatrix)) CALL FlagError("Interface equations matrix is already associated.",err,eror,*998)
    NULLIFY(interfaceConditionToSolverMatricesMap)
    CALL SolverMapping_InterfaceConditionToSolverMatricesMapGet(solverMapping,interfaceConditionIdx, &
      & interfaceConditionToSolverMatricesMap,err,error,*999)
    NULLIFY(interfaceMatricesToSolverMatrixMaps)
    CALL SolverMappingICToSMSMap_InterfaceMatricesToSolverMatrixMapGet(interfaceConditionToSolverMatricesMap,solverMatrixIdx, &
      & interfaceMatricesToSolverMatrixMaps,err,error,*999)
    NULLIFY(interfaceMatrixToSolverMatrixMap)
    CALL SolverMappingIMSToSMMap_InterfaceMatrixToSolverMatrixMapGet(interfaceMatricesToSolverMatrixMaps,interfaceMatrixIdx, &
      & interfaceMatrixToSolverMatrixMap,err,error,*999)
#endif    

    interfaceEquationsMatrix=>solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
      & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
      & interfaceMatrixToSolverMatrixMaps(interfaceMatrixIdx)%ptr%interfaceMatrix

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(interfaceEquationsMatrix)) THEN
      localError="The interface equations matrix is not associated for the interface matrix index of "// &
        & TRIM(NumberToVString(interfaceMatrixIdx,"*",err,error))// &
        & " of the interface matrix to solver matrix maps for the solver matrix index of "// &
        & TRIM(NumberToVString(solverMatrixIndex,"*",err,error))// &
        & " of the interface matrices to solver matrix maps for the interface condition index "// &
        & TRIM(NumberToVString(interfaceConditionIdx,"*",err,error))// &
        & " of the interface condition to solver matrices map of the solver mapping."
      CALL FlagError(localError,err,error,*999)    
    ENDIF
#endif    
      
    EXITS("SolverMapping_InterfaceEquationsMatrixToSolverMatrixGet")
    RETURN
999 NULLIFY(interfaceEquationsMatrix)
998 ERRORSEXITS("SolverMapping_InterfaceEquationsMatrixToSolverMatrixGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_InterfaceEquationsMatrixToSolverMatrixGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the interface matrix to solver matrix map for a solver mapping.
  SUBROUTINE SolverMapping_InterfaceMatrixToSolverMatrixMapGet(solverMapping,interfaceConditionIdx,interfaceMatrixIdx, &
    & solverMatrixIdx,interfaceMatrixToSolverMatrixMap,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the interface matrix to solver matrix map for
    INTEGER(INTG), INTENT(IN) :: interfaceConditionIdx !<The interface condition index in the solver mapping to get the interface matrix to solver matrix map for
    INTEGER(INTG), INTENT(IN) :: interfaceMatrixIdx !<The interface matrix index for the interface matrix to solver matrix map
    INTEGER(INTG), INTENT(IN) :: solverMatrixIdx !<The solver matrix index for the interface matrix to solver matrix map
    TYPE(InterfaceToSolverMapsType), POINTER :: interfaceMatrixToSolverMatrixMap !<On exit, a pointer to the specified interface matrix to solver matrix map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(InterfaceConditionToSolverMatricesMapType), POINTER :: interfaceConditionToSolverMatricesMap
    TYPE(InterfaceMatrixToSolverMatricesMapType), POINTER :: interfaceMatrixToSolverMatricesMap
#endif
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif
 
    ENTERS("SolverMapping_InterfaceMatrixToSolverMatrixMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(interfaceMatrixToSolverMatrixMap)) &
      & CALL FlagError("Interface matrix to solver matrix map is already associated.",err,error,*998)
    NULLIFY(interfaceConditionToSolverMatricesMap)
    CALL SolverMapping_InterfaceConditionToSolverMatricesMapGet(solverMapping,interfaceConditionIdx, &
      & interfaceConditionToSolverMatricesMap,err,error,*999)
    NULLIFY(interfaceMatrixToSolverMatricesMap)
    CALL SolverMappingICToSMSMap_InterfaceMatrixToSolverMatricesMapGet(interfaceMatrixToSolverMatricesMap,interfaceMatrixIdx, &
      & interfaceMatrixToSolverMatricesMap,err,error,*999)
    IF(solverMatrixIdx<1.OR.solverMatrixIdx>interfaceMatrixToSolverMatriceMap%numberOfSolverMatrices) THEN
      localError="The specified solver matrix index of "//TRIM(NumberToVString(solverMatrixIdx,"*",err,error))//&
        & " is invalid for the interface matrix index of "//TRIM(NumberToVString(interfaceMatrixIdx,"*",err,error))//&
        & " of the interface matrix to solver matrices maps for the interface condition index of "// &
        & TRIM(NumberToVString(interfaceConditionIdx,"*",err,error))//" of the interface condition to solver matrices maps "// &
        & "of the solver mapping. The solver matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(interfaceMatrixToSolverMatriceMap%numberOfSolverMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)      
    ENDIF
    IF(.NOT.ALLOCATED(interfaceMatrixToSolverMatriceMap%interfaceMatrixToSolverMatrixMaps)) THEN
      localError="The interface matrix to solver matrix maps is not allocated for the interface matrix index of "// &
        & TRIM(NumberToVString(interfaceMatrixIdx,"*",err,error))//&
        & " of the interface matrix to solver matrices maps for the interface condition index of "// &
        & TRIM(NumberToVString(interfaceConditionIdx,"*",err,error))//" of the interface condition to solver matrices maps "// &
        & "of the solver mapping."
    ENDIF
#endif    
    
    interfaceMatrixToSolverMatrixMap=>solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
      & interfaceMatrixToSolverMatricesMaps(interfaceMatrixIdx)%ptr% &
      & interfaceMatrixToSolverMatrixMaps(solverMatrixIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(interfaceMatrixToSolverMatrixMap)) THEN
      localError="The interface matrix to solver matrix map is not associated for the solver matrix index of "// &
        & TRIM(NumberToVString(solverMatrixIdx,"*",err,error))// &
        & " of the interface matrix to solver matrix maps of the interface matrix index of "// &
        & TRIM(NumberToVString(interfaceMatrixIdx,"*",err,error))// &
        & " of the interface matrix to solver matrices map of the interface condition index of "// &
        & TRIM(NumberToVString(interfaceConditionIdx,"*",err,error))// &
        & " of the interface condition to solver matrices map of the solver mapping."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("SolverMapping_InterfaceMatrixToSolverMatrixMapGet")
    RETURN
999 NULLIFY(interfaceMatrixToSolverMatrixMap)
998 ERRORS("SolverMapping_InterfaceMatrixToSolverMatrixMapGet",err,error)
    EXITS("SolverMapping_InterfaceMatrixToSolverMatrixMapGet")
    RETURN 1
    
  END SUBROUTINE SolverMapping_InterfaceMatrixToSolverMatrixMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the interface row to solver rows map for solver mapping.
  SUBROUTINE SolverMapping_InterfaceRowToSolverRowsMapGet(solverMapping,interfaceConditionIdx,interfaceMatrixIndex, &
    & interfaceRowToSolverRowsMap,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the interface row to solver rows map for
    INTEGER(INTG), INTENT(IN) :: interfaceConditionIdx !<The interface condition index in the solver mapping to get the interface row to solver rows map for
    INTEGER(INTG), INTENT(IN) :: interfaceMatrixIndex !<The solver matrix index in the solver mapping to get the interface row to solver rows map for
    TYPE(MatrixRowColCouplingType), POINTER :: interfaceRowToSolverRowsMap(:) !<On exit, a pointer to the specified interface row to solver rows map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(InterfaceConditionToSolverMatricesMapType), POINTER :: interfaceConditionToSolverMatricesMap
    TYPE(InterfaceMatricesToSolverMatrixMapType), POINTER :: interfaceMatrixToSolverMatricesMap
#endif
#ifdef WITH_POSTCHECKS
    TYPE(VARYING_STRING) :: localError
#endif
 
    ENTERS("SolverMapping_InterfaceRowToSolverRowsMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(interfaceRowToSolverRowsMap)) &
      & CALL FlagError("Interface row to solver rows map is already associated.",err,error,*998)
    NULLIFY(interfaceConditionToSolverMatricesMap)
    CALL SolverMapping_InterfaceConditionToSolverMatricesMapGet(solverMapping,interfaceConditionIdx, &
      & interfaceConditionToSolverMatricesMap,err,error,*999)
    NULLIFY(interfaceMatrixToSolverMatricesMap)
    CALL SolverMappingICToSMSMap_InterfaceMatrixToSolverMatricesMapGet(interfaceConditionToSolverMatricesMap, &
      & interfaceMatrixIdx,interfaceMatrixToSolverMatricesMap,err,error,*999)
#endif    

    interfaceRowToSolverRowsMap=>solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
      & interfaceMatrixToSolverMatricesMaps(interfaceMatrixIndex)%ptr%interfaceRowToSolverRowsMap

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(interfaceRowToSolverRowsMap)) THEN
      localError="The interface row to solver rows map for the specified interface matrix index of "// &
        & TRIM(NumberToVString(interfaceMatrixIdx,"*",err,error))//" of interface condition index "// &
        & TRIM(NumberToVString(interfaceConditionIdx,"*",err,error))//" is not associated."      
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("SolverMapping_InterfaceRowToSolverRowsMapGet")
    RETURN
999 NULLIFY(interfaceRowToSolverRowsMap)
998 ERRORS("SolverMapping_InterfaceRowToSolverRowsMapGet",err,error)
    EXITS("SolverMapping_InterfaceRowToSolverRowsMapGet")
    RETURN 1
    
  END SUBROUTINE SolverMapping_InterfaceRowToSolverRowsMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a Jacobian equations matrix in an equations set mapped to a solver matrix in a solver mapping.
  SUBROUTINE SolverMapping_JacobianEquationsMatrixToSolverMatrixGet(solverMapping,solverMatrixIndex,equationsSetIdx, &
    & jacobianMatrixIdx,jacobianMatrix,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the Jacobian equations matrix for
    INTEGER(INTG), INTENT(IN) :: solverMatrixIndex !<The solver matrix index in the solver mapping to get the Jacobian equations matrix for
    INTEGER(INTG), INTENT(IN) :: equationsSetIdx !<The equations set index containing the Jacobian matrices mapped to the solver matrix
    INTEGER(INTG), INTENT(IN) :: jacobianMatrixIdx !<The Jacobian matrix index of the Jacobian matrix to the solver matrix
    TYPE(EquationsMatrixIdx), POINTER :: jacobianEquationsMatrix !<On return the specified Jacobian equations matrix. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(EquationsMatricesToSolverMatrixMapType), POINTER :: equationsMatricesToSolverMatrixMap
    TYPE(EquationsSetToSolverMatricesMapType), POINTER :: equationsSetToSolverMatricesMap
    TYPE(JacobianMatrixToSolverMatrixMapType), POINTER :: jacobianMatrixToSolverMatrixMap
#endif
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMapping_JacobianEquationsMatrixToSolverMatrixGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(jacobianMatrix)) CALL FlagError("Jacobian equations matrix is already associated.",err,eror,*998)
    NULLIFY(equationsSetToSolverMatricesMap)
    CALL SolverMapping_EquationsSetToSolverMatricesMapGet(solverMapping,equationsSetIdx,equationsSetToSolverMatricesMap, &
      & err,error,*999)
    NULLIFY(equationsMatricesToSolverMatrixMap)
    CALL SolverMappingESToSMSMap_EquationsMatricesToSolverMatrixMapGet(equationsSetToSolverMatricesMap,solverMatrixIdx, &
      & equationsMatricesToSolverMatrixMap,err,error,*999)
    NULLIFY(jacobianMatrixToSolverMatrixMap)
    CALL SolverMappinEMSToSMMap_JacobianMatrixToSolverMatrixMapGet(equationsMatricesToSolverMatrixMap,jacobianMatrixIdx, &
      & jacobianMatrixToSolverMatrixMap,err,error,*999)
#endif
    
    jacobianEquationsMatrix=>solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr &
      & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
      & jacobianMatrixToSolverMatrixMaps(jacobianMatrixIdx)%ptr%jacobianMatrix

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(jacobianEquationsMatrix)) THEN
      localError="The Jacobian equations matrix is not associated for the Jacobian matrix index of "// &
        & TRIM(NumberToVString(jacobianMatrixIdx,"*",err,error))// &
        & " of the Jacobian matrix to solver matrix maps for the solver matrix index of "// &
        & TRIM(NumberToVString(solverMatrixIndex,"*",err,error))// &
        & " of the equations matrices to solver matrix maps for the equations set index "// &
        & TRIM(NumberToVString(equationsSetIdx,"*",err,error))// &
        & " of the equations set to solver matrices maps of the solver mapping."
      CALL FlagError(localError,err,error,*999)    
    ENDIF
#endif
    
    EXITS("SolverMapping_JacobianEquationsMatrixToSolverMatrixGet")
    RETURN
999 NULLIFY(jacobianEquationsMatrix)
998 ERRORSEXITS("SolverMapping_JacobianEquationsMatrixToSolverMatrixGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_JacobianEquationsMatrixToSolverMatrixGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the Jacobian matrix to solver matrix map for a solver mapping.
  SUBROUTINE SolverMapping_JacobianMatrixToSolverMatrixMapGet(solverMapping,equationsSetIndex,jacobianMatrixIndex, &
    & jacobianMatrixToSolverMatrixMap,err,error,*)
    
    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the Jacobian matrix to solver matrix map for
    INTEGER(INTG), INTENT(IN) :: equationsSetIndex !<The equations set index in the solver mapping to get the Jacobian matrix to solver matrix map for
    INTEGER(INTG), INTENT(IN) :: jacobianMatrixIndex !<The jacobian matrix index for the Jacobian matrix to solver matrix map
    TYPE(JacobianMatrixToSolverMatrixMapType), POINTER :: jacobianMatrixToSolverMatrixMap !<On exit, a pointer to the specified Jacobian matrix to solver matrix map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(EquationsSetToSolverMatricesMapType), POINTER :: equationsSetToSolverMatricesMap
#endif
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif
 
    ENTERS("SolverMapping_JacobianMatrixToSolverMatrixMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(jacobianMatrixToSolverMatrixMap)) &
      & CALL FlagError("Jacobian matrix to solver matrix map is already associated.",err,error,*998)
    NULLIFY(equationsSetToSolverMatricesMap)
    CALL SolverMapping_EquationsSetToSolverMatricesMapGet(solverMapping,equationsSetIdx,equationsSetToSolverMatricesMap, &
      & err,error,*999)    
    IF(jacobianMatrixIndex<1.OR.jacobianMatrixIndex>SIZE(equationsSetToSolverMatricesMap%jacobianMatrixToSolverMatrixMaps,1)) THEN
      localError="The specified Jacobian matrix index of "//TRIM(NumberToVString(jacobianMatrixIndex,"*",err,error))// &
        & " is invalid for the equations set index of "//TRIM(NumberToVString(equationsSetIndex,"*",err,error))// &
        & " of the equations set to solver matrices map of the solver mapping. The Jacobian matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(equationsSetToSolverMatricesMap%jacobianMatrixToSolverMatrixMaps,1),"*",err,error))//"."
      CALL FlagError(localError,err,eror,*999)
    ENDIF
    IF(.NOT.ALLOCATED(equationsSetToSolverMatricesMap%jacobianMatrixToSolverMatrixMaps)) THEN
      localError="The Jacobian matrix to solver matrix maps is not allocated for the equations set index of "// &
        & TRIM(NumberToVString(equationsSetIndex,"*",err,error))// &
        & " of the equations set to solver matrices map of the solver mapping."
      CALL FlagError(localError,err,eror,*999)
    ENDIF
#endif
   
    jacobianMatrixToSolverMatrixMap=>solverMapping%equationsSetToSolverMatricesMaps(equationsSetIndex)%ptr% &
      & jacobianMatrixToSolverMatrixMaps(jacobianMatrixIndex)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(jacobianMatrixToSolverMatrixMap)) THEN
      localError="The Jacobian matrix to solver matrix map is not associated for the Jacobian matrix index of "// &
        & TRIM(NumberToVString(jacobianMatrixIndex,"*",err,error))// &
        & " of the Jacobian matrix to solver matrix maps of the equations set index of "// &
        & TRIM(NumberToVString(equationsSetIndex,"*",err,error))// &
        & " of the equation set to solver matrices maps of the solver mapping."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("SolverMapping_JacobianMatrixToSolverMatrixMapGet")
    RETURN
999 NULLIFY(jacobianMatrixToSolverMatrixMap)
998 ERRORS("SolverMapping_JacobianMatrixToSolverMatrixMapGet",err,error)
    EXITS("SolverMapping_JacobianMatrixToSolverMatrixMapGet")
    RETURN 1
    
  END SUBROUTINE SolverMapping_JacobianMatrixToSolverMatrixMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a linear equations matrix in an equations set mapped to a solver matrix in a solver mapping.
  SUBROUTINE SolverMapping_LinearEquationsMatrixToSolverMatrixGet(solverMapping,solverMatrixIndex,equationsSetIdx, &
    & linearMatrixIdx,linearMatrix,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the linear equations matrix for
    INTEGER(INTG), INTENT(IN) :: solverMatrixIndex !<The solver matrix index in the solver mapping to get the linear equations matrix for
    INTEGER(INTG), INTENT(IN) :: equationsSetIdx !<The equations set index containing the linear matrices mapped to the solver matrix
    INTEGER(INTG), INTENT(IN) :: linearMatrixIdx !<The linear matrix index of the linear equations matrix to the solver matrix
    TYPE(EquationsMatrixIdx), POINTER :: linearEquationsMatrix !<On return the specified linear equations matrix. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(LinearEquationsMatrixToSolverMatrixMapType), POINTER :: linearEquationsMatrixToSolverMatrixMap
    TYPE(EquationsMatricesToSolverMatrixMapType), POINTER :: equationsMatricesToSolverMatrixMap
    TYPE(EquationsSetToSolverMatricesMapType), POINTER :: equationsSetToSolverMatricesMap
#endif
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMapping_LinearEquationsMatrixToSolverMatrixGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(linearMatrix)) CALL FlagError("Linear equations matrix is already associated.",err,eror,*998)
    NULLIFY(equationsSetToSolverMatricesMap)
    CALL SolverMapping_EquationsSetToSolverMatricesMapGet(solverMapping,equationsSetIdx,equationsSetToSolverMatricesMap, &
      & err,error,*999)
    NULLIFY(equationsMatricesToSolverMatrixMap)
    CALL SolverMappingESToSMSMap_EquationsMatricesToSolverMatrixMapGet(equationsSetToSolverMatricesMap,solverMatrixIdx, &
      & equationsMatricesToSolverMatrixMap,err,error,*999)
    NULLIFY(linearEquationsMatrixToSolverMatrixMap)
    CALL SolverMappingEMSToSMMap_LinearEquationsMatrixToSolverMatrixMapGet(equationsMatricesToSolverMatrixMap, &
      & linearMatrixIdx,linearEquationsMatrixToSolverMatrixMap,err,error,*999)
#endif    

    linearEquationsMatrix=>solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
      & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
      & linearEquationsMatrixToSolverMatrixMaps(linearMatrixIdx)%ptr%equationsMatrix

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(linearEquationsMatrix)) THEN
      localError="The linear equations matrix is not associated for the linear matrix index of "// &
        & TRIM(NumberToVString(linearMatrixIdx,"*",err,error))// &
        & " of the linear equations matrix to solver matrix maps for the solver matrix index of "// &
        & TRIM(NumberToVString(solverMatrixIndex,"*",err,error))// &
        & " of the equations matrices to solver matrix maps for the equations set index "// &
        & TRIM(NumberToVString(equationsSetIdx,"*",err,error))// &
        & " of the equations set to solver matrices map of the solver mapping."
      CALL FlagError(localError,err,error,*999)    
    ENDIF
#endif    
      
    EXITS("SolverMapping_LinearEquationsMatrixToSolverMatrixGet")
    RETURN
999 NULLIFY(linearEquationsMatrix)
998 ERRORSEXITS("SolverMapping_LinearEquationsMatrixToSolverMatrixGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_LinearEquationsMatrixToSolverMatrixGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the row DOFs mapping for solver mapping.
  SUBROUTINE SolverMapping_RowDOFSMappingGet(solverMapping,rowDOFSMapping,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the row dofs mapping for
    TYPE(DomainMappingType), POINTER :: rowDOFSMappng !<On exit, a pointer to the specified row DOFS mapping. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMapping_RowDOFSMappingGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(rowDOFSMapping)) CALL FlagError("Row DOFS mapping is already associated.",err,error,*998)    
    IF(.NOT.ASSOCIATED(solverMapping)) CALL FlagError("Solver mapping is not associated.",err,error,*999)
#endif

    rowDOFSMapping=>solverMapping%rowDOFsMapping

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(rowDOFSMapping)) &
      & CALL FlagError("The row DOFs mapping is not associated for the solver mapping.",err,error,*999)
#endif    
      
    EXITS("SolverMapping_RowDOFSMappingGet")
    RETURN
999 NULLIFY(rowDOFSMapping)
998 ERRORSEXITS("SolverMapping_RowDOFSMappingGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_RowDOFSMappingGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the solver matrix to equations map for solver mapping.
  SUBROUTINE SolverMapping_SolverMatrixToEquationsMapGet(solverMapping,solverMatrixIdx,solverMatrixToEquationsMap,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the solver matrix to equations map for
    INTEGER(INTG), INTENT(IN) :: solverMatrixIdx !<The solver matrix index to get the solver matrix to equations map for
    TYPE(SolverMatrixToEquationsMapType), POINTER :: solverMatrixToEquationsMap !<On exit, a pointer to the specified solver matrix to equations map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif
 
    ENTERS("SolverMapping_SolverMatrixToEquationsMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(solverMatrixToEquationsMap)) &
      & CALL FlagError("Solver matrix to equations map is already associated.",err,error,*998)  
    IF(.NOT.ASSOCIATED(solverMapping)) CALL FlagError("Solver mapping is not associated.",err,error,*999)
    IF(solverMatrixIdx<1.OR.solverMatrixIdx>solverMapping%numberOfSolverMatrices) THEN
      localError="The specified solver matrix index of "//TRIM(NumberToVString(solverMatrixIdx,"*",err,error))// &
        & " is invalid. The solver matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(solverMapping%numberOfSolverMatrices,"*",err,error))//"."
    ENDIF
    IF(.NOT.ALLOCATED(solverMapping%solverMatrixToEquationsMaps)) &
      & CALL FlagError("The solver matrix to equations maps is not allocated for the solver mapping.",err,error,*999)
#endif

    solverMatrixToEquationsMap=>solverMapping%solverMatrixToEquationsMaps(solverMatrix)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(solverMatrixToEquationsMap)) THEN
      localError="The solver matrix to equations map is not associated for solver matrix index "// &
        & TRIM(NumberToVString(solverMatrixIdx,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("SolverMapping_SolverMatrixToEquationsMapGet")
    RETURN
999 NULLIFY(solverMatrixToEquationsMap)
998 ERRORSEXITS("SolverMapping_SolverMatrixToEquationsMapGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_SolverMatrixToEquationsMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to a equations variable list for solver mapping create values cache.
  SUBROUTINE SolverMappingCreateValuesCache_EquationsVariableListGet(createValuesCache,solverMatrixIdx,equationsVariableList, &
    & err,error,*)

    !Argument variables
    TYPE(SolverMappingCreateValuesCacheType), POINTER :: createValuesCache !<A pointer to the solver mapping create values cache to get the equations variable list for
    INTEGER(INTG), INTENT(IN) :: solverMatrixIdx !<The solver matrix index to get the equations variable list for
    TYPE(ListType), POINTER :: equationsVariableList !<On exit, a pointer to the specified equations variable list. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif
 
    ENTERS("SolverMappingCreateValuesCache_EquationsVariableListGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(equationsVariableList)) &
      & CALL FlagError("Equations variable list is already associated.",err,error,*998)  
    IF(.NOT.ASSOCIATED(createValuesCache)) CALL FlagError("Solver mapping create values cache is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(createValuesCache%equationsVariableList)) &
      & CALL FlagError("The equations variable list is not associated for the solver mapping create values cache.",err,error,*999)
    IF(solverMatrixIdx<1.OR.solverMatrixIdx>SIZE(createValuesCache%equationsVariableList,1)) THEN
      localError="The specified solver matrix index of "//TRIM(NumberToVString(solverMatrixIdx,"*",err,error))// &
        & " is invalid. The solver matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(createValuesCache%equationsVariableList,1),"*",err,error))//"."
    ENDIF
#endif

    equationsVariableList=>createValuesCache%equationsVariableList(solverMatrix)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(equationsVariableList)) THEN
      localError="The equations variable list is not associated for solver matrix index "// &
        & TRIM(NumberToVString(solverMatrixIdx,"*",err,error))//" for the solver mapping create values cache."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("SolverMappingCreateValuesCache_EquationsVariableListGet")
    RETURN
999 NULLIFY(equationsVariableList)
998 ERRORS("SolverMappingCreateValuesCache_EquationsVariableListGet",err,error)
    EXITS("SolverMappingCreateValuesCache_EquationsVariableListGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingCreateValuesCache_EquationsVariableListGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to a interface variable list for solver mapping create values cache.
  SUBROUTINE SolverMappingCreateValuesCache_InterfaceVariableListGet(createValuesCache,solverMatrixIdx,interfaceVariableList, &
    & err,error,*)

    !Argument variables
    TYPE(SolverMappingCreateValuesCacheType), POINTER :: createValuesCache !<A pointer to the solver mapping create values cache to get the interface variable list for
    INTEGER(INTG), INTENT(IN) :: solverMatrixIdx !<The solver matrix index to get the interface variable list for
    TYPE(ListType), POINTER :: interfaceVariableList !<On exit, a pointer to the specified interface variable list. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif
 
    ENTERS("SolverMappingCreateValuesCache_InterfaceVariableListGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(interfaceVariableList)) &
      & CALL FlagError("Interface variable list is already associated.",err,error,*998)  
    IF(.NOT.ASSOCIATED(createValuesCache)) CALL FlagError("Solver mapping create values cache is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(createValuesCache%interfaceVariableList)) &
      & CALL FlagError("The interface variable list is not associated for the solver mapping create values cache.",err,error,*999)
    IF(solverMatrixIdx<1.OR.solverMatrixIdx>SIZE(createValuesCache%interfaceVariableList,1)) THEN
      localError="The specified solver matrix index of "//TRIM(NumberToVString(solverMatrixIdx,"*",err,error))// &
        & " is invalid. The solver matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(createValuesCache%interfacesVariableList,1),"*",err,error))//"."
    ENDIF
#endif

    interfaceVariableList=>createValuesCache%interfaceVariableList(solverMatrix)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(interfaceVariableList)) THEN
      localError="The interface variable list is not associated for solver matrix index "// &
        & TRIM(NumberToVString(solverMatrixIdx,"*",err,error))//" for the solver mapping create values cache."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("SolverMappingCreateValuesCache_InterfaceVariableListGet")
    RETURN
999 NULLIFY(interfaceVariableList)
998 ERRORS("SolverMappingCreateValuesCache_InterfaceVariableListGet",err,error)
    EXITS("SolverMappingCreateValuesCache_InterfaceVariableListGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingCreateValuesCache_InterfaceVariableListGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to a RHS variable list for solver mapping create values cache.
  SUBROUTINE SolverMappingCreateValuesCache_RHSVariableListGet(createValuesCache,rhsVariableList,err,error,*)

    !Argument variables
    TYPE(SolverMappingCreateValuesCacheType), POINTER :: createValuesCache !<A pointer to the solver mapping create values cache to get the RHS variable list for
    TYPE(ListType), POINTER :: rhsVariableList !<On exit, a pointer to the specified RHS variable list. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingCreateValuesCache_RHSVariableListGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(rhsVariableList)) &
      & CALL FlagError("RHS variable list is already associated.",err,error,*998)  
    IF(.NOT.ASSOCIATED(createValuesCache)) CALL FlagError("Solver mapping create values cache is not associated.",err,error,*999)
#endif

    rhsVariableList=>createValuesCache%equaitonsRHSVariableList

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(rhsVariableList)) &
      & CALL FlagError("The RHS variable list is not associated for the solver mapping create values cache.",err,error,*999)
#endif    
      
    EXITS("SolverMappingCreateValuesCache_RHSVariableListGet")
    RETURN
999 NULLIFY(rhsVariableList)
998 ERRORS("SolverMappingCreateValuesCache_RHSVariableListGet",err,error)
    EXITS("SolverMappingCreateValuesCache_RHSVariableListGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingCreateValuesCache_EquationsVariableListGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the dynamic equations matrix to solver matrix map for an equations matrices to solver matrix map.
  SUBROUTINE SolverMappingEMSToSMMap_DynamicEquationsMatrixToSolverMatrixMapGet(equationsMatricesToSolverMatrixMap, &
    & dynamicEquationsMatrixIdx,dynamicEquationsMatrixToSolverMatrixMap,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesToSolverMatrixMapType), POINTER :: equationsMatricesToSolverMatrixMap !<A pointer to the equations matrices to solver matrix map to get the dynamic equations matrix to solver matrix map for
    INTEGER(INTG), INTENT(IN) :: dynamicEquationsMatrixIdx !<The dynamic equations matrix index to get the dynamic equations matrix to solver matrix for
    TYPE(EquationsMatrixToSolverMatrixMapType), POINTER :: dynamicEquationsMatrixToSolverMatrixMap !<On exit, a pointer to the specified dynamic equations matrix to solver matrix map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMappingEMSToSMMap_DynamicEquationsMatrixToSolverMatrixGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dynamicEquationsMatrixToSolverMatrixMap)) &
      & CALL FlagError("Dynamic equations matrix to solver matrix map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsMatricesToSolverMatrixMap)) &
      & CALL FlagError("Equations matrices to solver matrix map is not associated.",err,error,*999)
    IF(dynamicEquationsMatrixIdx<1.OR. &
      & dynamicEquationsMatrixIdx>equationsMatricesToSolverMatrixMap%numberOfDynamicEquationsMatrices) THEN
      localError="The specified dynamic equations matrix index of "// &
        & TRIM(NumberToVString(dynamicEquationMatrixIdx,"*",err,error))//" is invalid for the equations matrices to solver "// &
        & "matrix map. The dynamic equations matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(equationsMatricesToSolverMatrixMap%numberOfDynamicEquationsMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.MOT.ALLOCATED(equationsMatricesToSolverMatrixMap%dynamicEquationsMatrixToSolverMatrixMaps)) &
      & CALL FlagError("The dynamic equations to solver matrix maps array is not allocated for the equations matrices to "// &
      & "solver matrix map.",err,error,*999)
#endif    

    dynamicEquationsMatrixToSolverMap=>equationsMatricesToSolverMatrixMap% &
      & dynamicEquationsMatrixToSolverMatrixMaps(dynamicEquationsMatrixIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(dynamicEquationsMatrix)) THEN
      localError="The dynamic equations matrix to solver map is not associated for dynamic equations matrix index "// &
        & TRIM(NumberToVString(dynamicEquationsMatrix,"*",err,error))// &
        & " of the equations matrices to solver matrix map."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("SolverMappingEMSToSMMap_DynamicEquationsMatrixToSolverMatrixMapGet")
    RETURN
999 NULLIFY(dynamicEquationsMatrix)
998 ERRORS("SolverMappingEMSToSMMap_DynamicEquationsMatrixToSolverMatrixMapGet",err,error)
    EXITS("SolverMappingEMSToSMMap_DynamicEquationsMatrixToSolverMatrixMapGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingEMSToSMMap_DynamicEquationsMatrixToSolverMatrixMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the Jacobian equations matrix to solver matrix map for an equations matrices to solver matrix map.
  SUBROUTINE SolverMappingEMSToSMMap_JacobianEquationsMatrixToSolverMatrixMapGet(equationsMatricesToSolverMatrixMap, &
    & jacobianEquationsMatrixIdx,jacobianEquationsMatrixToSolverMatrixMap,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesToSolverMatrixMapType), POINTER :: equationsMatricesToSolverMatrixMap !<A pointer to the equations matrices to solver matrix map to get the Jacobian equations matrix to solver matrix map for
    INTEGER(INTG), INTENT(IN) :: jacobianEquationsMatrixIdx !<The Jacobian equations matrix index to get the Jacobian equations matrix to solver matrix for
    TYPE(EquationsMatrixToSolverMatrixMapType), POINTER :: jacobianEquationsMatrixToSolverMatrixMap !<On exit, a pointer to the specified Jacobian equations matrix to solver matrix map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMappingEMSToSMMap_JacobianEquationsMatrixToSolverMatrixGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(jacobianEquationsMatrixToSolverMatrixMap)) &
      & CALL FlagError("Jacobian equations matrix to solver matrix map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsMatricesToSolverMatrixMap)) &
      & CALL FlagError("Equations matrices to solver matrix map is not associated.",err,error,*999)
    IF(jacobianEquationsMatrixIdx<1.OR. &
      & jacobianEquationsMatrixIdx>equationsMatricesToSolverMatrixMap%numberOfJacobianEquationsMatrices) THEN
      localError="The specified Jacobian equations matrix index of "// &
        & TRIM(NumberToVString(jacobianEquationMatrixIdx,"*",err,error))//" is invalid for the equations matrices to solver "// &
        & "matrix map. The Jacobian equations matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(equationsMatricesToSolverMatrixMap%numberOfJacobianEquationsMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.MOT.ALLOCATED(equationsMatricesToSolverMatrixMap%jacobianEquationsMatrixToSolverMatrixMaps)) &
      & CALL FlagError("The Jacobian equations to solver matrix maps array is not allocated for the equations matrices to "// &
      & "solver matrix map.",err,error,*999)
#endif    

    jacobianEquationsMatrixToSolverMap=>equationsMatricesToSolverMatrixMap% &
      & jacobianEquationsMatrixToSolverMatrixMaps(jacobianEquationsMatrixIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(jacobianEquationsMatrix)) THEN
      localError="The Jacobian equations matrix to solver map is not associated for Jacobian equations matrix index "// &
        & TRIM(NumberToVString(jacobianEquationsMatrix,"*",err,error))// &
        & " of the equations matrices to solver matrix map."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("SolverMappingEMSToSMMap_JacobianEquationsMatrixToSolverMatrixMapGet")
    RETURN
999 NULLIFY(jacobianEquationsMatrix)
998 ERRORS("SolverMappingEMSToSMMap_JacobianEquationsMatrixToSolverMatrixMapGet",err,error)
    EXITS("SolverMappingEMSToSMMap_JacobianEquationsMatrixToSolverMatrixMapGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingEMSToSMMap_JacobianEquationsMatrixToSolverMatrixMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the linear equations matrix to solver matrix map for an equations matrices to solver matrix map.
  SUBROUTINE SolverMappingEMSToSMMap_LinearEquationsMatrixToSolverMatrixMapGet(equationsMatricesToSolverMatrixMap, &
    & linearEquationsMatrixIdx,linearEquationsMatrixToSolverMatrixMap,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesToSolverMatrixMapType), POINTER :: equationsMatricesToSolverMatrixMap !<A pointer to the equations matrices to solver matrix map to get the linear equations matrix to solver matrix map for
    INTEGER(INTG), INTENT(IN) :: linearEquationsMatrixIdx !<The linear equations matrix index to get the linear equations matrix to solver matrix for
    TYPE(EquationsMatrixToSolverMatrixMapType), POINTER :: linearEquationsMatrixToSolverMatrixMap !<On exit, a pointer to the specified linear equations matrix to solver matrix map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMappingEMSToSMMap_LinearEquationsMatrixToSolverMatrixGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(linearEquationsMatrixToSolverMatrixMap)) &
      & CALL FlagError("Linear equations matrix to solver matrix map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsMatricesToSolverMatrixMap)) &
      & CALL FlagError("Equations matrices to solver matrix map is not associated.",err,error,*999)
    IF(linearEquationsMatrixIdx<1.OR. &
      & linearEquationsMatrixIdx>equationsMatricesToSolverMatrixMap%numberOfLinearEquationsMatrices) THEN
      localError="The specified linear equations matrix index of "// &
        & TRIM(NumberToVString(linearEquationMatrixIdx,"*",err,error))//" is invalid for the equations matrices to solver "// &
        & "matrix map. The linear equations matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(equationsMatricesToSolverMatrixMap%numberOfLinearEquationsMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.MOT.ALLOCATED(equationsMatricesToSolverMatrixMap%linearEquationsMatrixToSolverMatrixMaps)) &
      & CALL FlagError("The linear equations to solver matrix maps array is not allocated for the equations matrices to "// &
      & "solver matrix map.",err,error,*999)
#endif    

    linearEquationsMatrixToSolverMap=>equationsMatricesToSolverMatrixMap% &
      & linearEquationsMatrixToSolverMatrixMaps(linearEquationsMatrixIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(linearEquationsMatrix)) THEN
      localError="The linear equations matrix to solver map is not associated for linear equations matrix index "// &
        & TRIM(NumberToVString(linearEquationsMatrix,"*",err,error))// &
        & " of the equations matrices to solver matrix map."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("SolverMappingEMSToSMMap_LinearEquationsMatrixToSolverMatrixMapGet")
    RETURN
999 NULLIFY(linearEquationsMatrix)
998 ERRORS("SolverMappingEMSToSMMap_LinearEquationsMatrixToSolverMatrixMapGet",err,error)
    EXITS("SolverMappingEMSToSMMap_LinearEquationsMatrixToSolverMatrixMapGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingEMSToSMMap_LinearEquationsMatrixToSolverMatrixMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of dynamic matrices in an equations matrices to solver matrix map.
  SUBROUTINE SolverMappingEMSToSMMap_NumberOfDynamicMatricesGet(equationsMatricesToSolverMatrixMap,numberOfDynamicMatrices, &
    & err,error,*)

    !Argument variables
    TYPE(EquationsMatricesToSolverMatrixMapType), POINTER :: equationsMatricesToSolverMatrixMap !<A pointer to the equations matrices to solver matrix map to get the number of dynamic matrices for
    INTEGER(INTG), INTENT(OUT) :: numberOfDynamicMatrices !<On exit, the number of dynamic matrices mapped to the solver matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingEMSToSMMap_NumberOfDynamicMatricesGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(equationsMatricesToSolverMatrixMap)) &
      & CALL FlagError("Equations matrices to solver matrix map is not associated.",err,error,*999)
#endif    
    
    numberOfDynamicMatrices=equationsMatricesToSolverMatrixMap%numberOfDynamicEquationsMatrices
      
    EXITS("SolverMappingEMSToSMMap_NumberOfDynamicMatricesGet")
    RETURN
999 ERRORSEXITS("SolverMappingEMSToSMMap_NumberOfDynamicMatricesGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMappingEMSToSMMap_NumberOfDynamicMatricesGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of Jacobian matrices in an equations matrices to solver matrix map.
  SUBROUTINE SolverMappingEMSToSMMap_NumberOfJacobianMatricesGet(equationsMatricesToSolverMatrixMap,numberOfJacobianMatrices, &
    & err,error,*)

    !Argument variables
    TYPE(EquationsMatricesToSolverMatrixMapType), POINTER :: equationsMatricesToSolverMatrixMap !<A pointer to the equations matrices to solver matrix map to get the number of jacobian matrices for
    INTEGER(INTG), INTENT(OUT) :: numberOfJacobianMatrices !<On exit, the number of Jacobian matrices mapped to the solver matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingEMSToSMMap_NumberOfJacobianMatricesGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(equationsMatricesToSolverMatrixMap)) &
      & CALL FlagError("Equations matrices to solver matrix map is not associated.",err,error,*999)
#endif    
    
    numberOfJacobianMatrices=equationsMatricesToSolverMatrixMap%numberOfJacobianEquationsMatrices
      
    EXITS("SolverMappingEMSToSMMap_NumberOfJacobianMatricesGet")
    RETURN
999 ERRORSEXITS("SolverMappingEMSToSMMap_NumberOfJacobianMatricesGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMappingEMSToSMMap_NumberOfJacobianMatricesGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of linear matrices in an equations matrices to solver matrix map.
  SUBROUTINE SolverMappingEMSToSMMap_NumberOfLinearMatricesGet(equationsMatricesToSolverMatrixMap,numberOfLinearMatrices, &
    & err,error,*)

    !Argument variables
    TYPE(EquationsMatricesToSolverMatrixMapType), POINTER :: equationsMatricesToSolverMatrixMap !<A pointer to the equations matrices to solver matrix map to get the number of linear matrices for
    INTEGER(INTG), INTENT(OUT) :: numberOfLinearMatrices !<On exit, the number of linear matrices mapped to the solver matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingEMSToSMMap_NumberOfLinearMatricesGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(equationsMatricesToSolverMatrixMap)) &
      & CALL FlagError("Equations matrices to solver matrix map is not associated.",err,error,*999)
#endif    
    
    numberOfLinearMatrices=equationsMatricesToSolverMatrixMap%numberOfLinearEquationsMatrices
      
    EXITS("SolverMappingEMSToSMMap_NumberOfLinearMatricesGet")
    RETURN
999 ERRORSEXITS("SolverMappingEMSToSMMap_NumberOfLinearMatricesGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMappingEMSToSMMap_NumberOfLinearMatricesGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the equations col to solver cols map for an equations matrix to solver matrix map.
  SUBROUTINE SolverMappingEMToSMMap_EquationsColToSolverColsMapGet(equationsMatrixToSolverMatrixMap,
    & equationsColToSolverColsMap,err,error,*)

    !Argument variables
    TYPE(EquationsMatrixToSolverMatrixMapType), POINTER :: equationsMatrixToSolverMatrixMap !<A pointer to the equations matrix to solver matrix map to get the equations col to solver cols map for
    TYPE(MatrixRowColCouplingType), POINTER :: equationsColToSolverColsMap(:) !<On exit, a pointer to the specified equations col to solver cols map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverMappingEMToSMMap_EquationsColToSolverColsMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(equationsColToSolverColsMap)) &
      & CALL FlagError("Equations to solver cols map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsMatrixToSolverMatrixMap)) &
      & CALL FlagError("Equations to matrix solver matrix map is not associated.",err,error,*999)
#endif    

    equationsColToSolverColsMap=>equationsMatrixToSolverMatrixMap%equationsColToSolverColsMap

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(equationsColToSolverColsMap)) CALL FlagError("The equations col to solver cols map for the specified "// &
      & "equations matrix to solver matrix map is not associated.",err,error,*999)
#endif    
      
    EXITS("SolverMappingEMToSMMap_EquationsColToSolverColsMapGet")
    RETURN
999 NULLIFY(equationsColToSolverColsMap)
998 ERRORS("SolverMappingEMToSMMap_EquationsColToSolverColsMapGet",err,error)
    EXITS("SolverMappingEMToSMMap_EquationsColToSolverColsMapGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingEMToSMMap_EquationsColToSolverColsMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the equations matrix for an equations to solver map.
  SUBROUTINE SolverMappingEMToSMMap_EquationsMatrixGet(equationsMatrixToSolverMatrixMap,equationsMatrix,err,error,*)

    !Argument variables
    TYPE(EquationsMatrixToSolverMatrixMapType), POINTER :: equationsMatrixToSolverMatrixMap !<A pointer to the equations matrix to solver matrix map to get the equations matrix for
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix !<On exit, a pointer to the specified equations matrix. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingEMToSMMap_EquationsMatrixGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(equationsMatrix)) CALL FlagError("Equations matrix is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsMatrixToSolverMatrixMap)) &
      & CALL FlagError("Equations matrix to solver matrix map is not associated.",err,error,*999)
#endif    

    equationsMatrix=>equationsMatrixToSolverMatrixMap%equationsMatrix

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(equationsMatrix)) &
      & CALL FlagError("The Equations matrix for the specified equations matrix to solver matrix map is not associated.", &
      & err,error,*999)
#endif    
      
    EXITS("SolverMappingEMToSMMap_EquationsMatrixGet")
    RETURN
999 NULLIFY(equationsMatrix)
998 ERRORS("SolverMappingEMToSMMap_EquationsMatrixGet",err,error)
    EXITS("SolverMappingEMToSMMap_EquationsMatrixGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingEMToSMMap_EquationsMatrixGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the solver matrix for an equations matrix to solver matrix map.
  SUBROUTINE SolverMappingEMToSMMap_SolverMatrixGet(equationsMatrixToSolverMatrixMap,solverMatrix,err,error,*)

    !Argument variables
    TYPE(EquationsMatrixToSolverMatrixMapType), POINTER :: equationsMatrixToSolverMatrixMap !<A pointer to the equations matrix to solver matrix map to get the solver matrix for
    TYPE(SolverMatrixType), POINTER :: solverMatrix !<On exit, a pointer to the specified solver matrix. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("_SolverMatrixGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(solverMatrix)) CALL FlagError("Solver matrix is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsMatrixToSolverMatrixMap)) &
      & CALL FlagError("Equations matrix to solver matrix map is not associated.",err,error,*999)
#endif    

    solverMatrix=>equationsMatrixToSolverMatrixMap%solverMatrix

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(solverMatrix)) &
      & CALL FlagError("The solver matrix for the specified equations matrix to solver matrix map is not associated.", &
      & err,error,*999)
#endif
      
    EXITS("SolverMappingEMToSMMap_SolverMatrixGet")
    RETURN
999 NULLIFY(solverMatrix)
998 ERRORS("SolverMappingEMToSMMap_SolverMatrixGet",err,error)
    EXITS("SolverMappingEMToSMMap_SolverMatrixGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingEMToSMMap_SolverMatrixGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the equations matrices to solver matrix map for a solver mapping equations set to solver matrices map.
  SUBROUTINE SolverMappingESToSMSMap_EquationsMatricesToSolverMatrixMapGet(equationsSetToSolverMatricesMap,solverMatrixIdx, &
    & equationsMatricesToSolverMatrixMap,err,error,*)

    !Argument variables
    TYPE(EquationsSetToSolverMatricesMapType), POINTER :: equationsSetToSolverMatricesMap !<A pointer to the solver mapping equations set to solver matrices map to get the equations matrices to solver matrix map for
    INTEGER(INTG), INTENT(IN) :: solverMatrixIdx !<The solver matrix index to get the equations matrices to solver matrix map for
    TYPE(EquationsMatricesToSolverMatrixMapType), POINTER :: equationsMatricesToSolverMatrixMap !<On exit, a pointer to the specified equations matrices to solver matrix map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMappingESToSMSMap_EquationsMatricesToSolverMatrixMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(equationsMatricesToSolverMatrixMap)) &
      & CALL FlagError("Equations matrices to solver matrix map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsSetToSolverMatricesMap))  &
      & CALL FlagError("Equations set to solver matrices map is not associated.",err,error,*999)
    IF(solverMatrixIdx<1.OR.solverMatrixIdx>equationsSetToSolverMatricesMap%numberOfSolverMatrices) THEN
      localError="The specified solver matrix index of "//TRIM(NumberToVString(solverMatrixIdx,"*",err,error))// &
        & " is invalid for the equations set to solver matrices map. The solver matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(equationsSetToSolverMatricesMap%numberOfSolverMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,eror,*999)
    ENDIF
    IF(.NOT.ALLOCATED(equationsSetToSolverMatricesMap%equationsMatricesToSolverMatrixMaps)) &
      & CALL FlagError("The equations matrices to solver matrix maps is not allocated for the equations set to "// &
      & "solver matrices map.",err,error,*999)
#endif    
    
    equationsMatricesToSolverMatrixMap=>equationsSetToSolverMatricesMap% &
      & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr
    
#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(equationsMatricesToSolverMatrixMap)) THEN
      localError="The equations matrices to solver matrix map is not associated for the solver matrix index of "// &
        & TRIM(NumberToVString(solverMatrixIdx,"*",err,error))// &
        & " of the equations set to solver matrices map of the solver mapping."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
      
    EXITS("SolverMappingESToSMSMap_EquationsMatricesToSolverMatrixMapGet")
    RETURN
999 NULLIFY(equationsMatricesToSolverMatrixMap)
998 ERRORS("SolverMappingESToSMSMap_EquationsMatricesToSolverMatrixMapGet",err,error)
    EXITS("SolverMappingESToSMSMap_EquationsMatricesToSolverMatrixMapGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingESToSMSMap_EquationsMatricesToSolverMatrixMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the equations matrix to solver matrices map for a solver mapping equations set to solver matrices map.
  SUBROUTINE SolverMappingESToSMSMap_EquationsMatrixToSolverMatricesMapGet(equationsSetToSolverMatricesMap,equationsMatrixIdx, &
    & equationsMatrixToSolverMatricesMap,err,error,*)

    !Argument variables
    TYPE(EquationsSetToSolverMatricesMapType), POINTER :: equationsSetToSolverMatricesMap !<A pointer to the solver mapping equations set to solver matrices map to get the equations matrix to solver matrices map for
    INTEGER(INTG), INTENT(IN) :: equationsMatrixIdx !<The equations matrix index in to get the equations matrix to solver matrices map for
    TYPE(EquationsMatrixToSolverMatricesMapType), POINTER :: equationsMatrixToSolverMatricesMap !<On exit, a pointer to the specified equations matrix to solver matrices map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMappingESToSMSMap_EquationsMatrixToSolverMatricesMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(equationsMatrixToSolverMatricesMap)) &
      & CALL FlagError("Equations matrix to solver matrices map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsSetToSolverMatricesMap))  &
      & CALL FlagError("Equations set to solver matrices map is not associated.",err,error,*999)
    IF(equationsMatrixIdx<1.OR.equationsMatrixIdx>equationsSetToSolverMatricesMap%numberOfEquationsMatrices) THEN
      localError="The specified equations matrix index of "//TRIM(NumberToVString(equationsMatrixIdx,"*",err,error))// &
        & " is invalid for the equations set to solver matrices map. The equations matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(equationsSetToSolverMatricesMap%numberOfEquationsMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,eror,*999)
    ENDIF
    IF(.NOT.ALLOCATED(equationsSetToSolverMatricesMap%equationsMatrixToSolverMatricesMaps)) &
      & CALL FlagError("The equations matrix to solver matrices maps is not allocated for the equations set to "// &
      & "solver matrices map.",err,error,*999)
#endif    
    
    equationsMatrixToSolverMatricesMap=>equationsSetToSolverMatricesMap% &
      & equationsMatrixToSolverMatricesMaps(equationsMatrixIdx)%ptr
    
#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(equationsMatrixToSolverMatricesMap)) THEN
      localError="The equations matrix to solver matrices map is not associated for the equations matrix index of "// &
        & TRIM(NumberToVString(equationsMatrixIdx,"*",err,error))// &
        & " of the equations set to solver matrices map of the solver mapping."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
      
    EXITS("SolverMappingESToSMSMap_EquationsMatrixToSolverMatricesMapGet")
    RETURN
999 NULLIFY(equationsMatrixToSolverMatricesMap)
998 ERRORS("SolverMappingESToSMSMap_EquationsMatrixToSolverMatricesMapGet",err,error)
    EXITS("SolverMappingESToSMSMap_EquationsMatrixToSolverMatricesMapGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingESToSMSMap_EquationsMatrixToSolverMatricesMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the Jacobian matrix to solver matrix map for a solver mapping equations set to solver matrices map.
  SUBROUTINE SolverMappingESToSMSMap_JacobianMatrixToSolverMatrixMapGet(equationsSetToSolverMatricesMap,jacobianMatrixIdx, &
    & jacobianMatrixToSolverMatrixMap,err,error,*)

    !Argument variables
    TYPE(EquationsSetToSolverMatricesMapType), POINTER :: equationsSetToSolverMatricesMap !<A pointer to the solver mapping equations set to solver matrices map to get the Jacobian matrix to solver matrix map for
    INTEGER(INTG), INTENT(IN) :: jacobianMatrixIdx !<The Jacobian matrix index in to get the Jacobian matrix to solver matrix map for
    TYPE(JacobianMatrixToSolverMatrixMapType), POINTER :: jacobianMatrixToSolverMatrixMap !<On exit, a pointer to the specified Jacobian matrix to solver matrix map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMappingESToSMSMap_JacobianMatrixToSolverMatrixMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(jacobianMatrixToSolverMatrixMap)) &
      & CALL FlagError("Jacobian matrix to solver matrix map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsSetToSolverMatricesMap))  &
      & CALL FlagError("Equations set to solver matrices map is not associated.",err,error,*999)
    IF(jacobianMatrixIdx<1.OR.jacobianMatrixIdx>equationsSetToSolverMatricesMap%numberOfJacobianMatrices) THEN
      localError="The specified Jacobian matrix index of "//TRIM(NumberToVString(jacobianMatrixIdx,"*",err,error))// &
        & " is invalid for the equations set to solver matrices map. The Jacobian matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(equationsSetToSolverMatricesMap%numberOfJacobianMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,eror,*999)
    ENDIF
    IF(.NOT.ALLOCATED(equationsSetToSolverMatricesMap%jacobianMatrixToSolverMatrixMaps)) &
      & CALL FlagError("The Jacobian matrix to solver matrix maps is not allocated for the equations set to "// &
      & "solver matrices map.",err,error,*999)
#endif    
    
    jacobianMatrixToSolverMatrixMap=>equationsSetToSolverMatricesMap% &
      & jacobianMatrixToSolverMatrixMaps(jacobianMatrixIdx)%ptr
    
#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(jacobianMatrixToSolverMatrixMap)) THEN
      localError="The equations matrix to solver matrix map is not associated for the Jacobian matrix index of "// &
        & TRIM(NumberToVString(jacobianMatrixIdx,"*",err,error))// &
        & " of the equations set to solver matrices map of the solver mapping."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
      
    EXITS("SolverMappingESToSMSMap_JacobianMatrixToSolverMatrixMapGet")
    RETURN
999 NULLIFY(jacobianMatrixToSolverMatrixMap)
998 ERRORS("SolverMappingESToSMSMap_JacobianMatrixToSolverMatrixapGet",err,error)
    EXITS("SolverMappingESToSMSMap_JacobianMatrixToSolverMatrixMapGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingESToSMSMap_JacobianMatrixToSolverMatrixMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of equations matrices for a solver mapping equations set to solver matrices map.
  SUBROUTINE SolverMappingESToSMSMap_NumberOfEquationsMatricesGet(equationsSetToSolverMatricesMap,numberOfEquationsMatrices, &
    & err,error,*)

    !Argument variables
    TYPE(EquationsSetToSolverMatricesMapType), POINTER :: equationsSetToSolverMatricesMap !<A pointer to the solver mapping equations set to solver matrices map to get the number of equations matrices for
    INTEGER(INTG), INTENT(OUT) :: numberOfEquationsMatrices !<On exit, the number of equations matrices in the equations set to solver matrices map. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingESToSMSMap_NumberOfEquationsMatricesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(equationsSetToSolverMatricesMap))  &
      & CALL FlagError("Equations set to solver matrices map is not associated.",err,error,*999)
#endif    
    
    numberOfEquationsMatrices=equationsSetToSolverMatricesMap%numberOfEquationsMatrices
      
    EXITS("SolverMappingESToSMSMap_NumberOfEquationsMatricesGet")
    RETURN
999 ERRORS("SolverMappingESToSMSMap_NumberOfEquationsMatricesGet",err,error)
    EXITS("SolverMappingESToSMSMap_NumberOfEquationsMatricesGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingESToSMSMap_NumberOfEquationsMatricesGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of Jacobian matrices for a solver mapping jacobian set to solver matrices map.
  SUBROUTINE SolverMappingESToSMSMap_NumberOfJacobianMatricesGet(jacobianSetToSolverMatricesMap,numberOfJacobianMatrices, &
    & err,error,*)

    !Argument variables
    TYPE(JacobianSetToSolverMatricesMapType), POINTER :: jacobianSetToSolverMatricesMap !<A pointer to the solver mapping equations set to solver matrices map to get the number of Jacobian matrices for
    INTEGER(INTG), INTENT(OUT) :: numberOfJacobianMatrices !<On exit, the number of Jacobian matrices in the equations set to solver matrices map. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingESToSMSMap_NumberOfJacobianMatricesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(equationsSetToSolverMatricesMap))  &
      & CALL FlagError("Equations set to solver matrices map is not associated.",err,error,*999)
#endif    
    
    numberOfJacobianMatrices=equationsSetToSolverMatricesMap%numberOfJacobianMatrices
      
    EXITS("SolverMappingESToSMSMap_NumberOfJacobianMatricesGet")
    RETURN
999 ERRORS("SolverMappingESToSMSMap_NumberOfJacobianMatricesGet",err,error)
    EXITS("SolverMappingESToSMSMap_NumberOfJacobianMatricesGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingESToSMSMap_NumberOfJacobianMatricesGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of solver matrices for a solver mapping equations set to solver matrices map.
  SUBROUTINE SolverMappingESToSMSMap_NumberOfSolverMatricesGet(equationsSetToSolverMatricesMap,numberOfSolverMatrices, &
    & err,error,*)

    !Argument variables
    TYPE(EquationsSetToSolverMatricesMapType), POINTER :: equationsSetToSolverMatricesMap !<A pointer to the solver mapping equations set to solver matrices map to get the number of solver matrices for
    INTEGER(INTG), INTENT(OUT) :: numberOfSolverMatrices !<On exit, the number of solver matrices in the equations set to solver matrices map. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingESToSMSMap_NumberOfSolverMatricesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(equationsSetToSolverMatricesMap))  &
      & CALL FlagError("Equations set to solver matrices map is not associated.",err,error,*999)
#endif    
    
    numberOfSolverMatrices=equationsSetToSolverMatricesMap%numberOfSolverMatrices
      
    EXITS("SolverMappingESToSMSMap_NumberOfSolverMatricesGet")
    RETURN
999 ERRORS("SolverMappingESToSMSMap_NumberOfSolverMatricesGet",err,error)
    EXITS("SolverMappingESToSMSMap_NumberOfSolverMatricesGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingESToSMSMap_NumberOfSolverMatricesGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the interface matrix to solver matrices map for a solver mapping interface condition to solver matrices map.
  SUBROUTINE SolverMappingICToSMSMap_InterfaceMatrixToSolverMatricesMapGet(interfaceConditionToSolverMatricesMap, &
    & interfaceMatrixIdx,interfaceMatrixToSolverMatricesMap,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionToSolverMatricesMapType), POINTER :: interfaceConditionToSolverMatricesMap !<A pointer to the solver mapping interface condition to solver matrices map to get the interface matrix to solver matrices map for
    INTEGER(INTG), INTENT(IN) :: interfaceMatrixIdx !<The interface matrix index in to get the interface matrix to solver matrices map for
    TYPE(InterfaceMatrixToSolverMatricesMapType), POINTER :: interfaceMatrixToSolverMatricesMap !<On exit, a pointer to the specified interface matrix to solver matrices map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMappingICToSMSMap_InterfaceMatrixToSolverMatricesMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(interfaceMatrixToSolverMatricesMap)) &
      & CALL FlagError("Interface matrix to solver matrices map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceConditionToSolverMatricesMap))  &
      & CALL FlagError("Interface condition to solver matrices map is not associated.",err,error,*999)
    IF(interfaceMatrixIdx<1.OR.interfaceMatrixIdx>interfaceConditionToSolverMatricesMap%nuberOfInterfaceMatrices)) THEN
      localError="The specified interface matrix index of "//TRIM(NumberToVString(interfaceMatrixIdx,"*",err,error))// &
        & " is invalid for the interface condition to solver matrices map. The interface matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(interfaceConditionToSolverMatricesMap%numberOfInterfaceMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,eror,*999)
    ENDIF
    IF(.NOT.ALLOCATED(interfaceConditionToSolverMatricesMap%interfaceMatrixToSolverMatricesMaps)) &
      & CALL FlagError("The interface matrix to solver matrices maps is not allocated for the interface condition to "// &
      & "solver matrices map.",err,error,*999)
#endif    
    
    interfaceMatrixToSolverMatricesMap=>interfaceConditionToSolverMatricesMap% &
      & interfaceMatrixToSolverMatricesMaps(interfaceMatrixIndex)%ptr
    
#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(interfaceMatrixToSolverMatricesMap)) THEN
      localError="The interface matrix to solver matrices map is not associated for the interface matrix index of "// &
        & TRIM(NumberToVString(interfaceMatrixIdx,"*",err,error))// &
        & " of the interface condition to solver matrices map of the solver mapping."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
      
    EXITS("SolverMappingICToSMSMap_InterfaceMatrixToSolverMatricesMapGet")
    RETURN
999 NULLIFY(interfaceMatrixToSolverMatricesMap)
998 ERRORS("SolverMappingICToSMSMap_InterfaceMatrixToSolverMatricesMapGet",err,error)
    EXITS("SolverMappingICToSMSMap_InterfaceMatrixToSolverMatricesMapGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingICToSMSMap_InterfaceMatrixToSolverMatricesMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the interface matrix for an interface matrix to solver matrix map.
  SUBROUTINE SolverMappingIMToSMMap_InterfaceMatrixGet(interfaceMatrixToSolverMatrixMap,interfaceMatrix,err,error,*)

    !Argument variables
    TYPE(InterfaceMatrixToSolverMatrixMapType), POINTER :: interfaceMatrixToSolverMatrixMap !<A pointer to the interface matrix to solver matrix map to get the interface matrix for
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix !<On exit, a pointer to the specified interface matrix. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingIMToSMMap_InterfaceMatrixGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(interfaceMatrix)) CALL FlagError("Interface matrix is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMatrixToSolverMatrixMap)) &
      & CALL FlagError("Interface matrix to solver matrix map is not associated.",err,error,*999)
#endif    

    interfaceMatrix=>interfaceMatrixToSolverMatrixMap%interfaceMatrix

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(interfaceMatrix)) &
      & CALL FlagError("The Interface matrix for the specified interface matrix to solver matrix map is not associated.", &
      & err,error,*999)
#endif    
      
    EXITS("SolverMappingIMToSMMap_InterfaceMatrixGet")
    RETURN
999 NULLIFY(interfaceMatrix)
998 ERRORS("SolverMappingIMToSMMap_InterfaceMatrixGet",err,error)
    EXITS("SolverMappingIMToSMMap_InterfaceMatrixGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingIMToSMMap_InterfaceMatrixGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the interface row to solver cols map for an interface matrix to solver matrix map.
  SUBROUTINE SolverMappingIMToSMMap_InterfaceRowToSolverColsMapGet(interfaceMatrixToSolverMatrixMap,
    & interfaceRowToSolverColsMap,err,error,*)

    !Argument variables
    TYPE(InterfaceMatrixToSolverMatrixMapType), POINTER :: interfaceMatrixToSolverMatrixMap !<A pointer to the interface matrix to solver matrix map to get the interface row to solver cols map for
    TYPE(MatrixRowColCouplingType), POINTER :: interfaceRowToSolverColsMap(:) !<On exit, a pointer to the specified interface row to solver cols map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingIMToSMMap_InterfaceRowToSolverColsMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(interfaceRowToSolverColsMap)) &
      & CALL FlagError("Interface row to solver cols map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMatrixToSolverMatrixMap)) &
      & CALL FlagError("Interface matrix to solver matrix map is not associated.",err,error,*999)
#endif    

    interfaceRowToSolverColsMap=>interfaceMatrixToSolverMatrixMap%interfaceRowToSolverColsMap

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(interfaceRowToSolverColsMap)) &
      & CALL FlagError("The interface row to solver cols map for the specified interface to solver map is not associated.", &
      & err,error,*999)
#endif    
      
    EXITS("SolverMappingIMToSMMap_InterfaceRowToSolverColsMapGet")
    RETURN
999 NULLIFY(interfaceRowToSolverColsMap)
998 ERRORS("SolverMappingIMToSMMap_InterfaceRowToSolverColsMapGet",err,error)
    EXITS("SolverMappingIMToSMMap_InterfaceRowToSolverColsMapGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingIMToSMMap_InterfaceRowToSolverColsMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the solver matrix for an interface matrix to solver matrix map.
  SUBROUTINE SolverMappingIMToSMMap_SolverMatrixGet(interfaceToSolverMap,solverMatrix,err,error,*)

    !Argument variables
    TYPE(InterfaceMatrixToSolverMatrixMapType), POINTER :: interfaceMatrixToSolverMatrixMap !<A pointer to the interface matrix to solver matrix map to get the solver matrix for
    TYPE(SolverMatrixType), POINTER :: solverMatrix !<On exit, a pointer to the specified solver matrix. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingIMToSMMap_SolverMatrixGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(solverMatrix)) CALL FlagError("Solver matrix is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMatrixToSolverMatrixMap)) &
      & CALL FlagError("Interface to solver map is not associated.",err,error,*999)
#endif    

    solverMatrix=>interfaceMatrixToSolverMatrixMap%solverMatrix

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(solverMatrix)) &
      & CALL FlagError("The solver matrix for the specified interface matrix to solver matrix map is not associated.", &
      & err,error,*999)
#endif    
      
    EXITS("SolverMappingIMToSMMap_SolverMatrixGet")
    RETURN
999 NULLIFY(solverMatrix)
998 ERRORS("SolverMappingIMToSMMap_SolverMatrixGet",err,error)
    EXITS("SolverMappingIMToSMMap_SolverMatrixGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingIMToSMMap_SolverMatrixGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the Jacobian col to solver cols map for an Jacobian matrix to solver matrix map.
  SUBROUTINE SolverMappingJMToSMMap_JacobianColToSolverColsMapGet(jacobianMatrixToSolverMatrixMap,
    & jacobianColToSolverColsMap,err,error,*)

    !Argument variables
    TYPE(JacobianMatrixToSolverMatrixMapType), POINTER :: jacobianMatrixToSolverMatrixMap !<A pointer to the Jacobian matrix to solver matrix map to get the Jacobian col to solver cols map for
    TYPE(MatrixRowColCouplingType), POINTER :: jacobianColToSolverColsMap(:) !<On exit, a pointer to the specified Jacobian col to solver cols map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingJMToSMMap_JacobianColToSolverColsMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(jacobianColToSolverColsMap)) &
      & CALL FlagError("Jacobian to solver cols map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(jacobianMatrixToSolverMatrixMap)) &
      & CALL FlagError("Jacobian matrix to solver matrix map is not associated.",err,error,*999)
#endif    

    jacobianColToSolverColsMap=>jacobianMatrixToSolverMatrixMap%jacobianColToSolverColsMap

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(jacobianColToSolverColsMap)) CALL FlagError("The Jacobian col to solver cols map for the specified "// &
      & "Jacobian matrix to solver matrix map is not associated.",err,error,*999)
#endif    
      
    EXITS("SolverMappingJMToSMMap_JacobianColToSolverColsMapGet")
    RETURN
999 NULLIFY(jacobianColToSolverColsMap)
998 ERRORS("SolverMappingJMToSMMap_JacobianColToSolverColsMapGet",err,error)
    EXITS("SolverMappingJMToSMMap_JacobianColToSolverColsMapGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingJMToSMMap_JacobianColToSolverColsMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the Jacobian matrix for an Jacobian matrix to solver matrix map.
  SUBROUTINE SolverMappingJMToSMMap_JacobianMatrixGet(jacobianMatrixToSolverMatrixMap,jacobianMatrix,err,error,*)

    !Argument variables
    TYPE(JacobianMatrixToSolverMatrixMapType), POINTER :: jacobianMatrixToSolverMatrixMap !<A pointer to the Jacobian matrix to solver matrix map to get the Jacobian matrix for
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix !<On exit, a pointer to the specified Jacobian matrix. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingJMToSMMap_JacobianMatrixGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(jacobianMatrix)) CALL FlagError("Jacobian matrix is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(jacobianMatrixToSolverMatrixMap)) &
      & CALL FlagError("Jacobian matrix to solver matrix map is not associated.",err,error,*999)
#endif    

    jacobianMatrix=>jacobianMatrixToSolverMatrixMap%jacobianMatrix

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(jacobianMatrix)) CALL FlagError("The Jacobian matrix for the specified Jacobian matrix to solver "// &
      & "matrix map is not associated.",err,error,*999)
#endif    
      
    EXITS("SolverMappingJMToSMMap_JacobianMatrixGet")
    RETURN
999 NULLIFY(jacobianMatrix)
998 ERRORS("SolverMappingJMToSMMap_JacobianMatrixGet",err,error)
    EXITS("SolverMappingJMToSMMap_JacobianMatrixGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingJMToSMMap_JacobianMatrixGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the solver matrix for an Jacobian matrix to solver matrix map.
  SUBROUTINE SolverMappingJMToSMMap_SolverMatrixGet(jacobianMatrixToSolverMatrixMap,solverMatrix,err,error,*)

    !Argument variables
    TYPE(JacobianMatrixToSolverMatrixMapType), POINTER :: jacobianMatrixToSolverMatrixMap !<A pointer to the Jacobian matrix to solver matrix map to get the solver matrix for
    TYPE(SolverMatrixType), POINTER :: solverMatrix !<On exit, a pointer to the specified solver matrix. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingJMToSMMap_SolverMatrixGet",err,error,*998)

#ifdef WITH_PRECHECKS
    IF(ASSOCIATED(solverMatrix)) CALL FlagError("Solver matrix is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(jacobianMatrixToSolverMatrixMap)) &
      & CALL FlagError("Jacobian matrix to solver matrix map is not associated.",err,error,*999)
#endif    

    solverMatrix=>jacobianMatrixToSolverMatrixMap%solverMatrix

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(solverMatrix)) CALL FlagError("The solver matrix for the specified Jacobian matrix to solver matrix "// &
      & "map is not associated.",err,error,*999)
#endif    
      
    EXITS("SolverMappingJMToSMMap_SolverMatrixGet")
    RETURN
999 NULLIFY(solverMatrix)
998 ERRORS("SolverMappingJMToSMMap_SolverMatrixGet",err,error)
    EXITS("SolverMappingJMToSMMap_SolverMatrixGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingJMToSMMap_SolverMatrixGet
  
  !
  !================================================================================================================================
  !

END MODULE SolverMappingAccessRoutines
