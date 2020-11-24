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
  !> \brief Equations types for solver mapping
  !> \see SolverMappingRoutines
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_MAPPING_EQUATIONS_EQUATIONS_SET=1 !<The equations in the solver mapping is from an equations set \see SolverMappingRoutines_EquationsTypes,SolverMappingRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_MAPPING_EQUATIONS_INTERFACE_CONDITION=2 !<The equations in the solver mapping is from an interface condition \see SolverMappingRoutines_EquationsTypes,SolverMappingRoutines
  INTEGER(INTG), PARAMETER :: SOLVER_MAPPING_EQUATIONS_INTERFACE_TRANSPOSE=3 !<The equations in the solver mapping is from a transposed interface condition \see SolverMappingRoutines_EquationsTypes,SolverMappingRoutines
  !>@}
 
  !Module types

  !Module variables 

  !Interfaces

  PUBLIC SOLVER_MAPPING_EQUATIONS_DYNAMIC_MATRIX,SOLVER_MAPPING_EQUATIONS_LINEAR_MATRIX,SOLVER_MAPPING_EQUATIONS_NONLINEAR_MATRIX

  PUBLIC SOLVER_MAPPING_EQUATIONS_EQUATIONS_SET,SOLVER_MAPPING_EQUATIONS_INTERFACE_CONDITION, &
    & SOLVER_MAPPING_EQUATIONS_INTERFACE_TRANSPOSE

  PUBLIC SolverMapping_AssertIsFinished,SolverMapping_AssertNotFinished

  PUBLIC SolverMapping_ColumnDOFSMappingGet

  PUBLIC SolverMapping_CreateValuesCacheGet
  
  PUBLIC SolverMapping_DynamicMatrixToSolverMatrixGet
  
  PUBLIC SolverMapping_EMSToSMMapNumberOfDynamicMatricesGet

  !PUBLIC SolverMapping_IMSToSMMapNumberOfJacobianMatricesGet

  PUBLIC SolverMapping_EMSToSMMapNumberOfLinearMatricesGet

  PUBLIC SolverMapping_EMSToSMMapNumberOfVariablesGet

  PUBLIC SolverMapping_EMSToSMMapVariableGet

  PUBLIC SolverMapping_EMSToSMMapVariableDOFToSolverDOFsMapGet

  PUBLIC SolverMapping_EquationsSetGet

  PUBLIC SolverMapping_EquationsRowToSolverRowsMapGet

  PUBLIC SolverMapping_EquationsMatrixToSolverMatrixMapGet
  
  PUBLIC SolverMapping_EquationsSetToSolverMatricesMapGet

  PUBLIC SolverMapping_IMSToSMMapLagrangeVarDOFToSolverDOFsMapGet

  PUBLIC SolverMapping_IMSToSMMapDependentVarDOFToSolverDOFsMapGet

  PUBLIC SolverMapping_InterfaceColToSolverColsMapGet

  PUBLIC SolverMapping_InterfaceColToSolverRowsMapGet

  PUBLIC SolverMapping_InterfaceConditionGet

  PUBLIC SolverMapping_InterfaceConditionToSolverMatricesMapGet

  PUBLIC SolverMapping_InterfaceEquationsMatrixToSolverMatrixGet
  
  PUBLIC SolverMapping_InterfaceMatrixToSolverMatrixMapGet 

  PUBLIC SolverMapping_InterfaceRowToSolverRowsMapGet

  PUBLIC SolverMapping_JacobianEquationsMatrixToSolverMatrixGet
  
  PUBLIC SolverMapping_JacobianMatrixToSolverMatrixMapGet

  PUBLIC SolverMapping_LinearMatrixToSolverMatrixGet

  PUBLIC SolverMapping_NumberOfEquationsSetsGet

  PUBLIC SolverMapping_NumberOfInterfaceConditionsGet

  PUBLIC SolverMapping_NumberOfSolverMatricesGet

  PUBLIC SolverMapping_NumberOfRowsGet

  PUBLIC SolverMapping_NumberOfGlobalRowsGet

  PUBLIC SolverMapping_RHSVariablesListGet
  
  PUBLIC SolverMapping_RowDOFSMappingGet

  PUBLIC SolverMapping_SolverMatrixToEquationsMapGet

  PUBLIC SolverMapping_SolverRowToEquationsMapGet

  PUBLIC SolverMappingCVC_EquationsVariableListGet
  
  PUBLIC SolverMappingCVC_InterfaceVariableListGet
  
  PUBLIC SolverMappingCVC_RHSVariableListGet
  
  PUBLIC SolverMappingEMSToSMMap_DynamicMatrixToSolverMatrixMapGet

  PUBLIC SolverMappingEMSToSMMap_JacobianMatrixToSolverMatrixMapGet

  PUBLIC SolverMappingEMSToSMMap_LinearMatrixToSolverMatrixMapGet

  PUBLIC SolverMappingEMSToSMMap_NumberOfDynamicMatricesGet

  PUBLIC SolverMappingEMSToSMMap_NumberOfJacobianMatricesGet
  
  PUBLIC SolverMappingEMSToSMMap_NumberOfLinearMatricesGet

  PUBLIC SolverMappingEMSToSMMap_NumberOfVariablesGet

  PUBLIC SolverMappingEMSToSMMap_VariableGet

  PUBLIC SolverMappingEMSToSMMap_VariableDOFToSolverDOFsMapGet

  PUBLIC SolverMappingEMToSMMap_EquationsColToSolverColsMapGet

  PUBLIC SolverMappingEMToSMMap_EquationsMatrixGet

  PUBLIC SolverMappingEMToSMMap_SolverMatrixGet

  PUBLIC SolverMappingEMToSMSMap_EquationsMatrixToSolverMatrixMapGet

  PUBLIC SolverMappingEMToSMSMap_NumberOfSolverMatricesGet

  PUBLIC SolverMappingESToSMSMap_EquationsMatricesToSolverMatrixMapGet

  PUBLIC SolverMappingESToSMSMap_EquationsMatrixToSolverMatricesMapGet

  PUBLIC SolverMappingESToSMSMap_EquationsRowToSolverRowsMapGet

  PUBLIC SolverMappingESToSMSMap_JacobianMatrixToSolverMatrixMapGet

  PUBLIC SolverMappingESToSMSMap_NumberOfEquationsMatricesGet
  
  PUBLIC SolverMappingESToSMSMap_NumberOfJacobianMatricesGet

  PUBLIC SolverMappingESToSMSMap_NumberOfSolverMatricesGet

  PUBLIC SolverMappingESToSMSMap_SolverMappingGet

  PUBLIC SolverMappingICToSMSMap_InterfaceColToSolverRowsMapGet

  PUBLIC SolverMappingICToSMSMap_InterfaceMatricesToSolverMatrixMapGet

  PUBLIC SolverMappingIMToSMSMap_InterfaceMatrixToSolverMatrixMapGet

  PUBLIC SolverMappingIMToSMSMap_InterfaceRowToSolverRowsMapGet

  PUBLIC SolverMappingIMToSMSMap_NumberOfSolverMatricesGet

  PUBLIC SolverMappingIMSToSMMap_DependentVariableGet
  
  PUBLIC SolverMappingIMSToSMMap_DependentVarDOFToSolverDOFsMapGet
  
  PUBLIC SolverMappingIMSToSMMap_InterfaceColToSolverColsMapGet
  
  PUBLIC SolverMappingIMSToSMMap_InterfaceMatrixToSolverMatrixMapGet
  
  PUBLIC SolverMappingIMSToSMMap_LagrangeVariableGet

  PUBLIC SolverMappingIMSToSMMap_LagrangeVarDOFToSolverDOFsMapGet
  
  PUBLIC SolverMappingIMSToSMMap_NumberOfDependentVariablesGet
  
  PUBLIC SolverMappingIMSToSMMap_NumberOfInterfaceMatricesGet
  
  PUBLIC SolverMappingIMToSMMap_InterfaceMatrixGet

  PUBLIC SolverMappingIMToSMMap_InterfaceRowToSolverColsMapGet

  PUBLIC SolverMappingIMToSMMap_SolverMatrixGet

  PUBLIC SolverMappingJMToSMMap_JacobianColToSolverColsMapGet

  PUBLIC SolverMappingJMToSMMap_JacobianMatrixGet
  
  PUBLIC SolverMappingJMToSMMap_SolverMatrixGet

  PUBLIC SolverMappingSColToDEQSMap_DynamicCouplingInfoGet
  
  PUBLIC SolverMappingSColToDEQSMap_NumberOfDynamicMatricesGet
  
  PUBLIC SolverMappingSColToIEQSMap_InterfaceCouplingInfoGet
  
  PUBLIC SolverMappingSColToIEQSMap_NumberOfInterfaceMatricesGet
  
  PUBLIC SolverMappingSColToLEQSMap_LinearCouplingInfoGet
  
  PUBLIC SolverMappingSColToLEQSMap_NumberOfLinearMatricesGet
  
  PUBLIC SolverMappingSColToNLEQSMap_JacobianCouplingInfoGet
  
  PUBLIC SolverMappingSColToNLEQSMap_NumberOfJacobianMatricesGet
  
  PUBLIC SolverMappingSDOFToVDOFsMap_EquationsInfoGet

  PUBLIC SolverMappingSDOFToVDOFsMap_NumberOfEquationDOFsGet

  PUBLIC SolverMappingSDOFToVDOFsMap_VariableDOFCouplingGet

  PUBLIC SolverMappingSDOFToVDOFsMap_VariableGet

  PUBLIC SolverMappingSMToEQSMap_ColumnDOFsMappingGet

  PUBLIC SolverMappingSMToEQSMap_NumberOfColumnsGet

  PUBLIC SolverMappingSMToEQSMap_NumberOfDOFsGet
 
  PUBLIC SolverMappingSMToEQSMap_NumberOfEquationsSetsGet

  PUBLIC SolverMappingSMToEQSMap_NumberOfGlobalDOFsGet

  PUBLIC SolverMappingSMToEQSMap_NumberOfInterfaceConditionsGet

  PUBLIC SolverMappingSMToEQSMap_SolverDOFToVariableDOFsMapGet
  
  PUBLIC SolverMappingSMToEQSMap_SolverMappingGet

  PUBLIC SolverMappingSMToEQSMap_SolverMatrixGet

  PUBLIC SolverMappingSMToEQSMap_SolverMatrixToEquationsSetMapGet

  PUBLIC SolverMappingSMToEQSMap_SolverMatrixToInterfaceConditionMapGet

  PUBLIC SolverMappingSMToEQSMap_TotalNumberOfDOFsGet
 
  PUBLIC SolverMappingSMToEQSMap_VariablesListGet

  PUBLIC SolverMappingSMToESMap_EquationsGet
  
  PUBLIC SolverMappingSMToESMap_HaveEquationsGet

  PUBLIC SolverMappingSMToESMap_SolverColToDynamicEquationsMapGet

  PUBLIC SolverMappingSMToESMap_SolverColToLinearEquationsMapGet

  PUBLIC SolverMappingSMToESMap_SolverColToNonlinearEquationsMapGet

  PUBLIC SolverMappingSMToICMap_InterfaceEquationsGet

  PUBLIC SolverMappingSMToICMap_SolverColToInterfaceEquationsMapGet

  PUBLIC SolverMappingSRowToEQSMap_EquationsSetCouplingInfoGet

  PUBLIC SolverMappingSRowToEQSMap_InterfaceConditionCouplingInfoGet
  
  PUBLIC SolverMappingSRowToEQSMap_InterfaceConditionIndexGet

  PUBLIC SolverMappingSRowToEQSMap_NumberOfEquationsSetRowsGet

  PUBLIC SolverMappingVariable_EquationInfoGet

  PUBLIC SolverMappingVariable_EquationMatrixNumberGet

  PUBLIC SolverMappingVariable_EquationNumberOfMatricesGet

  PUBLIC SolverMappingVariable_FieldVariableGet

  PUBLIC SolverMappingVariable_NumberOfEquationsGet

  PUBLIC SolverMappingVariables_NumberOfVariablesGet

  PUBLIC SolverMappingVariables_VariableGet

  PUBLIC SolverMappingVariables_VariableInListCheck

  PUBLIC SolverMappingVDOFToSDOFsMap_SolverDOFCouplingGet

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
  SUBROUTINE SolverMapping_ColumnDOFSMappingGet(solverMapping,matrixIdx,columnDOFSMapping,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the column DOFs mapping for
    INTEGER(INTG), INTENT(IN) :: matrixIdx !<The matrix index in the solver mapping to get the column DOFs mapping for
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
    IF(matrixIdx<1.OR.matrixIdx>solverMapping%numberOfSolverMatrices) THEN
      localError="The specified matrix index of "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
        & " is invalid. The index must be >= 1 and <= "// &
        & TRIM(NumberToVString(solverMapping%numberOfSolverMatrices,"*",err,error))//"."      
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(ALLOCATED(solverMapping%solverMatricesToEquationsMaps)) &
      & CALL FlagError("The solver matrices to equations map is not allocated for the solver mapping.",err,error,*999)
    IF(.NOT.ASSOCIATED(solverMapping%solverMatricesToEquationsMaps(matrixIdx)%ptr)) THEN
      localError="The solver matrices to equations maps is not associated for matrix index "// &
        & TRIM(NumberToVString(matrixIdx,"*",err,error))//"."
    ENDIF
#endif    

    columnDOFSMapping=>solverMapping%solverMatricesToEquationsMaps(matrixIdx)%ptr%columnDOFSMapping

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(columnDOFSMapping)) THEN
      localError="The column DOFs mapping for the specified matrix index of "// &
        & TRIM(NumberToVString(matrixIdx,"*",err,error))//" is not associated."      
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
  
  !>Returns a dynamic matrix in an equations set mapped to a solver matrix in a solver mapping.
  SUBROUTINE SolverMapping_DynamicMatrixToSolverMatrixGet(solverMapping,solverMatrixIdx,equationsSetIdx, &
    & dynamicMatrixIdx,dynamicMatrix,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the dynamic matrix for
    INTEGER(INTG), INTENT(IN) :: solverMatrixIdx !<The solver matrix index in the solver mapping to get the dynamic matrix for
    INTEGER(INTG), INTENT(IN) :: equationsSetIdx !<The equations set index containing the dynamic matrices mapped to the solver matrix
    INTEGER(INTG), INTENT(IN) :: dynamicMatrixIdx !<The dynamic matrix index of the equations matrix to the solver matrix
    TYPE(EquationsMatrixType), POINTER :: dynamicMatrix !<On return the specified dynamic matrix. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(EquationsMatrixToSolverMatrixMapType), POINTER :: dynamicMatrixToSolverMatrixMap
    TYPE(EquationsMatricesToSolverMatrixMapType), POINTER :: equationsMatricesToSolverMatrixMap
    TYPE(EquationsSetToSolverMatricesMapType), POINTER :: equationsSetToSolverMatricesMap
#endif
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMapping_DynamicMatrixToSolverMatrixGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dynamicMatrix)) CALL FlagError("Dynamic matrix is already associated.",err,error,*998)
    NULLIFY(equationsSetToSolverMatricesMap)
    CALL SolverMapping_EquationsSetToSolverMatricesMapGet(solverMapping,equationsSetIdx,equationsSetToSolverMatricesMap, &
      & err,error,*999)
    NULLIFY(equationsMatricesToSolverMatrixMap)
    CALL SolverMappingESToSMSMap_EquationsMatricesToSolverMatrixMapGet(equationsSetToSolverMatricesMap,solverMatrixIdx, &
      & equationsMatricesToSolverMatrixMap,err,error,*999)
    NULLIFY(dynamicMatrixToSolverMatrixMap)
    CALL SolverMappingEMSToSMMap_DynamicMatrixToSolverMatrixMapGet(equationsMatricesToSolverMatrixMap, &
      & dynamicMatrixIdx,dynamicMatrixToSolverMatrixMap,err,error,*999)
#endif    

    dynamicMatrix=>solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
      & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
      & dynamicMatrixToSolverMatrixMaps(dynamicMatrixIdx)%ptr%equationsMatrix

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(dynamicMatrix)) THEN
      localError="The dynamic matrix is not associated for the dynamic matrix index of "// &
        & TRIM(NumberToVString(dynamicMatrixIdx,"*",err,error))// &
        & " of the dynamic matrix to solver matrix maps for the solver matrix index of "// &
        & TRIM(NumberToVString(solverMatrixIdx,"*",err,error))// &
        & " of the equations matrices to solver matrix maps for the equations set index "// &
        & TRIM(NumberToVString(equationsSetIdx,"*",err,error))// &
        & " of the equations set to solver matrices map of the solver mapping."
      CALL FlagError(localError,err,error,*999)    
    ENDIF
#endif    
      
    EXITS("SolverMapping_DynamicMatrixToSolverMatrixGet")
    RETURN
999 NULLIFY(dynamicMatrix)
998 ERRORSEXITS("SolverMapping_DynamicMatrixToSolverMatrixGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_DynamicMatrixToSolverMatrixGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of dynamic matrices in an equations set mapped to a solver matrix in a solver mapping.
  SUBROUTINE SolverMapping_EMSToSMMapNumberOfDynamicMatricesGet(solverMapping,solverMatrixIdx,equationsSetIdx, &
    & numberOfDynamicMatrices,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the number of dynamic matrices for
    INTEGER(INTG), INTENT(IN) :: solverMatrixIdx !<The solver matrix index in the solver mapping to get the number of dynamic matrices for
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
      & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%numberOfDynamicMatrices
      
    EXITS("SolverMapping_EMSToSMMapNumberOfDynamicMatricesGet")
    RETURN
999 ERRORSEXITS("SolverMapping_EMSToSMMapNumberOfDynamicMatricesGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_EMSToSMMapNumberOfDynamicMatricesGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of interface matrices in an interface condition mapped to a solver matrix in a solver mapping.
  SUBROUTINE SolverMapping_IMSToSMMapNumberOfInterfaceMatricesGet(solverMapping,solverMatrixIdx,interfaceConditionIdx, &
    & numberOfInterfaceMatrices,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the number of interface matrices for
    INTEGER(INTG), INTENT(IN) :: solverMatrixIdx !<The solver matrix index in the solver mapping to get the number of interface matrices for
    INTEGER(INTG), INTENT(IN) :: interfaceConditionIdx !<The interface condition index containing the interface matrices mapped to the solver matrix
    INTEGER(INTG), INTENT(OUT) :: numberOfInterfaceMatrices !<On exit, the number of interface matrices mapped to the solver matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(InterfaceMatricesToSolverMatrixMapType), POINTER :: interfaceMatricesToSolverMatrixMap
    TYPE(InterfaceConditionToSolverMatricesMapType), POINTER :: interfaceConditionToSolverMatricesMap
#endif
 
    ENTERS("SolverMapping_IMSToSMMapNumberOfInterfaceMatricesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    NULLIFY(interfaceConditionToSolverMatricesMap)
    CALL SolverMapping_InterfaceConditionToSolverMatricesMapGet(solverMapping,interfaceConditionIdx, &
      & interfaceConditionToSolverMatricesMap,err,error,*999)
    NULLIFY(interfaceMatricesToSolverMatrixMap)
    CALL SolverMappingESToSMSMap_InterfaceMatricesToSolverMatrixMapGet(interfaceConditionToSolverMatricesMap,solverMatrixIdx, &
      & interfaceMatricesToSolverMatrixMap,err,error,*999)
#endif    
    
    numberOfInterfaceMatrices=solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
      & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%numberOfInterfaceMatrices
      
    EXITS("SolverMapping_IMSToSMMapNumberOfInterfaceMatricesGet")
    RETURN
999 ERRORSEXITS("SolverMapping_IMSToSMMapNumberOfInterfaceMatricesGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_IMSToSMMapNumberOfInterfaceMatricesGet
  
  !
  !================================================================================================================================
  !
  
  !!>Returns the number of Jacobian matrices in an interface condition mapped to a solver matrix in a solver mapping.
  !SUBROUTINE SolverMapping_IMSToSMMapNumberOfJacobianMatricesGet(solverMapping,solverMatrixIdx,interfaceConditionIdx, &
  !  & numberOfJacobianMatrices,err,error,*)

  !  !Argument variables
  !  TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the number of Jacobian matrices for
  !  INTEGER(INTG), INTENT(IN) :: solverMatrixIdx !<The solver matrix index in the solver mapping to get the number of Jacobian matrices for
  !  INTEGER(INTG), INTENT(IN) :: interfaceConditionIdx !<The interface condition index containing the Jacobian matrices mapped to the solver matrix
  !  INTEGER(INTG), INTENT(OUT) :: numberOfJacobianMatrices !<On exit, the number of Jacobian matrices mapped to the solver matrix
  !  INTEGER(INTG), INTENT(OUT) :: err !<The error code
  !  TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
!#ifdef WITH_PRECHECKS
  !  TYPE(InterfaceMatricesToSolverMatrixMapType), POINTER :: interfaceMatricesToSolverMatrixMap
  !  TYPE(InterfaceConditionToSolverMatricesMapType), POINTER :: interfaceConditionToSolverMatricesMap
!#endif
 
  !  ENTERS("SolverMapping_IMSToSMMapNumberOfJacobianMatricesGet",err,error,*999)

!#ifdef WITH_PRECHECKS    
  !  NULLIFY(interfaceConditionToSolverMatricesMap)
  !  CALL SolverMapping_InterfaceConditionToSolverMatricesMapGet(solverMapping,interfaceConditionidx, &
  !    & interfaceConditionToSolverMatricesMap,err,error,*999)
  !  NULLIFY(interfaceMatricesToSolverMatrixMap)
  !  CALL SolverMappingICToSMSMap_InterfaceMatricesToSolverMatrixMapGet(interfaceConditionToSolverMatricesMap,solverMatrixIdx, &
  !    & interfaceMatricesToSolverMatrixMap,err,error,*999)
!#endif    
    
  !  numberOfJacobianMatrices=solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
  !    & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%numberOfEquationsJacobians
      
  !  EXITS("SolverMapping_IMSToSMMapNumberOfJacobianMatricesGet")
  !  RETURN
!999 ERRORSEXITS("SolverMapping_IMSToSMMapNumberOfJacobianMatricesGet",err,error)
  !  RETURN 1
    
  !END SUBROUTINE SolverMapping_IMSToSMMapNumberOfJacobianMatricesGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of linear matrices in an equations set mapped to a solver matrix in a solver mapping.
  SUBROUTINE SolverMapping_EMSToSMMapNumberOfLinearMatricesGet(solverMapping,solverMatrixIdx,equationsSetIdx, &
    & numberOfLinearMatrices,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the number of linear matrices for
    INTEGER(INTG), INTENT(IN) :: solverMatrixIdx !<The solver matrix index in the solver mapping to get the number of linear matrices for
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
      & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%numberOfLinearMatrices
      
    EXITS("SolverMapping_EMSToSMMapNumberOfLinearMatricesGet")
    RETURN
999 ERRORSEXITS("SolverMapping_EMSToSMMapNumberOfLinearMatricesGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_EMSToSMMapNumberOfLinearMatricesGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of variables in an equations set mapped to a solver matrix in a solver mapping.
  SUBROUTINE SolverMapping_EMSToSMMapNumberOfVariablesGet(solverMapping,solverMatrixIdx,equationsSetIdx,numberOfVariables, &
    & err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the number of variables for
    INTEGER(INTG), INTENT(IN) :: solverMatrixIdx !<The solver matrix index in the solver mapping to get the number of variables for
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
  SUBROUTINE SolverMapping_EMSToSMMapVariableGet(solverMapping,solverMatrixIdx,equationsSetIdx,variableIdx,fieldVariable, &
    & err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the variable for
    INTEGER(INTG), INTENT(IN) :: solverMatrixIdx !<The solver matrix index in the solver mapping to get the variable for
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
        & " is invalid for the solver matrix index of "//TRIM(NumberToVString(solverMatrixIdx,"*",err,error))// &
        & " of the equations matrices to solver matrix maps for the equations set index "// &
        & TRIM(NumberToVString(equationsSetIdx,"*",err,error))// &
        & " of the equations set to solver matrices map of the solver mapping. The variable index should be >= 1 and <= "// &
        & TRIM(NumberToVString(equationsMatricesToSolverMatrixMap%numberOfVariables,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(equationsMatricesToSolverMatrixMap%variables)) THEN
      localError="The variables array is not allocated for the solver matrix index of "// &
        & TRIM(NumberToVString(solverMatrixIdx,"*",err,error))// &
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
        & TRIM(NumberToVString(solverMatrixIdx,"*",err,error))// &
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
  SUBROUTINE SolverMapping_EMSToSMMapVariableDOFToSolverDOFsMapGet(solverMapping,solverMatrixIdx,equationsSetIdx,variableIdx, &
    & varDOFToSolverDOFsMap,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the variable to solver col map for
    INTEGER(INTG), INTENT(IN) :: solverMatrixIdx !<The solver matrix index in the solver mapping to get the variable to solver col map for
    INTEGER(INTG), INTENT(IN) :: equationsSetIdx !<The equations set index containing the variables mapped to the solver matrix
    INTEGER(INTG), INTENT(IN) :: variableIdx !<The variable index in the equations set mapped to the solver matrix
    TYPE(VariableDOFToSolverDOFsMapType), POINTER :: varDOFToSolverDOFsMap !<On exit, a pointer to the specified variable to solver col map for the solver matrix. Must not be associated on entry
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
 
    ENTERS("SolverMapping_EMSToSMMapVariableDOFToSolverDOFsMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(varDOFToSolverDOFsMap)) CALL FlagError("Variable DOF to solver DOFs map is already associated.",err,error,*998)
    NULLIFY(equationsSetToSolverMatricesMap)
    CALL SolverMapping_EquationsSetToSolverMatricesMapGet(solverMapping,equationsSetIdx,equationsSetToSolverMatricesMap, &
      & err,error,*999)
    NULLIFY(equationsMatricesToSolverMatrixMap)
    CALL SolverMappingESToSMSMap_EquationsMatricesToSolverMatrixMapGet(equationsSetToSolverMatricesMap,solverMatrixIdx, &
      & equationsMatricesToSolverMatrixMap,err,error,*999)
    IF(variableIdx<1.OR.variableIdx>equationsMatricesToSolverMatrixMap%numberOfVariables) THEN
      localError="The specified variable index of "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
        & " is invalid for the solver matrix index of "//TRIM(NumberToVString(solverMatrixIdx,"*",err,error))// &
        & " of the equations matrices to solver matrix maps for the equations set index "// &
        & TRIM(NumberToVString(equationsSetIdx,"*",err,error))// &
        & " of the equations set to solver matrices map of the solver mapping. The variable index should be >= 1 and <= "// &
        & TRIM(NumberToVString(equationsMatricesToSolverMatrixMap%numberOfVariables,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(equationsMatricesToSolverMatrixMap%varDOFToSolverDOFsMaps)) THEN
      localError="The variable DOF to solver DOFs maps array is not allocated for the solver matrix index of "// &
        & TRIM(NumberToVString(solverMatrixIdx,"*",err,error))// &
        & " of the equations matrices to solver matrix maps sm the equations set index "// &
        & TRIM(NumberToVString(equationsSetIdx,"*",err,error))// &
        & " of the equations set to solver matrices map of the solver mapping."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
     
   varDOFToSolverDOFsMap=>solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
     & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
     & varDOFToSolverDOFsMaps(variableIdx)%ptr

#ifdef WITH_POSTCHECKS   
    IF(.NOT.ASSOCIATED(varDOFToSolverDOFsMap)) THEN
      localError="The variable DOF to solver DOFs map is not associated for the specified variable index of "// &
        & TRIM(NumberToVString(variableIdx,"*",err,error))//" for the solver matrix index of "// &
        & TRIM(NumberToVString(solverMatrixIdx,"*",err,error))// &
        & " of the equations matrices to solver matrix maps for the equations set index "// &
        & TRIM(NumberToVString(equationsSetIdx,"*",err,error))// &
        & " of the equations set to solver matrices map of the solver mapping."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("SolverMapping_EMSToSMMapVariableDOFToSolverDOFsMapGet")
    RETURN
999 NULLIFY(varDOFToSolverDOFsMap)
998 ERRORSEXITS("SolverMapping_EMSToSMMapVariableDOFToSolverDOFsMapGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_EMSToSMMapVariableDOFToSolverDOFsMapGet
  
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
  SUBROUTINE SolverMapping_EquationsMatrixToSolverMatrixMapGet(solverMapping,equationsSetIdx,equationsMatrixIdx, &
    & solverMatrixIdx,equationsMatrixToSolverMatrixMap,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the equations matrix to solver matrix map for
    INTEGER(INTG), INTENT(IN) :: equationsSetIdx !<The equations set index in the solver mapping to get the equations matrix to solver matrix map for
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
    CALL SolverMapping_EquationsSetToSolverMatricesMapGet(solverMapping,equationsSetIdx,equationsSetToSolverMatricesMap, &
      & err,error,*999)
    NULLIFY(equationsMatrixToSolverMatricesMap)
    CALL SolverMappingESToSMSMap_EquationsMatrixToSolverMatricesMapGet(equationsSetToSolverMatricesMap,equationsMatrixIdx, &
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
    IF(.NOT.ALLOCATED(equationsMatrixToSolverMatricesMap%equationsMatrixToSolverMatrixMaps)) THEN
      localError="The equations matrix to solver matrix maps is not allocated for the equations matrix index of "// &
        & TRIM(NumberToVString(equationsMatrixIdx,"*",err,error))//&
        & " of the equations matrix to solver matrices maps for the equations set index of "// &
        & TRIM(NumberToVString(equationsSetIdx,"*",err,error))//" of the equations set to solver matrices maps "// &
        & "of the solver mapping."
    ENDIF
#endif    
    
    equationsMatrixToSolverMatrixMap=>solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
      & equationsMatrixToSolverMatricesMaps(equationsMatrixIdx)%ptr% &
      & equationsMatrixToSolverMatrixMaps(solverMatrixIdx)%ptr
    
#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(equationsMatrixToSolverMatrixMap)) THEN
      localError="The equations matrix to solver matrix map is not associated for the solver matrix index of "// &
        & TRIM(NumberToVString(solverMatrixIdx,"*",err,error))// &
        & " of the equations matrix to solver matrix maps of the equations matrix index of "// &
        & TRIM(NumberToVString(equationsMatrixIdx,"*",err,error))// &
        & " of the equations matrix to solver matrices maps for the equations set index of "// &
        & TRIM(NumberToVString(equationsSetIdx,"*",err,error))// &
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
    IF(equationsSetIdx.OR.equationsSetIdx>solverMapping%numberOfEquationsSets) THEN
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
  SUBROUTINE SolverMapping_IMSToSMMapLagrangeVarDOFToSolverDOFsMapGet(solverMapping,solverMatrixIdx,interfaceConditionIdx, &
    & lagrangeVarDOFToSolverDOFsMap,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the Lagrange variable to solver col map for
    INTEGER(INTG), INTENT(IN) :: solverMatrixIdx !<The solver matrix index in the solver mapping to get the Lagrange variable to solver col map for
    INTEGER(INTG), INTENT(IN) :: interfaceConditionIdx !<The interface condition index containing the Lagrange variable mapped to the solver matrix
    TYPE(VariableDOFToSolverDOFsMapType), POINTER :: lagrangeVarDOFToSolverDOFsMap !<On exit, a pointer to the specified Lagrange variable DOF to solver DOFs map for the solver matrix. Must not be associated on entry
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
 
    ENTERS("SolverMapping_IMSToSMMapLagrangeVarDOFToSolverDOFsMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(lagrangeVarDOFToSolverDOFsMap)) &
      & CALL FlagError("Lagrange variable DOF to solver DOFs map is already associated.",err,error,*998)
    NULLIFY(interfaceConditionToSolverMatricesMap)
    CALL SolverMapping_InterfaceConditionToSolverMatricesMapGet(solverMapping,interfaceConditionIdx, &
      & interfaceConditionToSolverMatricesMap,err,error,*999)
    NULLIFY(interfaceMatricesToSolverMatrixMap)
    CALL SolverMappingICToSMSMap_InterfaceMatricesToSolverMatrixMapGet(interfaceConditionToSolverMatricesMap,solverMatrixIdx, &
      & interfaceMatricesToSolverMatrixMap,err,error,*999)
#endif    
     
   lagrangeVarDOFToSolverDOFsMap=>solverMapping%InterfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
     & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
     & lagrangeVarDOFToSolverDOFsMap

#ifdef WITH_POSTCHECKS   
    IF(.NOT.ASSOCIATED(lagrangeVarDOFToSolverDOFsMap)) THEN
      localError="The Lagrange variable DOF to solver DOFs map is not associated for the solver matrix index of "// &
        & TRIM(NumberToVString(solverMatrixIdx,"*",err,error))// &
        & " of the interface matrices to solver matrix maps for the interface condition index "// &
        & TRIM(NumberToVString(interfaceConditionIdx,"*",err,error))// &
        & " of the interface condition to solver matrices map of the solver mapping."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("SolverMapping_IMSToSMMapLagrangeVarDOFToSolverDOFsMapGet")
    RETURN
999 NULLIFY(lagrangeVarDOFToSolverDOFsMap)
998 ERRORS("SolverMapping_IMSToSMMapLagrangeVarDOFToSolverDOFsMapGet",err,error)
    EXITS("SolverMapping_IMSToSMMapLagrangeVarDOFToSolverDOFsMapGet")
    RETURN 1
    
  END SUBROUTINE SolverMapping_IMSToSMMapLagrangeVarDOFToSolverDOFsMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the dependent variable to solver col map in an equations set mapped to a solver matrix in a solver mapping.
  SUBROUTINE SolverMapping_IMSToSMMapDependentVarDOFToSolverDOFsMapGet(solverMapping,solverMatrixIdx,interfaceConditionIdx, &
    & variableIdx,dependentVarDOFToSolverDOFsMap,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the dependent variable to solver col map for
    INTEGER(INTG), INTENT(IN) :: solverMatrixIdx !<The solver matrix index in the solver mapping to get the dependent variable to solver col map for
    INTEGER(INTG), INTENT(IN) :: interfaceConditionIdx !<The interface condition index containing the dependent variable mapped to the solver matrix
    INTEGER(INTG), INTENT(IN) :: variableIdx !<The variable index of the dependent variable mapped to the solver matrix
    TYPE(VariableDOFToSolverDOFsMapType), POINTER :: dependentVarDOFToSolverDOFsMap !<On exit, a pointer to the specified dependent variable DOF to solver DOFs map for the solver matrix. Must not be associated on entry
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
 
    ENTERS("SolverMapping_IMSToSMMapDependentVarDOFToSolverDOFsMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dependentVarDOFToSolverDOFsMap)) &
      & CALL FlagError("Dependent variable DOF to solver DOFs map is already associated.",err,error,*998)
    NULLIFY(interfaceConditionToSolverMatricesMap)
    CALL SolverMapping_InterfaceConditionToSolverMatricesMapGet(solverMapping,interfaceConditionIdx, &
      & interfaceConditionToSolverMatricesMap,err,error,*999)
    NULLIFY(interfaceMatricesToSolverMatrixMap)
    CALL SolverMappingICToSMSMap_InterfaceMatricesToSolverMatrixMapGet(interfaceConditionToSolverMatricesMap,solverMatrixIdx, &
      & interfaceMatricesToSolverMatrixMap,err,error,*999)
    IF(variableIdx<1.OR.variableIdx>interfaceMatricesToSolverMatrixMap%numberOfDependentVariables) THEN
      localError="The specified dependent variable index of "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
        & " is invalid for the solver matrix index of "//TRIM(NumberToVString(solverMatrixIdx,"*",err,error))// &
        & " of the interface matrices to solver matrix maps for the interface condition index "// &
        & TRIM(NumberToVString(interfaceConditionIdx,"*",err,error))//" of the interface condition to solver matrices map of "// &
        & "the solver mapping. The dependent variable index should be >= 1 and <= "// &
        & TRIM(NumberToVString(interfaceMatricesToSolverMatrixMap%numberOfDependentVariables,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(interfaceMatricesToSolverMatrixMap%dependentVarDOFToSolverDOFsMaps)) THEN
      localError="The dependent variable DOF to solver DOFs maps array is not allocated for the solver matrix index of "// &
        & TRIM(NumberToVString(solverMatrixIdx,"*",err,error))// &
        & " of the interface matrices to solver matrix maps of the interface condition index "// &
        & TRIM(NumberToVString(interfaceConditionIdx,"*",err,error))// &
        & " of the interface condition to solver matrices map of the solver mapping."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
     
   dependentVarDOFToSolverDOFsMap=>solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
     & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
     & dependentVarDOFToSolverDOFsMaps(variableIdx)%ptr

#ifdef WITH_POSTCHECKS   
    IF(.NOT.ASSOCIATED(dependentVarDOFToSolverDOFsMap)) THEN
      localError="The dependent variable DOF to solver DOFs map is not associated for the specified variable index of "// &
        & TRIM(NumberToVString(variableIdx,"*",err,error))//" for the solver matrix index of "// &
        & TRIM(NumberToVString(solverMatrixIdx,"*",err,error))// &
        & " of the interface matrices to solver matrix maps for the interface condition index "// &
        & TRIM(NumberToVString(interfaceConditionIdx,"*",err,error))// &
        & " of the interface condition to solver matrices map of the solver mapping."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("SolverMapping_IMSToSMMapDependentVarDOFToSolverDOFsMapGet")
    RETURN
999 NULLIFY(dependentVarDOFToSolverDOFsMap)
998 ERRORS("SolverMapping_IMSToSMMapDependentVarDOFToSolverDOFsMapGet",err,error)
    EXITS("SolverMapping_IMSToSMMapDependentVarDOFToSolverDOFsMapGet")
    RETURN 1
    
  END SUBROUTINE SolverMapping_IMSToSMMapDependentVarDOFToSolverDOFsMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the interface col to solver cols map for solver mapping.xs
  SUBROUTINE SolverMapping_InterfaceColToSolverColsMapGet(solverMapping,interfaceConditionIdx,solverMatrixIdx, &
    & interfaceColToSolverColsMap,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the interface col to solver cols map for
    INTEGER(INTG), INTENT(IN) :: interfaceConditionIdx !<The interface condition index in the solver mapping to get the interface col to solver cols map for
    INTEGER(INTG), INTENT(IN) :: solverMatrixIdx !<The solver matrix index in the solver mapping to get the interface col to solver cols map for
    TYPE(MatrixRowColCouplingType), POINTER :: interfaceColToSolverColsMap(:) !<On exit, a pointer to the specified interface col to solver cols map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(InterfaceConditionToSolverMatricesMapType), POINTER :: interfaceConditionToSolverMatricesMap
    TYPE(InterfaceMatricesToSolverMatrixMapType), POINTER :: interfaceMatricesToSolverMatrixMap
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
      & interfaceConditionToSolverMatricesMap,err,error,*999)
    NULLIFY(interfaceMatricesToSolverMatrixMap)
    CALL SolverMappingICToSMSMap_InterfaceMatricesToSolverMatrixMapGet(interfaceConditionToSolverMatricesMap,solverMatrixIdx, &
      & interfaceMatricesToSolverMatrixMap,err,error,*999)
#endif

    interfaceColToSolverColsMap=>solverMapping%interfaceConditionToSolverMatricesMaps(interfaceConditionIdx)%ptr% &
      & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr%interfaceColToSolverColsMap

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
#ifdef WITH_CHECKS    
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
  SUBROUTINE SolverMapping_InterfaceEquationsMatrixToSolverMatrixGet(solverMapping,solverMatrixIdx,interfaceConditionIdx, &
    & interfaceMatrixIdx,interfaceEquationsMatrix,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the interface equations matrix for
    INTEGER(INTG), INTENT(IN) :: solverMatrixIdx !<The solver matrix index in the solver mapping to get the interface equations matrix for
    INTEGER(INTG), INTENT(IN) :: interfaceConditionIdx !<The interface condition index containing the interface matrices mapped to the solver matrix
    INTEGER(INTG), INTENT(IN) :: interfaceMatrixIdx !<The interface matrix index of the interface equations matrix to the solver matrix
    TYPE(InterfaceMatrixType), POINTER :: interfaceEquationsMatrix !<On return the specified interface equations matrix. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(InterfaceConditionToSolverMatricesMapType), POINTER :: interfaceConditionToSolverMatricesMap
    TYPE(InterfaceMatricesToSolverMatrixMapType), POINTER :: interfaceMatricesToSolverMatrixMap
    TYPE(InterfaceMatrixToSolverMatrixMapType), POINTER :: interfaceMatrixToSolverMatrixMap
#endif
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMapping_InterfaceEquationsMatrixToSolverMatrixGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(interfaceEquationsMatrix)) CALL FlagError("Interface equations matrix is already associated.",err,error,*998)
    NULLIFY(interfaceConditionToSolverMatricesMap)
    CALL SolverMapping_InterfaceConditionToSolverMatricesMapGet(solverMapping,interfaceConditionIdx, &
      & interfaceConditionToSolverMatricesMap,err,error,*999)
    NULLIFY(interfaceMatricesToSolverMatrixMap)
    CALL SolverMappingICToSMSMap_InterfaceMatricesToSolverMatrixMapGet(interfaceConditionToSolverMatricesMap,solverMatrixIdx, &
      & interfaceMatricesToSolverMatrixMap,err,error,*999)
    NULLIFY(interfaceMatrixToSolverMatrixMap)
    CALL SolverMappingIMSToSMMap_InterfaceMatrixToSolverMatrixMapGet(interfaceMatricesToSolverMatrixMap,interfaceMatrixIdx, &
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
        & TRIM(NumberToVString(solverMatrixIdx,"*",err,error))// &
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
    TYPE(InterfaceMatrixToSolverMatrixMapType), POINTER :: interfaceMatrixToSolverMatrixMap !<On exit, a pointer to the specified interface matrix to solver matrix map. Must not be associated on entry.
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
    CALL SolverMappingICToSMSMap_InterfaceMatrixToSolverMatricesMapGet(interfaceConditionToSolverMatricesMap,interfaceMatrixIdx, &
      & interfaceMatrixToSolverMatricesMap,err,error,*999)
    IF(solverMatrixIdx<1.OR.solverMatrixIdx>interfaceMatrixToSolverMatricesMap%numberOfSolverMatrices) THEN
      localError="The specified solver matrix index of "//TRIM(NumberToVString(solverMatrixIdx,"*",err,error))//&
        & " is invalid for the interface matrix index of "//TRIM(NumberToVString(interfaceMatrixIdx,"*",err,error))//&
        & " of the interface matrix to solver matrices maps for the interface condition index of "// &
        & TRIM(NumberToVString(interfaceConditionIdx,"*",err,error))//" of the interface condition to solver matrices maps "// &
        & "of the solver mapping. The solver matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(interfaceMatrixToSolverMatricesMap%numberOfSolverMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)      
    ENDIF
    IF(.NOT.ALLOCATED(interfaceMatrixToSolverMatricesMap%interfaceMatrixToSolverMatrixMaps)) THEN
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
  SUBROUTINE SolverMapping_InterfaceRowToSolverRowsMapGet(solverMapping,interfaceConditionIdx,interfaceMatrixIdx, &
    & interfaceRowToSolverRowsMap,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the interface row to solver rows map for
    INTEGER(INTG), INTENT(IN) :: interfaceConditionIdx !<The interface condition index in the solver mapping to get the interface row to solver rows map for
    INTEGER(INTG), INTENT(IN) :: interfaceMatrixIdx !<The solver matrix index in the solver mapping to get the interface row to solver rows map for
    TYPE(MatrixRowColCouplingType), POINTER :: interfaceRowToSolverRowsMap(:) !<On exit, a pointer to the specified interface row to solver rows map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(InterfaceConditionToSolverMatricesMapType), POINTER :: interfaceConditionToSolverMatricesMap
    TYPE(InterfaceMatrixToSolverMatricesMapType), POINTER :: interfaceMatrixToSolverMatricesMap
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
      & interfaceMatrixToSolverMatricesMaps(interfaceMatrixIdx)%ptr%interfaceRowToSolverRowsMap

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
  SUBROUTINE SolverMapping_JacobianEquationsMatrixToSolverMatrixGet(solverMapping,solverMatrixIdx,equationsSetIdx, &
    & jacobianMatrixIdx,jacobianEquationsMatrix,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the Jacobian equations matrix for
    INTEGER(INTG), INTENT(IN) :: solverMatrixIdx !<The solver matrix index in the solver mapping to get the Jacobian equations matrix for
    INTEGER(INTG), INTENT(IN) :: equationsSetIdx !<The equations set index containing the Jacobian matrices mapped to the solver matrix
    INTEGER(INTG), INTENT(IN) :: jacobianMatrixIdx !<The Jacobian matrix index of the Jacobian matrix to the solver matrix
    TYPE(JacobianMatrixType), POINTER :: jacobianEquationsMatrix !<On return the specified Jacobian equations matrix. Must not be associated on entry.
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
    IF(ASSOCIATED(jacobianEquationsMatrix)) CALL FlagError("Jacobian equations matrix is already associated.",err,error,*998)
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
    
    jacobianEquationsMatrix=>solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
      & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
      & jacobianMatrixToSolverMatrixMaps(jacobianMatrixIdx)%ptr%jacobianMatrix

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(jacobianEquationsMatrix)) THEN
      localError="The Jacobian equations matrix is not associated for the Jacobian matrix index of "// &
        & TRIM(NumberToVString(jacobianMatrixIdx,"*",err,error))// &
        & " of the Jacobian matrix to solver matrix maps for the solver matrix index of "// &
        & TRIM(NumberToVString(solverMatrixIdx,"*",err,error))// &
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
  SUBROUTINE SolverMapping_JacobianMatrixToSolverMatrixMapGet(solverMapping,equationsSetIdx,jacobianMatrixIdx, &
    & jacobianMatrixToSolverMatrixMap,err,error,*)
    
    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the Jacobian matrix to solver matrix map for
    INTEGER(INTG), INTENT(IN) :: equationsSetIdx !<The equations set index in the solver mapping to get the Jacobian matrix to solver matrix map for
    INTEGER(INTG), INTENT(IN) :: jacobianMatrixIdx !<The jacobian matrix index for the Jacobian matrix to solver matrix map
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
    IF(jacobianMatrixIdx<1.OR.jacobianMatrixIdx>SIZE(equationsSetToSolverMatricesMap%jacobianMatrixToSolverMatrixMaps,1)) THEN
      localError="The specified Jacobian matrix index of "//TRIM(NumberToVString(jacobianMatrixIdx,"*",err,error))// &
        & " is invalid for the equations set index of "//TRIM(NumberToVString(equationsSetIdx,"*",err,error))// &
        & " of the equations set to solver matrices map of the solver mapping. The Jacobian matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(equationsSetToSolverMatricesMap%jacobianMatrixToSolverMatrixMaps,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(equationsSetToSolverMatricesMap%jacobianMatrixToSolverMatrixMaps)) THEN
      localError="The Jacobian matrix to solver matrix maps is not allocated for the equations set index of "// &
        & TRIM(NumberToVString(equationsSetIdx,"*",err,error))// &
        & " of the equations set to solver matrices map of the solver mapping."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
   
    jacobianMatrixToSolverMatrixMap=>solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
      & jacobianMatrixToSolverMatrixMaps(jacobianMatrixIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(jacobianMatrixToSolverMatrixMap)) THEN
      localError="The Jacobian matrix to solver matrix map is not associated for the Jacobian matrix index of "// &
        & TRIM(NumberToVString(jacobianMatrixIdx,"*",err,error))// &
        & " of the Jacobian matrix to solver matrix maps of the equations set index of "// &
        & TRIM(NumberToVString(equationsSetIdx,"*",err,error))// &
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
  
  !>Returns a linear matrix in an equations set mapped to a solver matrix in a solver mapping.
  SUBROUTINE SolverMapping_LinearMatrixToSolverMatrixGet(solverMapping,solverMatrixIdx,equationsSetIdx, &
    & linearMatrixIdx,linearMatrix,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the linear matrix for
    INTEGER(INTG), INTENT(IN) :: solverMatrixIdx !<The solver matrix index in the solver mapping to get the linear matrix for
    INTEGER(INTG), INTENT(IN) :: equationsSetIdx !<The equations set index containing the linear matrices mapped to the solver matrix
    INTEGER(INTG), INTENT(IN) :: linearMatrixIdx !<The linear matrix index of the linear matrix to the solver matrix
    TYPE(EquationsMatrixType), POINTER :: linearMatrix !<On return the specified linear matrix. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(EquationsMatrixToSolverMatrixMapType), POINTER :: linearMatrixToSolverMatrixMap
    TYPE(EquationsMatricesToSolverMatrixMapType), POINTER :: equationsMatricesToSolverMatrixMap
    TYPE(EquationsSetToSolverMatricesMapType), POINTER :: equationsSetToSolverMatricesMap
#endif
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMapping_LinearMatrixToSolverMatrixGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(linearMatrix)) CALL FlagError("Linear equations matrix is already associated.",err,error,*998)
    NULLIFY(equationsSetToSolverMatricesMap)
    CALL SolverMapping_EquationsSetToSolverMatricesMapGet(solverMapping,equationsSetIdx,equationsSetToSolverMatricesMap, &
      & err,error,*999)
    NULLIFY(equationsMatricesToSolverMatrixMap)
    CALL SolverMappingESToSMSMap_EquationsMatricesToSolverMatrixMapGet(equationsSetToSolverMatricesMap,solverMatrixIdx, &
      & equationsMatricesToSolverMatrixMap,err,error,*999)
    NULLIFY(linearMatrixToSolverMatrixMap)
    CALL SolverMappingEMSToSMMap_LinearMatrixToSolverMatrixMapGet(equationsMatricesToSolverMatrixMap, &
      & linearMatrixIdx,linearMatrixToSolverMatrixMap,err,error,*999)
#endif    

    linearMatrix=>solverMapping%equationsSetToSolverMatricesMaps(equationsSetIdx)%ptr% &
      & equationsMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr% &
      & linearMatrixToSolverMatrixMaps(linearMatrixIdx)%ptr%equationsMatrix

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(linearMatrix)) THEN
      localError="The linear matrix is not associated for the linear matrix index of "// &
        & TRIM(NumberToVString(linearMatrixIdx,"*",err,error))// &
        & " of the linear matrix to solver matrix maps for the solver matrix index of "// &
        & TRIM(NumberToVString(solverMatrixIdx,"*",err,error))// &
        & " of the equations matrices to solver matrix maps for the equations set index "// &
        & TRIM(NumberToVString(equationsSetIdx,"*",err,error))// &
        & " of the equations set to solver matrices map of the solver mapping."
      CALL FlagError(localError,err,error,*999)    
    ENDIF
#endif    
      
    EXITS("SolverMapping_LinearMatrixToSolverMatrixGet")
    RETURN
999 NULLIFY(linearMatrix)
998 ERRORSEXITS("SolverMapping_LinearMatrixToSolverMatrixGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_LinearMatrixToSolverMatrixGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of equations sets in a solver mapping.
  SUBROUTINE SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the number of equations sets for
    INTEGER(INTG), INTENT(OUT) :: numberOfEquationsSets !<On exit, the number of equations sets in the solver mapping.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMapping_NumberOfEquationsSetsGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solverMapping)) CALL FlagError("Solver mapping is not associated.",err,error,*999)
#endif

    numberOfEquationsSets=solverMapping%numberOfEquationsSets
      
    EXITS("SolverMapping_NumberOfEquationsSetsGet")
    RETURN
999 ERRORSEXITS("SolverMapping_NumberOfEquationsSetsGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_NumberOfEquationsSetsGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of interface conditions in a solver mapping.
  SUBROUTINE SolverMapping_NumberOfInterfaceConditionsGet(solverMapping,numberOfInterfaceConditions,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the number of interface conditions for
    INTEGER(INTG), INTENT(OUT) :: numberOfInterfaceConditions !<On exit, the number of interface conditions in the solver mapping.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMapping_NumberOfInterfaceConditionsGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solverMapping)) CALL FlagError("Solver mapping is not associated.",err,error,*999)
#endif

    numberOfInterfaceConditions=solverMapping%numberOfInterfaceConditions
      
    EXITS("SolverMapping_NumberOfInterfaceConditionsGet")
    RETURN
999 ERRORSEXITS("SolverMapping_NumberOfInterfaceConditionsGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_NumberOfInterfaceConditionsGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of solver matrices in a solver mapping.
  SUBROUTINE SolverMapping_NumberOfSolverMatricesGet(solverMapping,numberOfSolverMatrices,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the number of solver matrices for
    INTEGER(INTG), INTENT(OUT) :: numberOfSolverMatrices !<On exit, the number of solver matrices in the solver mapping.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMapping_NumberOfSolverMatricesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solverMapping)) CALL FlagError("Solver mapping is not associated.",err,error,*999)
#endif

    numberOfSolverMatrices=solverMapping%numberOfSolverMatrices
      
    EXITS("SolverMapping_NumberOfSolverMatricesGet")
    RETURN
999 ERRORSEXITS("SolverMapping_NumberOfSolverMatricesGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_NumberOfSolverMatricesGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of rows in a solver mapping.
  SUBROUTINE SolverMapping_NumberOfRowsGet(solverMapping,numberOfRows,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the number of rows for
    INTEGER(INTG), INTENT(OUT) :: numberOfRows !<On exit, the number of rows in the solver mapping.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMapping_NumberOfRowsGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solverMapping)) CALL FlagError("Solver mapping is not associated.",err,error,*999)
#endif

    numberOfRows=solverMapping%numberOfRows
      
    EXITS("SolverMapping_NumberOfRowsGet")
    RETURN
999 ERRORSEXITS("SolverMapping_NumberOfRowsGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_NumberOfRowsGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of global rows in a solver mapping.
  SUBROUTINE SolverMapping_NumberOfGlobalRowsGet(solverMapping,numberOfGlobalRows,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the number of global rows for
    INTEGER(INTG), INTENT(OUT) :: numberOfGlobalRows !<On exit, the number of global rows in the solver mapping.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMapping_NumberOfGlobalRowsGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solverMapping)) CALL FlagError("Solver mapping is not associated.",err,error,*999)
#endif

    numberOfGlobalRows=solverMapping%numberOfGlobalRows
      
    EXITS("SolverMapping_NumberOfGlobalRowsGet")
    RETURN
999 ERRORSEXITS("SolverMapping_NumberOfGlobalRowsGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_NumberOfGlobalRowsGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the RHS variables list for solver mapping.
  SUBROUTINE SolverMapping_RHSVariablesListGet(solverMapping,rhsVariablesList,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the RHS variables list for
    TYPE(SolverMappingVariablesType), POINTER :: rhsVariablesList !<On exit, a pointer to the RHS variables list for the solver mapping. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMapping_RHSVariablesListGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(rhsVariablesList)) CALL FlagError("RHS variables list is already associated.",err,error,*998)    
    IF(.NOT.ASSOCIATED(solverMapping)) CALL FlagError("Solver mapping is not associated.",err,error,*999)
#endif

    rhsVariablesList=>solverMapping%rhsVariablesList

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(rhsVariablesList)) &
      & CALL FlagError("The RHS variables list is not associated for the solver mapping.",err,error,*999)
#endif    
      
    EXITS("SolverMapping_RHSVariablesListGet")
    RETURN
999 NULLIFY(rhsVariablesList)
998 ERRORSEXITS("SolverMapping_RHSVariablesListGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_RHSVariablesListGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the row DOFs mapping for solver mapping.
  SUBROUTINE SolverMapping_RowDOFSMappingGet(solverMapping,rowDOFSMapping,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the row dofs mapping for
    TYPE(DomainMappingType), POINTER :: rowDOFSMapping !<On exit, a pointer to the specified row DOFS mapping. Must not be associated on entry.
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
    IF(.NOT.ALLOCATED(solverMapping%solverMatricesToEquationsMaps)) &
      & CALL FlagError("The solver matrices to equations maps is not allocated for the solver mapping.",err,error,*999)
#endif

    solverMatrixToEquationsMap=>solverMapping%solverMatricesToEquationsMaps(solverMatrixIdx)%ptr

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
  
  !>Returns a pointer to the solver row to equations map for a solver mapping.
  SUBROUTINE SolverMapping_SolverRowToEquationsMapGet(solverMapping,rowIdx,solverRowToEquationsMap,err,error,*)

    !Argument variables
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping to get the solver row to equations map for
    INTEGER(INTG), INTENT(IN) :: rowIdx !<The row index to get the solver row to equations map for
    TYPE(SolverRowToEquationsMapType), POINTER :: solverRowToEquationsMap !<On exit, a pointer to the specified solver row to equations map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif
 
    ENTERS("SolverMapping_SolverRowToEquationsMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(solverRowToEquationsMap)) &
      & CALL FlagError("The solver row to equations map is already associated.",err,error,*998)  
    IF(.NOT.ASSOCIATED(solverMapping)) CALL FlagError("Solver mapping is not associated.",err,error,*999)
    IF(rowIdx<1.OR.rowIdx>solverMapping%numberOfRows) THEN
      localError="The specified solver row index of "//TRIM(NumberToVString(rowIdx,"*",err,error))// &
        & " is invalid. The solver row index should be >= 1 and <= "// &
        & TRIM(NumberToVString(solverMapping%numberOfRows,"*",err,error))//"."
    ENDIF
    IF(.NOT.ALLOCATED(solverMapping%solverRowToEquationsMaps)) &
      & CALL FlagError("The solver row to equations maps is not allocated for the solver mapping.",err,error,*999)
#endif

    solverRowToEquationsMap=>solverMapping%solverRowToEquationsMaps(rowIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(solverRowToEquationsMap)) THEN
      localError="The solver row to equations map is not associated for row index "// &
        & TRIM(NumberToVString(rowIdx,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("SolverMapping_SolverRowToEquationsMapGet")
    RETURN
999 NULLIFY(solverRowToEquationsMap)
998 ERRORSEXITS("SolverMapping_SolverRowToEquationsMapGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMapping_SolverRowToEquationsMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to a equations variable list for solver mapping create values cache.
  SUBROUTINE SolverMappingCVC_EquationsVariableListGet(createValuesCache,solverMatrixIdx,equationsVariableList, &
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
 
    ENTERS("SolverMappingCVC_EquationsVariableListGet",err,error,*998)

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

    equationsVariableList=>createValuesCache%equationsVariableList(solverMatrixIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(equationsVariableList)) THEN
      localError="The equations variable list is not associated for solver matrix index "// &
        & TRIM(NumberToVString(solverMatrixIdx,"*",err,error))//" for the solver mapping create values cache."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("SolverMappingCVC_EquationsVariableListGet")
    RETURN
999 NULLIFY(equationsVariableList)
998 ERRORS("SolverMappingCVC_EquationsVariableListGet",err,error)
    EXITS("SolverMappingCVC_EquationsVariableListGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingCVC_EquationsVariableListGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to a interface variable list for solver mapping create values cache.
  SUBROUTINE SolverMappingCVC_InterfaceVariableListGet(createValuesCache,solverMatrixIdx,interfaceVariableList, &
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
 
    ENTERS("SolverMappingCVC_InterfaceVariableListGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(interfaceVariableList)) &
      & CALL FlagError("Interface variable list is already associated.",err,error,*998)  
    IF(.NOT.ASSOCIATED(createValuesCache)) CALL FlagError("Solver mapping create values cache is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(createValuesCache%interfaceVariableList)) &
      & CALL FlagError("The interface variable list is not associated for the solver mapping create values cache.",err,error,*999)
    IF(solverMatrixIdx<1.OR.solverMatrixIdx>SIZE(createValuesCache%interfaceVariableList,1)) THEN
      localError="The specified solver matrix index of "//TRIM(NumberToVString(solverMatrixIdx,"*",err,error))// &
        & " is invalid. The solver matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(createValuesCache%interfaceVariableList,1),"*",err,error))//"."
    ENDIF
#endif

    interfaceVariableList=>createValuesCache%interfaceVariableList(solverMatrixIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(interfaceVariableList)) THEN
      localError="The interface variable list is not associated for solver matrix index "// &
        & TRIM(NumberToVString(solverMatrixIdx,"*",err,error))//" for the solver mapping create values cache."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("SolverMappingCVC_InterfaceVariableListGet")
    RETURN
999 NULLIFY(interfaceVariableList)
998 ERRORS("SolverMappingCVC_InterfaceVariableListGet",err,error)
    EXITS("SolverMappingCVC_InterfaceVariableListGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingCVC_InterfaceVariableListGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to a RHS variable list for solver mapping create values cache.
  SUBROUTINE SolverMappingCVC_RHSVariableListGet(createValuesCache,rhsVariableList,err,error,*)

    !Argument variables
    TYPE(SolverMappingCreateValuesCacheType), POINTER :: createValuesCache !<A pointer to the solver mapping create values cache to get the RHS variable list for
    TYPE(ListType), POINTER :: rhsVariableList !<On exit, a pointer to the specified RHS variable list. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingCVC_RHSVariableListGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(rhsVariableList)) &
      & CALL FlagError("RHS variable list is already associated.",err,error,*998)  
    IF(.NOT.ASSOCIATED(createValuesCache)) CALL FlagError("Solver mapping create values cache is not associated.",err,error,*999)
#endif

    rhsVariableList=>createValuesCache%equationsRHSVariableList

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(rhsVariableList)) &
      & CALL FlagError("The RHS variable list is not associated for the solver mapping create values cache.",err,error,*999)
#endif    
      
    EXITS("SolverMappingCVC_RHSVariableListGet")
    RETURN
999 NULLIFY(rhsVariableList)
998 ERRORS("SolverMappingCVC_RHSVariableListGet",err,error)
    EXITS("SolverMappingCVC_RHSVariableListGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingCVC_RHSVariableListGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the dynamic matrix to solver matrix map for an equations matrices to solver matrix map.
  SUBROUTINE SolverMappingEMSToSMMap_DynamicMatrixToSolverMatrixMapGet(equationsMatricesToSolverMatrixMap, &
    & dynamicMatrixIdx,dynamicMatrixToSolverMatrixMap,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesToSolverMatrixMapType), POINTER :: equationsMatricesToSolverMatrixMap !<A pointer to the equations matrices to solver matrix map to get the dynamic matrix to solver matrix map for
    INTEGER(INTG), INTENT(IN) :: dynamicMatrixIdx !<The dynamic equations matrix index to get the dynamic matrix to solver matrix for
    TYPE(EquationsMatrixToSolverMatrixMapType), POINTER :: dynamicMatrixToSolverMatrixMap !<On exit, a pointer to the specified dynamic matrix to solver matrix map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMappingEMSToSMMap_DynamicMatrixToSolverMatrixMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dynamicMatrixToSolverMatrixMap)) &
      & CALL FlagError("Dynamic matrix to solver matrix map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsMatricesToSolverMatrixMap)) &
      & CALL FlagError("Equations matrices to solver matrix map is not associated.",err,error,*999)
    IF(dynamicMatrixIdx<1.OR.dynamicMatrixIdx>equationsMatricesToSolverMatrixMap%numberOfDynamicMatrices) THEN
      localError="The specified dynamic matrix index of "// &
        & TRIM(NumberToVString(dynamicMatrixIdx,"*",err,error))//" is invalid for the equations matrices to solver "// &
        & "matrix map. The dynamic equations matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(equationsMatricesToSolverMatrixMap%numberOfDynamicMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(equationsMatricesToSolverMatrixMap%dynamicMatrixToSolverMatrixMaps)) &
      & CALL FlagError("The dynamic equations to solver matrix maps array is not allocated for the equations matrices to "// &
      & "solver matrix map.",err,error,*999)
#endif    

    dynamicMatrixToSolverMatrixMap=>equationsMatricesToSolverMatrixMap%dynamicMatrixToSolverMatrixMaps(dynamicMatrixIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(dynamicMatrixToSolverMatrixMap)) THEN
      localError="The dynamic matrix to solver map is not associated for dynamic matrix index "// &
        & TRIM(NumberToVString(dynamicMatrixIdx,"*",err,error))//" of the equations matrices to solver matrix map."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("SolverMappingEMSToSMMap_DynamicMatrixToSolverMatrixMapGet")
    RETURN
999 NULLIFY(dynamicMatrixToSolverMatrixMap)
998 ERRORS("SolverMappingEMSToSMMap_DynamicMatrixToSolverMatrixMapGet",err,error)
    EXITS("SolverMappingEMSToSMMap_DynamicMatrixToSolverMatrixMapGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingEMSToSMMap_DynamicMatrixToSolverMatrixMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the Jacobian matrix to solver matrix map for an equations matrices to solver matrix map.
  SUBROUTINE SolverMappingEMSToSMMap_JacobianMatrixToSolverMatrixMapGet(equationsMatricesToSolverMatrixMap, &
    & jacobianMatrixIdx,jacobianMatrixToSolverMatrixMap,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesToSolverMatrixMapType), POINTER :: equationsMatricesToSolverMatrixMap !<A pointer to the equations matrices to solver matrix map to get the Jacobian matrix to solver matrix map for
    INTEGER(INTG), INTENT(IN) :: jacobianMatrixIdx !<The Jacobian matrix index to get the Jacobian matrix to solver matrix for
    TYPE(JacobianMatrixToSolverMatrixMapType), POINTER :: jacobianMatrixToSolverMatrixMap !<On exit, a pointer to the specified Jacobian matrix to solver matrix map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMappingEMSToSMMap_JacobianMatrixToSolverMatrixGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(jacobianMatrixToSolverMatrixMap)) &
      & CALL FlagError("Jacobian matrix to solver matrix map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsMatricesToSolverMatrixMap)) &
      & CALL FlagError("Equations matrices to solver matrix map is not associated.",err,error,*999)
    IF(jacobianMatrixIdx<1.OR.jacobianMatrixIdx>equationsMatricesToSolverMatrixMap%numberOfJacobianMatrices) THEN
      localError="The specified Jacobian matrix index of "// &
        & TRIM(NumberToVString(jacobianMatrixIdx,"*",err,error))//" is invalid for the equations matrices to solver "// &
        & "matrix map. The Jacobian matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(equationsMatricesToSolverMatrixMap%numberOfJacobianMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(equationsMatricesToSolverMatrixMap%jacobianMatrixToSolverMatrixMaps)) &
      & CALL FlagError("The Jacobian equations to solver matrix maps array is not allocated for the equations matrices to "// &
      & "solver matrix map.",err,error,*999)
#endif    

    jacobianMatrixToSolverMatrixMap=>equationsMatricesToSolverMatrixMap%jacobianMatrixToSolverMatrixMaps(jacobianMatrixIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(jacobianMatrixToSolverMatrixMap)) THEN
      localError="The Jacobian matrix to solver map is not associated for Jacobian equations matrix index "// &
        & TRIM(NumberToVString(jacobianMatrixIdx,"*",err,error))//" of the equations matrices to solver matrix map."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("SolverMappingEMSToSMMap_JacobianMatrixToSolverMatrixMapGet")
    RETURN
999 NULLIFY(jacobianMatrixToSolverMatrixMap)
998 ERRORS("SolverMappingEMSToSMMap_JacobianMatrixToSolverMatrixMapGet",err,error)
    EXITS("SolverMappingEMSToSMMap_JacobianMatrixToSolverMatrixMapGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingEMSToSMMap_JacobianMatrixToSolverMatrixMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the linear matrix to solver matrix map for an equations matrices to solver matrix map.
  SUBROUTINE SolverMappingEMSToSMMap_LinearMatrixToSolverMatrixMapGet(equationsMatricesToSolverMatrixMap, &
    & linearMatrixIdx,linearMatrixToSolverMatrixMap,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesToSolverMatrixMapType), POINTER :: equationsMatricesToSolverMatrixMap !<A pointer to the equations matrices to solver matrix map to get the linear matrix to solver matrix map for
    INTEGER(INTG), INTENT(IN) :: linearMatrixIdx !<The linear matrix index to get the linear matrix to solver matrix for
    TYPE(EquationsMatrixToSolverMatrixMapType), POINTER :: linearMatrixToSolverMatrixMap !<On exit, a pointer to the specified linear matrix to solver matrix map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMappingEMSToSMMap_LinearMatrixToSolverMatrixGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(linearMatrixToSolverMatrixMap)) &
      & CALL FlagError("Linear matrix to solver matrix map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsMatricesToSolverMatrixMap)) &
      & CALL FlagError("Equations matrices to solver matrix map is not associated.",err,error,*999)
    IF(linearMatrixIdx<1.OR.linearMatrixIdx>equationsMatricesToSolverMatrixMap%numberOfLinearMatrices) THEN
      localError="The specified linear matrix index of "//TRIM(NumberToVString(linearMatrixIdx,"*",err,error))// &
        & " is invalid for the equations matrices to solver matrix map. "// &
        & "The linear equations matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(equationsMatricesToSolverMatrixMap%numberOfLinearMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(equationsMatricesToSolverMatrixMap%linearMatrixToSolverMatrixMaps)) &
      & CALL FlagError("The linear matrix to solver matrix maps array is not allocated for the equations matrices to "// &
      & "solver matrix map.",err,error,*999)
#endif    

    linearMatrixToSolverMatrixMap=>equationsMatricesToSolverMatrixMap%linearMatrixToSolverMatrixMaps(linearMatrixIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(linearMatrixToSolverMatrixMap)) THEN
      localError="The linear matrix to solver map is not associated for linear matrix index "// &
        & TRIM(NumberToVString(linearMatrixIdx,"*",err,error))//" of the equations matrices to solver matrix map."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("SolverMappingEMSToSMMap_LinearMatrixToSolverMatrixMapGet")
    RETURN
999 NULLIFY(linearMatrixToSolverMatrixMap)
998 ERRORS("SolverMappingEMSToSMMap_LinearMatrixToSolverMatrixMapGet",err,error)
    EXITS("SolverMappingEMSToSMMap_LinearMatrixToSolverMatrixMapGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingEMSToSMMap_LinearMatrixToSolverMatrixMapGet
  
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
    
    numberOfDynamicMatrices=equationsMatricesToSolverMatrixMap%numberOfDynamicMatrices
      
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
    
    numberOfJacobianMatrices=equationsMatricesToSolverMatrixMap%numberOfJacobianMatrices
      
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
    
    numberOfLinearMatrices=equationsMatricesToSolverMatrixMap%numberOfLinearMatrices
      
    EXITS("SolverMappingEMSToSMMap_NumberOfLinearMatricesGet")
    RETURN
999 ERRORSEXITS("SolverMappingEMSToSMMap_NumberOfLinearMatricesGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMappingEMSToSMMap_NumberOfLinearMatricesGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of variables in an equations matrices to solver matrix map.
  SUBROUTINE SolverMappingEMSToSMMap_NumberOfVariablesGet(equationsMatricesToSolverMatrixMap,numberOfVariables,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesToSolverMatrixMapType), POINTER :: equationsMatricesToSolverMatrixMap !<A pointer to the equations matrices to solver matrix map to get the number of variables for
    INTEGER(INTG), INTENT(OUT) :: numberOfVariables !<On exit, the number of variables  mapped to the solver matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingEMSToSMMap_NumberOfVariablesGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(equationsMatricesToSolverMatrixMap)) &
      & CALL FlagError("Equations matrices to solver matrix map is not associated.",err,error,*999)
#endif    
    
    numberOfVariables=equationsMatricesToSolverMatrixMap%numberOfVariables
      
    EXITS("SolverMappingEMSToSMMap_NumberOfVariablesGet")
    RETURN
999 ERRORSEXITS("SolverMappingEMSToSMMap_NumberOfVariablesGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMappingEMSToSMMap_NumberOfVariablesGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the variableIdx'th field variable in an equations matrices to solver matrix map.
  SUBROUTINE SolverMappingEMSToSMMap_VariableGet(equationsMatricesToSolverMatrixMap,variableIdx,fieldVariable,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesToSolverMatrixMapType), POINTER :: equationsMatricesToSolverMatrixMap !<A pointer to the equations matrices to solver matrix map to get the field variable for
    INTEGER(INTG), INTENT(IN) :: variableIdx !<The variable index to get the field variable for
    TYPE(FieldVariableType), POINTER :: fieldVariable !<On exit, a pointer to the specified field variable. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMappingEMSToSMMap_VariableGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsMatricesToSolverMatrixMap)) &
      & CALL FlagError("Equations matrices to solver matrix map is not associated.",err,error,*999)
    IF(variableIdx<1.OR.variableIdx>equationsMatricesToSolverMatrixMap%numberOfVariables) THEN
      localError="The specified variable index of "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
        & " is invalid for the equations matrices to solver matrix map. The variable index should be >= 1 and <= "// &
        & TRIM(NumberToVString(equationsMatricesToSolverMatrixMap%numberOfVariables,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(equationsMatricesToSolverMatrixMap%variables)) &
      & CALL FlagError("The variables array is not allocated for the equations matrices to solver matrix map.",err,error,*999)
#endif    

    fieldVariable=>equationsMatricesToSolverMatrixMap%variables(variableIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(fieldVariable)) THEN
      localError="The variable is not associated for variable index "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
        & " of the equations matrices to solver matrix map."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("SolverMappingEMSToSMMap_VariableGet")
    RETURN
999 NULLIFY(fieldVariable)
998 ERRORS("SolverMappingEMSToSMMap_VariableGet",err,error)
    EXITS("SolverMappingEMSToSMMap_VariableGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingEMSToSMMap_VariableGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the variableIdx'th variable to solver column map in an equations matrices to solver matrix map.
  SUBROUTINE SolverMappingEMSToSMMap_VariableDOFToSolverDOFsMapGet(equationsMatricesToSolverMatrixMap,variableIdx, &
    & varDOFToSolverDOFsMap,err,error,*)

    !Argument variables
    TYPE(EquationsMatricesToSolverMatrixMapType), POINTER :: equationsMatricesToSolverMatrixMap !<A pointer to the equations matrices to solver matrix map to get the variable to solver col map for
    INTEGER(INTG), INTENT(IN) :: variableIdx !<The variable index to get the variable to solver column map for
    TYPE(VariableDOFToSolverDOFsMapType), POINTER :: varDOFToSolverDOFsMap !<On exit, a pointer to the specified variable DOF to solver DOFs map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMappingEMSToSMMap_VariableDOFToSolverDOFsMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(varDOFToSolverDOFsMap)) CALL FlagError("Variable DOF to solver DOFs map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsMatricesToSolverMatrixMap)) &
      & CALL FlagError("Equations matrices to solver matrix map is not associated.",err,error,*999)
    IF(variableIdx<1.OR.variableIdx>equationsMatricesToSolverMatrixMap%numberOfVariables) THEN
      localError="The specified variable index of "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
        & " is invalid for the equations matrices to solver matrix map. The variable index should be >= 1 and <= "// &
        & TRIM(NumberToVString(equationsMatricesToSolverMatrixMap%numberOfVariables,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(equationsMatricesToSolverMatrixMap%varDOFToSolverDOFsMaps)) &
      & CALL FlagError("The variable DOF to solver DOFs map is not allocated for the equations matrices to solver matrix map.", &
      & err,error,*999)
#endif    

    varDOFToSolverDOFsMap=>equationsMatricesToSolverMatrixMap%varDOFToSolverDOFsMaps(variableIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(varDOFToSolverDOFsMap)) THEN
      localError="The variable DOF to solver DOFs map is not associated for variable index "// &
        & TRIM(NumberToVString(variableIdx,"*",err,error))//" of the equations matrices to solver matrix map."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("SolverMappingEMSToSMMap_VariableDOFToSolverDOFsMapGet")
    RETURN
999 NULLIFY(varDOFToSolverDOFsMap)
998 ERRORS("SolverMappingEMSToSMMap_VariableDOFToSolverDOFsMapGet",err,error)
    EXITS("SolverMappingEMSToSMMap_VariableDOFToSolverDOFsMapGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingEMSToSMMap_VariableDOFToSolverDOFsMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the equations col to solver cols map for an equations matrix to solver matrix map.
  SUBROUTINE SolverMappingEMToSMMap_EquationsColToSolverColsMapGet(equationsMatrixToSolverMatrixMap, &
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
 
    ENTERS("SolverMappingEMToSMMap_SolverMatrixGet",err,error,*998)

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
  
  !>Returns a pointer to the equations matrix to solver matrix map for a solver mapping equations matrix to solver matrices map.
  SUBROUTINE SolverMappingEMToSMSMap_EquationsMatrixToSolverMatrixMapGet(equationsMatrixToSolverMatricesMap,solverMatrixIdx, &
    & equationsMatrixToSolverMatrixMap,err,error,*)

    !Argument variables
    TYPE(EquationsMatrixToSolverMatricesMapType), POINTER :: equationsMatrixToSolverMatricesMap !<A pointer to the solver mapping equations matrix to solver matrices map to get the equations matrix to solver matrix map for
    INTEGER(INTG), INTENT(IN) :: solverMatrixIdx !<The solver matrix index to get the equations matrix to solver matrix map for
    TYPE(EquationsMatrixToSolverMatrixMapType), POINTER :: equationsMatrixToSolverMatrixMap !<On exit, a pointer to the specified equations matrix to solver matrix map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMappingEMToSMSMap_EquationsMatrixToSolverMatrixMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(equationsMatrixToSolverMatrixMap)) &
      & CALL FlagError("Equations matrix to solver matrix map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsMatrixToSolverMatricesMap))  &
      & CALL FlagError("Equations matrix to solver matrices map is not associated.",err,error,*999)
    IF(solverMatrixIdx<1.OR.solverMatrixIdx>equationsMatrixToSolverMatricesMap%numberOfSolverMatrices) THEN
      localError="The specified solver matrix index of "//TRIM(NumberToVString(solverMatrixIdx,"*",err,error))// &
        & " is invalid for the equations matrix to solver matrices map. The solver matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(equationsMatrixToSolverMatricesMap%numberOfSolverMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(equationsMatrixToSolverMatricesMap%equationsMatrixToSolverMatrixMaps)) &
      & CALL FlagError("The equations matrix to solver matrix maps is not allocated for the equations matrix to "// &
      & "solver matrices map.",err,error,*999)
#endif    
    
    equationsMatrixToSolverMatrixMap=>equationsMatrixToSolverMatricesMap%equationsMatrixToSolverMatrixMaps(solverMatrixIdx)%ptr
    
#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(equationsMatrixToSolverMatrixMap)) THEN
      localError="The equations matrix to solver matrix map is not associated for the solver matrix index of "// &
        & TRIM(NumberToVString(solverMatrixIdx,"*",err,error))// &
        & " of the equations matrix to solver matrices map of the solver mapping."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
      
    EXITS("SolverMappingEMToSMSMap_EquationsMatrixToSolverMatrixMapGet")
    RETURN
999 NULLIFY(equationsMatrixToSolverMatrixMap)
998 ERRORS("SolverMappingEMToSMSMap_EquationsMatrixToSolverMatrixMapGet",err,error)
    EXITS("SolverMappingEMToSMSMap_EquationsMatrixToSolverMatrixMapGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingEMToSMSMap_EquationsMatrixToSolverMatrixMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of solver matrices for a solver mapping equations matrix to solver matrices map.
  SUBROUTINE SolverMappingEMToSMSMap_NumberOfSolverMatricesGet(equationsMatrixToSolverMatricesMap,numberOfSolverMatrices, &
    & err,error,*)

    !Argument variables
    TYPE(EquationsMatrixToSolverMatricesMapType), POINTER :: equationsMatrixToSolverMatricesMap !<A pointer to the solver mapping equations matrix to solver matrices map to get the number of solver matrices for
    INTEGER(INTG), INTENT(OUT) :: numberOfSolverMatrices !<On exit, the number of solver matrices for the equations matrix to solver matrices map for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingEMToSMSMap_NumberOfSolverMatriceset",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(equationsMatrixToSolverMatricesMap))  &
      & CALL FlagError("Equations matrix to solver matrices map is not associated.",err,error,*999)
#endif    
    
    numberOfSolverMatrices=equationsMatrixToSolverMatricesMap%numberOfSolverMatrices
          
    EXITS("SolverMappingEMToSMSMap_NumberOfSolverMatricesGet")
    RETURN
999 ERRORS("SolverMappingEMToSMSMap_NumberOfSolverMatricesGet",err,error)
    EXITS("SolverMappingEMToSMSMap_NumberOfSolverMatricesGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingEMToSMSMap_NumberOfSolverMatricesGet
  
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
      CALL FlagError(localError,err,error,*999)
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
      CALL FlagError(localError,err,error,*999)
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
  
  !>Returns a pointer to the equations row to solver rows map for a solver mapping equations set to solver matrices map.
  SUBROUTINE SolverMappingESToSMSMap_EquationsRowToSolverRowsMapGet(equationsSetToSolverMatricesMap,equationsRowToSolverRowsMap, &
    & err,error,*)

    !Argument variables
    TYPE(EquationsSetToSolverMatricesMapType), POINTER :: equationsSetToSolverMatricesMap !<A pointer to the solver mapping equations set to solver matrices map to get the equations row tos solver rows map for
    TYPE(MatrixRowColCouplingType), POINTER :: equationsRowToSolverRowsMap(:) !<On exit, a pointer to the equations row to solver rows map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingESToSMSMap_EquationsRowToSolverRowsMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(equationsRowToSolverRowsMap)) &
      & CALL FlagError("Equations row to solver rows map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsSetToSolverMatricesMap))  &
      & CALL FlagError("Equations set to solver matrices map is not associated.",err,error,*999)
#endif    
    
    equationsRowToSolverRowsMap=>equationsSetToSolverMatricesMap%equationsRowToSolverRowsMap
    
#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(equationsRowToSolverRowsMap)) &
      & CALL FlagError("The equations row to solver rows map for the equations set to solver matrices map is not associated.", &
      & err,error,*999)
#endif
      
    EXITS("SolverMappingESToSMSMap_EquationsRowToSolverRowsMapGet")
    RETURN
999 NULLIFY(equationsRowToSolverRowsMap)
998 ERRORS("SolverMappingESToSMSMap_EquationsRowToSolverRowsMapGet",err,error)
    EXITS("SolverMappingESToSMSMap_EquationsRowToSolverRowsMapGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingESToSMSMap_EquationsRowToSolverRowsMapGet
  
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
      CALL FlagError(localError,err,error,*999)
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
  
  !>Returns the number of Jacobian matrices for a solver mapping equations set to solver matrices map.
  SUBROUTINE SolverMappingESToSMSMap_NumberOfJacobianMatricesGet(equationsSetToSolverMatricesMap,numberOfJacobianMatrices, &
    & err,error,*)

    !Argument variables
    TYPE(EquationsSetToSolverMatricesMapType), POINTER :: equationsSetToSolverMatricesMap !<A pointer to the solver mapping equations set to solver matrices map to get the number of Jacobian matrices for
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
  
  !>Returns a pointer to the solver mapping for a solver mapping equations set to solver matrices map.
  SUBROUTINE SolverMappingESToSMSMap_SolverMappingGet(equationsSetToSolverMatricesMap,solverMapping,err,error,*)

    !Argument variables
    TYPE(EquationsSetToSolverMatricesMapType), POINTER :: equationsSetToSolverMatricesMap !<A pointer to the solver mapping equations set to solver matrices map to get the solver mapping for
    TYPE(SolverMappingType), POINTER :: solverMapping !<On exit, a pointer to the solver mapping for the equations set to solver matrices map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingESToSMSMap_SolverMappingGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(solverMapping)) &
      & CALL FlagError("Solver mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsSetToSolverMatricesMap))  &
      & CALL FlagError("Equations set to solver matrices map is not associated.",err,error,*999)
#endif    
    
    solverMapping=>equationsSetToSolverMatricesMap%solverMapping
    
#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(solverMapping)) &
      & CALL FlagError("The solver mapping for the equations set to solver matrices map is not associated.",err,error,*999)
#endif
      
    EXITS("SolverMappingESToSMSMap_SolverMappingGet")
    RETURN
999 NULLIFY(solverMapping)
998 ERRORS("SolverMappingESToSMSMap_SolverMappingGet",err,error)
    EXITS("SolverMappingESToSMSMap_SolverMappingGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingESToSMSMap_SolverMappingGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the interface column to solver rows map for a solver mapping interface condition to solver matrices map.
  SUBROUTINE SolverMappingICToSMSMap_InterfaceColToSolverRowsMapGet(interfaceConditionToSolverMatricesMap, &
    & interfaceColToSolverRowsMap,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionToSolverMatricesMapType), POINTER :: interfaceConditionToSolverMatricesMap !<A pointer to the solver mapping interface condition to solver matrices map to get the interface column to solver rows map for
    TYPE(MatrixRowColCouplingType), POINTER :: interfaceColToSolverRowsMap(:) !<On exit, a pointer to the specified interface column to sovler rows map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMappingICToSMSMap_InterfaceColToSolverRowsMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(interfaceColToSolverRowsMap)) &
      & CALL FlagError("Interface column to solver rows map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceConditionToSolverMatricesMap))  &
      & CALL FlagError("Interface condition to solver matrices map is not associated.",err,error,*999)
#endif    
    
    interfaceColToSolverRowsMap=>interfaceConditionToSolverMatricesMap%interfaceColToSolverRowsMap
    
#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(interfaceColToSolverRowsMap)) THEN
      localError="The interface column to solver rows map is not associated for the interface condition to solver "// &
        & "matrices map of the solver mapping."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
      
    EXITS("SolverMappingICToSMSMap_InterfaceColToSolverRowsMapGet")
    RETURN
999 NULLIFY(interfaceColToSolverRowsMap)
998 ERRORS("SolverMappingICToSMSMap_InterfaceColToSolverRowsMapGet",err,error)
    EXITS("SolverMappingICToSMSMap_InterfaceColToSolverRowsMapGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingICToSMSMap_InterfaceColToSolverRowsMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the interface matrices to solver matrix map for a solver mapping interface condition to solver matrices map.
  SUBROUTINE SolverMappingICToSMSMap_InterfaceMatricesToSolverMatrixMapGet(interfaceConditionToSolverMatricesMap, &
    & solverMatrixIdx,interfaceMatricesToSolverMatrixMap,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionToSolverMatricesMapType), POINTER :: interfaceConditionToSolverMatricesMap !<A pointer to the solver mapping interface condition to solver matrices map to get the interface matrices to solver matrix map for
    INTEGER(INTG), INTENT(IN) :: solverMatrixIdx !<The solver matrix index in to get the interface matrices to solver matrix map for
    TYPE(InterfaceMatricesToSolverMatrixMapType), POINTER :: interfaceMatricesToSolverMatrixMap !<On exit, a pointer to the specified interface matrices to solver matrix map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMappingICToSMSMap_InterfaceMatricesToSolverMatrixMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(interfaceMatricesToSolverMatrixMap)) &
      & CALL FlagError("Interface matrices to solver matrix map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceConditionToSolverMatricesMap))  &
      & CALL FlagError("Interface condition to solver matrices map is not associated.",err,error,*999)
    IF(solverMatrixIdx<1.OR.solverMatrixIdx>interfaceConditionToSolverMatricesMap%numberOfSolverMatrices) THEN
      localError="The specified solver matrix index of "//TRIM(NumberToVString(solverMatrixIdx,"*",err,error))// &
        & " is invalid for the interface condition to solver matrices map. The solver matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(interfaceConditionToSolverMatricesMap%numberOfSolverMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(interfaceConditionToSolverMatricesMap%interfaceMatricesToSolverMatrixMaps)) &
      & CALL FlagError("The interface matrices to solver matrix maps is not allocated for the interface condition to "// &
      & "solver matrices map.",err,error,*999)
#endif    
    
    interfaceMatricesToSolverMatrixMap=>interfaceConditionToSolverMatricesMap% &
      & interfaceMatricesToSolverMatrixMaps(solverMatrixIdx)%ptr
    
#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(interfaceMatricesToSolverMatrixMap)) THEN
      localError="The interface matrices to solver matrix map is not associated for the solver matrix index of "// &
        & TRIM(NumberToVString(solverMatrixIdx,"*",err,error))// &
        & " of the interface condition to solver matrices map of the solver mapping."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
      
    EXITS("SolverMappingICToSMSMap_InterfaceMatricesToSolverMatrixMapGet")
    RETURN
999 NULLIFY(interfaceMatricesToSolverMatrixMap)
998 ERRORS("SolverMappingICToSMSMap_InterfaceMatricesToSolverMatrixMapGet",err,error)
    EXITS("SolverMappingICToSMSMap_InterfaceMatricesToSolverMatrixMapGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingICToSMSMap_InterfaceMatricesToSolverMatrixMapGet
  
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
    IF(interfaceMatrixIdx<1.OR.interfaceMatrixIdx>interfaceConditionToSolverMatricesMap%numberOfInterfaceMatrices) THEN
      localError="The specified interface matrix index of "//TRIM(NumberToVString(interfaceMatrixIdx,"*",err,error))// &
        & " is invalid for the interface condition to solver matrices map. The interface matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(interfaceConditionToSolverMatricesMap%numberOfInterfaceMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(interfaceConditionToSolverMatricesMap%interfaceMatrixToSolverMatricesMaps)) &
      & CALL FlagError("The interface matrix to solver matrices maps is not allocated for the interface condition to "// &
      & "solver matrices map.",err,error,*999)
#endif    
    
    interfaceMatrixToSolverMatricesMap=>interfaceConditionToSolverMatricesMap% &
      & interfaceMatrixToSolverMatricesMaps(interfaceMatrixIdx)%ptr
    
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
  SUBROUTINE SolverMappingIMToSMMap_InterfaceRowToSolverColsMapGet(interfaceMatrixToSolverMatrixMap, &
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
  SUBROUTINE SolverMappingIMToSMMap_SolverMatrixGet(interfaceMatrixToSolverMatrixMap,solverMatrix,err,error,*)

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
  
  !>Returns a pointer to the variableIdx'th dependent variable for a solver mapping interface matrices to solver matrix map.
  SUBROUTINE SolverMappingIMSToSMMap_DependentVariableGet(interfaceMatricesToSolverMatrixMap,variableIdx,dependentVariable, &
    & err,error,*)

    !Argument variables
    TYPE(InterfaceMatricesToSolverMatrixMapType), POINTER :: interfaceMatricesToSolverMatrixMap !<A pointer to the solver mapping interface matrices to solver matrix map to get the dependent variable for
    INTEGER(INTG), INTENT(IN) :: variableIdx !<The dependent variable index in to get
    TYPE(FieldVariableType), POINTER :: dependentVariable !<On exit, a pointer to the specified dependent variable. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMappingIMSToSMMap_DependentVariableGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dependentVariable)) &
      & CALL FlagError("Dependent variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMatricesToSolverMatrixMap))  &
      & CALL FlagError("Interface matrices to solver matrix map is not associated.",err,error,*999)
    IF(variableIdx<1.OR.variableIdx>interfaceMatricesToSolverMatrixMap%numberOfDependentVariables) THEN
      localError="The specified dependent variable index of "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
        & " is invalid for the interface matrices to solver matrix map. The dependent variable index should be >= 1 and <= "// &
        & TRIM(NumberToVString(interfaceMatricesToSolverMatrixMap%numberOfDependentVariables,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(interfaceMatricesToSolverMatrixMap%dependentVariables)) &
      & CALL FlagError("The dependent variables array is not allocated for the interface matrices to solver matrix map.", &
      & err,error,*999)
#endif    
    
    dependentVariable=>interfaceMatricesToSolverMatrixMap%dependentVariables(variableIdx)%ptr
    
#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(dependentVariable)) THEN
      localError="The dependent variable is not associated for the dependent variable index of "// &
        & TRIM(NumberToVString(variableIdx,"*",err,error))// &
        & " of the interface matrices to solver matrix map of the solver mapping."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
      
    EXITS("SolverMappingIMSToSMMap_DependentVariableGet")
    RETURN
999 NULLIFY(dependentVariable)
998 ERRORS("SolverMappingIMSToSMMap_DependentVariableGet",err,error)
    EXITS("SolverMappingIMSToSMMap_DependentVariableGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingIMSToSMMap_DependentVariableGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the dependent variable DOF to solver DOFs map for a solver mapping interface matrices to solver matrix map.
  SUBROUTINE SolverMappingIMSToSMMap_DependentVarDOFToSolverDOFsMapGet(interfaceMatricesToSolverMatrixMap,interfaceMatrixIdx, &
    & dependentVarDOFToSolverDOFsMap,err,error,*)

    !Argument variables
    TYPE(InterfaceMatricesToSolverMatrixMapType), POINTER :: interfaceMatricesToSolverMatrixMap !<A pointer to the solver mapping interface matrices to solver matrix map to get the dependent variable DOF to solver DOFs map for
    INTEGER(INTG), INTENT(IN) :: interfaceMatrixIdx !<The interface matrix index to get the dependent variable DOF to solver DOFs map
    TYPE(VariableDOFToSolverDOFsMapType), POINTER :: dependentVarDOFToSolverDOFsMap !<On exit, a pointer to the dependent variable DOF to solver DOFs map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif
 
    ENTERS("SolverMappingIMSToSMMap_DependentVarDOFToSolverDOFsMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dependentVarDOFToSolverDOFsMap)) &
      & CALL FlagError("Dependent variable DOF to solver DOFs map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMatricesToSolverMatrixMap))  &
      & CALL FlagError("Interface matrices to solver matrix map is not associated.",err,error,*999)
    IF(interfaceMatrixIdx<1.OR.interfaceMatrixIdx>interfaceMatricesToSolverMatrixMap%numberOfInterfaceMatrices) THEN
      localError="The specified interface matrix index of "//TRIM(NumberToVString(interfaceMatrixIdx,"*",err,error))// &
        & " is invalid. The interface matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(interfaceMatricesToSolverMatrixMap%numberOfInterfaceMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)      
    ENDIF
    IF(.NOT.ALLOCATED(interfaceMatricesToSolverMatrixMap%dependentVarDOFToSolverDOFsMaps)) &
      & CALL FlagError("The dependent variable DOF to solver DOFs map is not allocated for the interface matrices to solver"// &
      & " matrix map.",err,error,*999)
#endif    
    
    dependentVarDOFToSolverDOFsMap=>interfaceMatricesToSolverMatrixMap%dependentVarDOFToSolverDOFsMaps(interfaceMatrixIdx)%ptr
    
#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(dependentVarDOFToSolverDOFsMap)) THEN
      localError="The dependent variable DOF to solver DOFs map is not associated for "// &
        & " interface matrix index "//TRIM(NumberToVString(interfaceMatrixIdx,"*",err,error))// &
        & " of the interface matrices to solver matrix map."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
      
    EXITS("SolverMappingIMSToSMMap_DependentVarDOFToSolverDOFsMapGet")
    RETURN
999 NULLIFY(dependentVarDOFToSolverDOFsMap)
998 ERRORS("SolverMappingIMSToSMMap_DependentVarDOFToSolverDOFsMapGet",err,error)
    EXITS("SolverMappingIMSToSMMap_DependentVarDOFToSolverDOFsMapGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingIMSToSMMap_DependentVarDOFToSolverDOFsMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the interface column to solver columns map for a solver mapping interface matrices to solver matrix map.
  SUBROUTINE SolverMappingIMSToSMMap_InterfaceColToSolverColsMapGet(interfaceMatricesToSolverMatrixMap, &
    & interfaceColToSolverColsMap,err,error,*)

    !Argument variables
    TYPE(InterfaceMatricesToSolverMatrixMapType), POINTER :: interfaceMatricesToSolverMatrixMap !<A pointer to the solver mapping interface matrices to solver matrix map to get the interface column to solver columns map for
    TYPE(MatrixRowColCouplingType), POINTER :: interfaceColToSolverColsMap(:) !<On exit, a pointer to the interface column to solver columns map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMappingIMSToSMMap_InterfaceColToSolverColsMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(interfaceColToSolverColsMap)) &
      & CALL FlagError("Interface column to solver columns map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMatricesToSolverMatrixMap))  &
      & CALL FlagError("Interface matrices to solver matrix map is not associated.",err,error,*999)
#endif    
    
    interfaceColToSolverColsMap=>interfaceMatricesToSolverMatrixMap%interfaceColToSolverColsMap
    
#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(interfaceColToSolverColsMap)) THEN
      localError="The interface column to solver columns map is not associated for the interface matrices to solver "// &
        & "matrix map of the solver mapping."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
      
    EXITS("SolverMappingIMSToSMMap_InterfaceColToSolverColsMapGet")
    RETURN
999 NULLIFY(interfaceColToSolverColsMap)
998 ERRORS("SolverMappingIMSToSMMap_InterfaceColToSolverColsMapGet",err,error)
    EXITS("SolverMappingIMSToSMMap_InterfaceColToSolverColsMapGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingIMSToSMMap_InterfaceColToSolverColsMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the interface matrix to solver matrix map for a solver mapping interface matrices to solver matrix map.
  SUBROUTINE SolverMappingIMSToSMMap_InterfaceMatrixToSolverMatrixMapGet(interfaceMatricesToSolverMatrixMap, &
    & interfaceMatrixIdx,interfaceMatrixToSolverMatrixMap,err,error,*)

    !Argument variables
    TYPE(InterfaceMatricesToSolverMatrixMapType), POINTER :: interfaceMatricesToSolverMatrixMap !<A pointer to the solver mapping interface matrices to solver matrix map to get the interface matrix to solver matrix map for
    INTEGER(INTG), INTENT(IN) :: interfaceMatrixIdx !<The interface matrix index in to get the interface matrix to solver matrix map for
    TYPE(InterfaceMatrixToSolverMatrixMapType), POINTER :: interfaceMatrixToSolverMatrixMap !<On exit, a pointer to the specified interface matrix to solver matrix map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMappingIMSToSMMap_InterfaceMatrixToSolverMatrixMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(interfaceMatrixToSolverMatrixMap)) &
      & CALL FlagError("Interface matrix to solver matrix map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMatricesToSolverMatrixMap))  &
      & CALL FlagError("Interface matrices to solver matrix map is not associated.",err,error,*999)
    IF(interfaceMatrixIdx<1.OR.interfaceMatrixIdx>interfaceMatricesToSolverMatrixMap%numberOfInterfaceMatrices) THEN
      localError="The specified interface matrix index of "//TRIM(NumberToVString(interfaceMatrixIdx,"*",err,error))// &
        & " is invalid for the interface matrices to solver matrix map. The interface matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(interfaceMatricesToSolverMatrixMap%numberOfInterfaceMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(interfaceMatricesToSolverMatrixMap%interfaceMatrixToSolverMatrixMaps)) &
      & CALL FlagError("The interface matrix to solver matrix maps is not allocated for the interface matrices to "// &
      & "solver matrix map.",err,error,*999)
#endif    
    
    interfaceMatrixToSolverMatrixMap=>interfaceMatricesToSolverMatrixMap% &
      & interfaceMatrixToSolverMatrixMaps(interfaceMatrixIdx)%ptr
    
#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(interfaceMatrixToSolverMatrixMap)) THEN
      localError="The interface matrix to solver matrix map is not associated for the interface matrix index of "// &
        & TRIM(NumberToVString(interfaceMatrixIdx,"*",err,error))// &
        & " of the interface matrices to solver matrix map of the solver mapping."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
      
    EXITS("SolverMappingIMSToSMMap_InterfaceMatrixToSolverMatrixMapGet")
    RETURN
999 NULLIFY(interfaceMatrixToSolverMatrixMap)
998 ERRORS("SolverMappingIMSToSMMap_InterfaceMatrixToSolverMatrixMapGet",err,error)
    EXITS("SolverMappingIMSToSMMap_InterfaceMatrixToSolverMatrixMapGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingIMSToSMMap_InterfaceMatrixToSolverMatrixMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the Lagrange variable for a solver mapping interface matrices to solver matrix map.
  SUBROUTINE SolverMappingIMSToSMMap_LagrangeVariableGet(interfaceMatricesToSolverMatrixMap,lagrangeVariable,err,error,*)

    !Argument variables
    TYPE(InterfaceMatricesToSolverMatrixMapType), POINTER :: interfaceMatricesToSolverMatrixMap !<A pointer to the solver mapping interface matrices to solver matrix map to get the Lagrange variable for
    TYPE(FieldVariableType), POINTER :: lagrangeVariable !<On exit, a pointer to the Lagrange variable. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingIMSToSMMap_LagrangeVariableGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(lagrangeVariable)) CALL FlagError("Lagrange variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMatricesToSolverMatrixMap))  &
      & CALL FlagError("Interface matrices to solver matrix map is not associated.",err,error,*999)
#endif    
    
    lagrangeVariable=>interfaceMatricesToSolverMatrixMap%lagrangeVariable
    
#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(lagrangeVariable)) &
      & CALL FlagError("The Lagrange variable is not associated for the interface matrices to solver matrix map.",err,error,*999)
#endif
      
    EXITS("SolverMappingIMSToSMMap_LagrangeVariableGet")
    RETURN
999 NULLIFY(lagrangeVariable)
998 ERRORS("SolverMappingIMSToSMMap_LagrangeVariableGet",err,error)
    EXITS("SolverMappingIMSToSMMap_LagrangeVariableGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingIMSToSMMap_LagrangeVariableGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the Lagrange variable DOF to solver DOFs map for a solver mapping interface matrices to solver matrix map.
  SUBROUTINE SolverMappingIMSToSMMap_LagrangeVarDOFToSolverDOFsMapGet(interfaceMatricesToSolverMatrixMap, &
    & lagrangeVarDOFToSolverDOFsMap,err,error,*)

    !Argument variables
    TYPE(InterfaceMatricesToSolverMatrixMapType), POINTER :: interfaceMatricesToSolverMatrixMap !<A pointer to the solver mapping interface matrices to solver matrix map to get the Lagrange variable to solver column map for
    TYPE(VariableDOFToSolverDOFsMapType), POINTER :: lagrangeVarDOFToSolverDOFsMap !<On exit, a pointer to the Lagrange variable DOF to solver DOFs map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingIMSToSMMap_LagrangeVarDOFToSolverDOFsMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(lagrangeVarDOFToSolverDOFsMap)) &
      & CALL FlagError("Lagrange variable DOF to solver DOFs map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMatricesToSolverMatrixMap))  &
      & CALL FlagError("Interface matrices to solver matrix map is not associated.",err,error,*999)
#endif    
    
    lagrangeVarDOFToSolverDOFsMap=>interfaceMatricesToSolverMatrixMap%lagrangeVarDOFToSolverDOFsMap
    
#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(lagrangeVarDOFToSolverDOFsMap)) &
      & CALL FlagError("The Lagrange variable DOF to solver DOFs map is not associated for the interface matrices to solver "// &
      & "matrix map.",err,error,*999)
#endif
      
    EXITS("SolverMappingIMSToSMMap_LagrangeVarDOFToSolverDOFsMapGet")
    RETURN
999 NULLIFY(lagrangeVarDOFToSolverDOFsMap)
998 ERRORS("SolverMappingIMSToSMMap_LagrangeVarDOFToSolverDOFsMapGet",err,error)
    EXITS("SolverMappingIMSToSMMap_LagrangeVarDOFToSolverDOFsMapGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingIMSToSMMap_LagrangeVarDOFToSolverDOFsMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of interface matrices for a solver mapping interface matrices to solver matrix map.
  SUBROUTINE SolverMappingIMSToSMMap_NumberOfInterfaceMatricesGet(interfaceMatricesToSolverMatrixMap,numberOfInterfaceMatrices, &
    & err,error,*)

    !Argument variables
    TYPE(InterfaceMatricesToSolverMatrixMapType), POINTER :: interfaceMatricesToSolverMatrixMap !<A pointer to the solver mapping interface matrices to solver matrix map to get the number of interface matrices for
    INTEGER(INTG), INTENT(OUT) :: numberOfInterfaceMatrices !<On exit, the number of interface matrices for the interface matrices to solver matrix map
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingIMSToSMMap_NumberOfInterfaceMatricesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(interfaceMatricesToSolverMatrixMap))  &
      & CALL FlagError("Interface matrices to solver matrix map is not associated.",err,error,*999)
#endif    
    
    numberOfInterfaceMatrices=interfaceMatricesToSolverMatrixMap%numberOfInterfaceMatrices
          
    EXITS("SolverMappingIMSToSMMap_NumberOfInterfaceMatricesGet")
    RETURN
999 ERRORS("SolverMappingIMSToSMMap_NumberOfInterfaceMatricesGet",err,error)
    EXITS("SolverMappingIMSToSMMap_NumberOfInterfaceMatricesGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingIMSToSMMap_NumberOfInterfaceMatricesGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of dependent variables for a solver mapping interface matrices to solver matrix map.
  SUBROUTINE SolverMappingIMSToSMMap_NumberOfDependentVariablesGet(interfaceMatricesToSolverMatrixMap,numberOfDependentVariables, &
    & err,error,*)

    !Argument variables
    TYPE(InterfaceMatricesToSolverMatrixMapType), POINTER :: interfaceMatricesToSolverMatrixMap !<A pointer to the solver mapping interface matrices to solver matrix map to get the number of dependent variables for
    INTEGER(INTG), INTENT(OUT) :: numberOfDependentVariables !<On exit, the number of dependent variables for the interface matrices to solver matrix map
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingIMSToSMMap_NumberOfDependentVariablesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(interfaceMatricesToSolverMatrixMap))  &
      & CALL FlagError("Interface matrices to solver matrix map is not associated.",err,error,*999)
#endif    
    
    numberOfDependentVariables=interfaceMatricesToSolverMatrixMap%numberOfInterfaceMatrices
          
    EXITS("SolverMappingIMSToSMMap_NumberOfDependentVariablesGet")
    RETURN
999 ERRORS("SolverMappingIMSToSMMap_NumberOfDependentVariablesGet",err,error)
    EXITS("SolverMappingIMSToSMMap_NumberOfDependentVariablesGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingIMSToSMMap_NumberOfDependentVariablesGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the interface matrix to solver matrix map for a solver matrix index in an interface matrix to a solver matrices map.
  SUBROUTINE SolverMappingIMToSMSMap_InterfaceMatrixToSolverMatrixMapGet(interfaceMatrixToSolverMatricesMap, &
    & solverMatrixIdx,interfaceMatrixToSolverMatrixMap,err,error,*)

    !Argument variables
    TYPE(InterfaceMatrixToSolverMatricesMapType), POINTER :: interfaceMatrixToSolverMatricesMap !<A pointer to the interface matrix to solver matrices map to get the interface matrix to solver matrix map for
    INTEGER(INTG), INTENT(IN) :: solverMatrixIdx !<The solver matrix index to get the interface matrix to solver matrix map for.
    TYPE(InterfaceMatrixToSolverMatrixMapType), POINTER :: interfaceMatrixToSolverMatrixMap !<On exit, the interface matrix to solver matrix map for the specified solver matrix index for the interface matrix to solver matrices map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMappingIMToSMSMap_InterfaceMatrixToSolverMatrixMapGet",err,error,*998)

#ifdef WITH_PRECHECKS
    IF(ASSOCIATED(interfaceMatrixToSolverMatrixMap)) &
      & CALL FlagError("Interface matrix to solver matrix map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMatrixToSolverMatricesMap)) &
      & CALL FlagError("Interface matrix to solver matrices map is not associated.",err,error,*999)
    IF(solverMatrixIdx<1.OR.solverMatrixIdx>interfaceMatrixToSolverMatricesMap%numberOfSolverMatrices) THEN
      localError="The specified solver matrix index of "//TRIM(NumberToVString(solverMatrixIdx,"*",err,error))// &
        & " is invalid for the interface matrices to solver matrix map. The solver matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(interfaceMatrixToSolverMatricesMap%numberOfSolverMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(interfaceMatrixToSolverMatricesMap%interfaceMatrixToSolverMatrixMaps)) &
      & CALL FlagError("The interface matrix to solver matrix maps array is not allocated for the interface matrix "// &
      & "to solver matrices map.",err,error,*999)
#endif    
    
    interfaceMatrixToSolverMatrixMap=>interfaceMatrixToSolverMatricesMap%interfaceMatrixToSolverMatrixMaps(solverMatrixIdx)%ptr

#ifdef WITH_POSTCHECKS
    IF(.NOT.ASSOCIATED(interfaceMatrixToSolverMatrixMap)) THEN
      localError="The interface matrix to solver matrix map is not associated for solver matrix index "// &
        & TRIM(NumberToVString(solverMatrixIdx,"*",err,error))//" of the interface matrix to solver matrices map."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("SolverMappingIMToSMSMap_InterfaceMatrixToSolverMatrixMapGet")
    RETURN
999 NULLIFY(interfaceMatrixToSolverMatrixMap)
998 ERRORSEXITS("SolverMappingIMToSMSMap_InterfaceMatrixToSolverMatrixMapGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMappingIMToSMSMap_InterfaceMatrixToSolverMatrixMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the interface row to solver rows map in an interface matrix to a solver matrices map.
  SUBROUTINE SolverMappingIMToSMSMap_InterfaceRowToSolverRowsMapGet(interfaceMatrixToSolverMatricesMap, &
    & interfaceRowToSolverRowsMap,err,error,*)

    !Argument variables
    TYPE(InterfaceMatrixToSolverMatricesMapType), POINTER :: interfaceMatrixToSolverMatricesMap !<A pointer to the interface matrix to solver matrices map to get the interface row to solver rows map for
    TYPE(MatrixRowColCouplingType), POINTER :: interfaceRowToSolverRowsMap(:) !<On exit, the interface row to solver rows map for the interface matrix to solver matrices map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingIMToSMSMap_InterfaceRowToSolverRowsMapGet",err,error,*998)

#ifdef WITH_PRECHECKS
    IF(ASSOCIATED(interfaceRowToSolverRowsMap)) &
      & CALL FlagError("Interface row to solver rows map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceMatrixToSolverMatricesMap)) &
      & CALL FlagError("Interface matrix to solver matrices map is not associated.",err,error,*999)
#endif    
    
    interfaceRowToSolverRowsMap=>interfaceMatrixToSolverMatricesMap%interfaceRowToSolverRowsMap

#ifdef WITH_POSTCHECKS
    IF(.NOT.ASSOCIATED(interfaceRowToSolverRowsMap)) &
      CALL FlagError("The interface row to solver rows map is not associated for the interface matrices to solver matrix map.", &
      & err,error,*999)
#endif    
      
    EXITS("SolverMappingIMToSMSMap_InterfaceRowToSolverRowsMapGet")
    RETURN
999 NULLIFY(interfaceRowToSolverRowsMap)
998 ERRORSEXITS("SolverMappingIMToSMSMap_InterfaceRowToSolverRowsMapGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMappingIMToSMSMap_InterfaceRowToSolverRowsMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of solver matrices in an interface matrix to a solver matrices map.
  SUBROUTINE SolverMappingIMToSMSMap_NumberOfSolverMatricesGet(interfaceMatrixToSolverMatricesMap,numberOfSolverMatrices, &
    & err,error,*)

    !Argument variables
    TYPE(InterfaceMatrixToSolverMatricesMapType), POINTER :: interfaceMatrixToSolverMatricesMap !<A pointer to the interface matrix to solver matrices map to get the number of solver matrices for
    INTEGER(INTG), INTENT(OUT) :: numberOfSolverMatrices !<On exit, the number of solver matrices in the interface matrix to solver matrices mapping.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingIMToSMSMap_NumberOfSolverMatricesGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(interfaceMatrixToSolverMatricesMap)) &
      & CALL FlagError("Interface matrix to solver matrices map is not associated.",err,error,*999)
#endif    
    
    numberOfSolverMatrices=interfaceMatrixToSolverMatricesMap%numberOfSolverMatrices
    
    EXITS("SolverMappingIMToSMSMap_NumberOfSolverMatricesGet")
    RETURN
999 ERRORSEXITS("SolverMappingIMToSMSMap_NumberOfSolverMatricesGet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverMappingIMToSMSMap_NumberOfSolverMatricesGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the Jacobian col to solver cols map for an Jacobian matrix to solver matrix map.
  SUBROUTINE SolverMappingJMToSMMap_JacobianColToSolverColsMapGet(jacobianMatrixToSolverMatrixMap, &
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
  
  !>Returns the dynamic matrix coupling information for a dynamic matrix in solver column to dynamic equations map.
  SUBROUTINE SolverMappingSColToDEQSMap_DynamicCouplingInfoGet(solverColToDynamicEquationsMap,dynamicMatrixIdx, &
    & dynamicMatrixNumber,columnNumber,couplingCoefficient,err,error,*)

    !Argument variables
    TYPE(SolverColToDynamicEquationsMapType), POINTER :: solverColToDynamicEquationsMap !<A pointer to the solver column to dynamic equations map to get the dynamic coupling information for
    INTEGER(INTG), INTENT(IN) :: dynamicMatrixIdx !<The dynamic matrix index to get the coupling information for
    INTEGER(INTG), INTENT(OUT) :: dynamicMatrixNumber !<On exit, the number of the dynamic matrix in the solver column to dynamic equations map.
    INTEGER(INTG), INTENT(OUT) :: columnNumber !<On exit, the column number in the dynamic matrix the solver column is coupled to
    REAL(DP), INTENT(OUT) :: couplingCoefficient !<On exit, the coupling coefficient applied to the solver column to dynamic matrix column
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif
    
    ENTERS("SolverMappingSColToDEQSMap_DynamicCouplingInfoGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(solverColToDynamicEquationsMap)) &
      & CALL FlagError("The solver column to dynamic equations map is not associated.",err,error,*999)
    IF(dynamicMatrixIdx<1.OR.dynamicMatrixIdx>solverColToDynamicEquationsMap%numberOfDynamicMatrices) THEN
      localError="The specified dynamic matrix index of "//TRIM(NumberToVString(dynamicMatrixIdx,"*",err,error))// &
        & " is invalid. The dynamic matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(solverColToDynamicEquationsMap%numberOfDynamicMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(solverColToDynamicEquationsMap%equationsMatrixNumbers)) &
      & CALL FlagError("The equations matrix numbers array is not allocated for the solver column to dynamic matrices map.", &
      & err,error,*999)
    IF(.NOT.ALLOCATED(solverColToDynamicEquationsMap%equationsColumnNumbers)) &
      & CALL FlagError("The equations column numbers array is not allocated for the solver column to dynamic matrices map.", &
      & err,error,*999)
    IF(.NOT.ALLOCATED(solverColToDynamicEquationsMap%couplingCoefficients)) &
      & CALL FlagError("The coupling coefficients array is not allocated for the solver column to dynamic matrices map.", &
      & err,error,*999)
#endif    

    dynamicMatrixNumber=solverColToDynamicEquationsMap%equationsMatrixNumbers(dynamicMatrixIdx)
    columnNumber=solverColToDynamicEquationsMap%equationsColumnNumbers(dynamicMatrixIdx)
    couplingCoefficient=solverColToDynamicEquationsMap%couplingCoefficients(dynamicMatrixIdx)
      
    EXITS("SolverMappingSColToDEQSMap_DynamicCouplingInfoGet")
    RETURN
999 dynamicMatrixNumber=0
    columnNumber=0
    couplingCoefficient=0.0_DP
    ERRORS("SolverMappingSColToDEQSMap_DynamicCouplingInfoGet",err,error)
    EXITS("SolverMappingSColToDEQSMap_DynamicCouplingInfoGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingSColToDEQSMap_DynamicCouplingInfoGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of dynamic matrices for a solver column to dynamic equations map.
  SUBROUTINE SolverMappingSColToDEQSMap_NumberOfDynamicMatricesGet(solverColToDynamicEquationsMap,numberOfDynamicMatrices, &
    & err,error,*)

    !Argument variables
    TYPE(SolverColToDynamicEquationsMapType), POINTER :: solverColToDynamicEquationsMap !<A pointer to the solver column to dynamic equations map to get the number of dynamic matrices for
    INTEGER(INTG), INTENT(OUT) :: numberOfDynamicMatrices !<On exit, the number of dynamic equations matrices in the solver column to dynamic equations map.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingSColToDEQSMap_NumberOfDynamicMatricesGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(solverColToDynamicEquationsMap)) &
      & CALL FlagError("The solver column to dynamic equations map is not associated.",err,error,*999)
#endif    

    numberOfDynamicMatrices=solverColToDynamicEquationsMap%numberOfDynamicMatrices
      
    EXITS("SolverMappingSColToDEQSMap_NumberOfDynamicMatricesGet")
    RETURN
999 numberOfDynamicMatrices=0
    ERRORS("SolverMappingSColToDEQSMap_NumberOfDynamicMatricesGet",err,error)
    EXITS("SolverMappingSColToDEQSMap_NumberOfDynamicMatricesGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingSColToDEQSMap_NumberOfDynamicMatricesGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the interface matrix coupling information for a interface matrix in solver column to interface equations map.
  SUBROUTINE SolverMappingSColToIEQSMap_InterfaceCouplingInfoGet(solverColToInterfaceEquationsMap,interfaceMatrixIdx, &
    & interfaceMatrixNumber,columnNumber,couplingCoefficient,err,error,*)

    !Argument variables
    TYPE(SolverColToInterfaceEquationsMapType), POINTER :: solverColToInterfaceEquationsMap !<A pointer to the solver column to interface equations map to get the interface coupling information for
    INTEGER(INTG), INTENT(IN) :: interfaceMatrixIdx !<The interface matrix index to get the coupling information for
    INTEGER(INTG), INTENT(OUT) :: interfaceMatrixNumber !<On exit, the number of the interface matrix in the solver column to interface equations map.
    INTEGER(INTG), INTENT(OUT) :: columnNumber !<On exit, the column number in the interface matrix the solver column is coupled to
    REAL(DP), INTENT(OUT) :: couplingCoefficient !<On exit, the coupling coefficient applied to the solver column to interface matrix column
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif
    
    ENTERS("SolverMappingSColToIEQSMap_InterfaceCouplingInfoGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(solverColToInterfaceEquationsMap)) &
      & CALL FlagError("The solver column to interface equations map is not associated.",err,error,*999)
    IF(interfaceMatrixIdx<1.OR.interfaceMatrixIdx>solverColToInterfaceEquationsMap%numberOfInterfaceMatrices) THEN
      localError="The specified interface matrix index of "//TRIM(NumberToVString(interfaceMatrixIdx,"*",err,error))// &
        & " is invalid. The interface matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(solverColToInterfaceEquationsMap%numberOfInterfaceMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(solverColToInterfaceEquationsMap%interfaceMatrixNumbers)) &
      & CALL FlagError("The interface matrix numbers array is not allocated for the solver column to interface matrices map.", &
      & err,error,*999)
    IF(.NOT.ALLOCATED(solverColToInterfaceEquationsMap%interfaceColumnNumbers)) &
      & CALL FlagError("The interface column numbers array is not allocated for the solver column to interface matrices map.", &
      & err,error,*999)
    IF(.NOT.ALLOCATED(solverColToInterfaceEquationsMap%couplingCoefficients)) &
      & CALL FlagError("The coupling coefficients array is not allocated for the solver column to interface matrices map.", &
      & err,error,*999)
#endif    

    interfaceMatrixNumber=solverColToInterfaceEquationsMap%interfaceMatrixNumbers(interfaceMatrixIdx)
    columnNumber=solverColToInterfaceEquationsMap%interfaceColumnNumbers(interfaceMatrixIdx)
    couplingCoefficient=solverColToInterfaceEquationsMap%couplingCoefficients(interfaceMatrixIdx)
      
    EXITS("SolverMappingSColToIEQSMap_InterfaceCouplingInfoGet")
    RETURN
999 interfaceMatrixNumber=0
    columnNumber=0
    couplingCoefficient=0.0_DP
    ERRORS("SolverMappingSColToIEQSMap_InterfaceCouplingInfoGet",err,error)
    EXITS("SolverMappingSColToIEQSMap_InterfaceCouplingInfoGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingSColToIEQSMap_InterfaceCouplingInfoGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of interface matrices for a solver column to interface equations map.
  SUBROUTINE SolverMappingSColToIEQSMap_NumberOfInterfaceMatricesGet(solverColToInterfaceEquationsMap,numberOfInterfaceMatrices, &
    & err,error,*)

    !Argument variables
    TYPE(SolverColToInterfaceEquationsMapType), POINTER :: solverColToInterfaceEquationsMap !<A pointer to the solver column to interface equations map to get the number of interface matrices for
    INTEGER(INTG), INTENT(OUT) :: numberOfInterfaceMatrices !<On exit, the number of interface matrices in the solver column to interface equations map.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingSColToIEQSMap_NumberOfInterfaceMatricesGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(solverColToInterfaceEquationsMap)) &
      & CALL FlagError("The solver column to interface equations map is not associated.",err,error,*999)
#endif    

    numberOfInterfaceMatrices=solverColToInterfaceEquationsMap%numberOfInterfaceMatrices
      
    EXITS("SolverMappingSColToIEQSMap_NumberOfInterfaceMatricesGet")
    RETURN
999 numberOfInterfaceMatrices=0
    ERRORS("SolverMappingSColToIEQSMap_NumberOfInterfaceMatricesGet",err,error)
    EXITS("SolverMappingSColToIEQSMap_NumberOfInterfaceMatricesGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingSColToIEQSMap_NumberOfInterfaceMatricesGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the linear matrix coupling information for a linear matrix in solver column to linear equations map.
  SUBROUTINE SolverMappingSColToLEQSMap_LinearCouplingInfoGet(solverColToLinearEquationsMap,linearMatrixIdx, &
    & linearMatrixNumber,columnNumber,couplingCoefficient,err,error,*)

    !Argument variables
    TYPE(SolverColToLinearEquationsMapType), POINTER :: solverColToLinearEquationsMap !<A pointer to the solver column to linear equations map to get the linear coupling information for
    INTEGER(INTG), INTENT(IN) :: linearMatrixIdx !<The linear matrix index to get the coupling information for
    INTEGER(INTG), INTENT(OUT) :: linearMatrixNumber !<On exit, the number of the linear matrix in the solver column to linear equations map.
    INTEGER(INTG), INTENT(OUT) :: columnNumber !<On exit, the column number in the linear matrix the solver column is coupled to
    REAL(DP), INTENT(OUT) :: couplingCoefficient !<On exit, the coupling coefficient applied to the solver column to linear matrix column
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif
    
    ENTERS("SolverMappingSColToLEQSMap_LinearCouplingInfoGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(solverColToLinearEquationsMap)) &
      & CALL FlagError("The solver column to linear equations map is not associated.",err,error,*999)
    IF(linearMatrixIdx<1.OR.linearMatrixIdx>solverColToLinearEquationsMap%numberOfLinearMatrices) THEN
      localError="The specified linear matrix index of "//TRIM(NumberToVString(linearMatrixIdx,"*",err,error))// &
        & " is invalid. The linear matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(solverColToLinearEquationsMap%numberOfLinearMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(solverColToLinearEquationsMap%equationsMatrixNumbers)) &
      & CALL FlagError("The equations matrix numbers array is not allocated for the solver column to linear matrices map.", &
      & err,error,*999)
    IF(.NOT.ALLOCATED(solverColToLinearEquationsMap%equationsColumnNumbers)) &
      & CALL FlagError("The equations column numbers array is not allocated for the solver column to linear matrices map.", &
      & err,error,*999)
    IF(.NOT.ALLOCATED(solverColToLinearEquationsMap%couplingCoefficients)) &
      & CALL FlagError("The coupling coefficients array is not allocated for the solver column to linear matrices map.", &
      & err,error,*999)
#endif    

    linearMatrixNumber=solverColToLinearEquationsMap%equationsMatrixNumbers(linearMatrixIdx)
    columnNumber=solverColToLinearEquationsMap%equationsColumnNumbers(linearMatrixIdx)
    couplingCoefficient=solverColToLinearEquationsMap%couplingCoefficients(linearMatrixIdx)
      
    EXITS("SolverMappingSColToLEQSMap_LinearCouplingInfoGet")
    RETURN
999 linearMatrixNumber=0
    columnNumber=0
    couplingCoefficient=0.0_DP
    ERRORS("SolverMappingSColToLEQSMap_LinearCouplingInfoGet",err,error)
    EXITS("SolverMappingSColToLEQSMap_LinearCouplingInfoGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingSColToLEQSMap_LinearCouplingInfoGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of linear matrices for a solver column to linear equations map.
  SUBROUTINE SolverMappingSColToLEQSMap_NumberOfLinearMatricesGet(solverColToLinearEquationsMap,numberOfLinearMatrices, &
    & err,error,*)

    !Argument variables
    TYPE(SolverColToLinearEquationsMapType), POINTER :: solverColToLinearEquationsMap !<A pointer to the solver column to linear equations map to get the number of linear matrices for
    INTEGER(INTG), INTENT(OUT) :: numberOfLinearMatrices !<On exit, the number of linear matrices in the solver column to linear equations map.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingSColToLEQSMap_NumberOfLinearMatricesGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(solverColToLinearEquationsMap)) &
      & CALL FlagError("The solver column to linear equations map is not associated.",err,error,*999)
#endif    

    numberOfLinearMatrices=solverColToLinearEquationsMap%numberOfLinearMatrices
      
    EXITS("SolverMappingSColToLEQSMap_NumberOfLinearMatricesGet")
    RETURN
999 numberOfLinearMatrices=0
    ERRORS("SolverMappingSColToLEQSMap_NumberOfLinearMatricesGet",err,error)
    EXITS("SolverMappingSColToLEQSMap_NumberOfLinearMatricesGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingSColToLEQSMap_NumberOfLinearMatricesGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the Jacobian matrix coupling information for a Jacobian matrix in solver column to nonlinear equations map.
  SUBROUTINE SolverMappingSColToNLEQSMap_JacobianCouplingInfoGet(solverColToNonlinearEquationsMap,jacobianMatrixIdx, &
    & jacobianMatrixNumber,columnNumber,couplingCoefficient,err,error,*)

    !Argument variables
    TYPE(SolverColToNonlinearEquationsMapType), POINTER :: solverColToNonlinearEquationsMap !<A pointer to the solver column to nonlinear equations map to get the Jacobian coupling information for
    INTEGER(INTG), INTENT(IN) :: jacobianMatrixIdx !<The Jacobian matrix index to get the coupling information for
    INTEGER(INTG), INTENT(OUT) :: jacobianMatrixNumber !<On exit, the number of the Jacobian matrix in the solver column to nonlinear equations map.
    INTEGER(INTG), INTENT(OUT) :: columnNumber !<On exit, the column number in the Jacobian matrix the solver column is coupled to
    REAL(DP), INTENT(OUT) :: couplingCoefficient !<On exit, the coupling coefficient applied to the solver column to nonlinear matrix column
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif
    
    ENTERS("SolverMappingSColToNLEQSMap_JacobianCouplingInfoGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(solverColToNonlinearEquationsMap)) &
      & CALL FlagError("The solver column to nonlinear equations map is not associated.",err,error,*999)
    IF(jacobianMatrixIdx<1.OR.jacobianMatrixIdx>solverColToNonlinearEquationsMap%numberOfJacobianMatrices) THEN
      localError="The specified Jacobian matrix index of "//TRIM(NumberToVString(jacobianMatrixIdx,"*",err,error))// &
        & " is invalid. The Jacobian matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(solverColToNonlinearEquationsMap%numberOfJacobianMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(solverColToNonlinearEquationsMap%jacobianMatrixNumbers)) &
      & CALL FlagError("The Jacobian matrix numbers array is not allocated for the solver column to nonlinear matrices map.", &
      & err,error,*999)
    IF(.NOT.ALLOCATED(solverColToNonlinearEquationsMap%jacobianColumnNumbers)) &
      & CALL FlagError("The Jacobian column numbers array is not allocated for the solver column to nonlinear matrices map.", &
      & err,error,*999)
    IF(.NOT.ALLOCATED(solverColToNonlinearEquationsMap%couplingCoefficients)) &
      & CALL FlagError("The coupling coefficients array is not allocated for the solver column to nonlinear matrices map.", &
      & err,error,*999)
#endif    

    jacobianMatrixNumber=solverColToNonlinearEquationsMap%jacobianMatrixNumbers(jacobianMatrixIdx)
    columnNumber=solverColToNonlinearEquationsMap%jacobianColumnNumbers(jacobianMatrixIdx)
    couplingCoefficient=solverColToNonlinearEquationsMap%couplingCoefficients(jacobianMatrixIdx)
      
    EXITS("SolverMappingSColToNLEQSMap_JacobianCouplingInfoGet")
    RETURN
999 jacobianMatrixNumber=0
    columnNumber=0
    couplingCoefficient=0.0_DP
    ERRORS("SolverMappingSColToNLEQSMap_JacobianCouplingInfoGet",err,error)
    EXITS("SolverMappingSColToNLEQSMap_JacobianCouplingInfoGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingSColToNLEQSMap_JacobianCouplingInfoGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of Jacobian matrices for a solver column to nonlinear equations map.
  SUBROUTINE SolverMappingSColToNLEQSMap_NumberOfJacobianMatricesGet(solverColToNonlinearEquationsMap,numberOfJacobianMatrices, &
    & err,error,*)

    !Argument variables
    TYPE(SolverColToNonlinearEquationsMapType), POINTER :: solverColToNonlinearEquationsMap !<A pointer to the solver column to nonlinear equations map to get the number of Jacobian matrices for
    INTEGER(INTG), INTENT(OUT) :: numberOfJacobianMatrices !<On exit, the number of Jacobian equations matrices in the solver column to nonlinear equations map.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingSColToNLEQSMap_NumberOfJacobianMatricesGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(solverColToNonlinearEquationsMap)) &
      & CALL FlagError("The solver column to nonlinear equations map is not associated.",err,error,*999)
#endif    

    numberOfJacobianMatrices=solverColToNonlinearEquationsMap%numberOfJacobianMatrices
      
    EXITS("SolverMappingSColToNLEQSMap_NumberOfJacobianMatricesGet")
    RETURN
999 numberOfJacobianMatrices=0
    ERRORS("SolverMappingSColToNLEQSMap_NumberOfJacobianMatricesGet",err,error)
    EXITS("SolverMappingSColToNLEQSMap_NumberOfJacobianMatricesGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingSColToNLEQSMap_NumberOfJacobianMatricesGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns equations information for solver DOF mapped to variable DOFs.
  SUBROUTINE SolverMappingSDOFToVDOFsMap_EquationsInfoGet(solverDOFToVariableDOFsMap,equationDOFIdx,equationType,equationIdx, &
    & err,error,*)

    !Argument variables
    TYPE(SolverDOFToVariableDOFsMapType), POINTER :: solverDOFToVariableDOFsMap !<A pointer to the solver DOF to variable DOFs map to get the equations information for
    INTEGER(INTG), INTENT(IN) :: equationDOFIdx !<The equation DOF index to get the equations information for.
    INTEGER(INTG), INTENT(OUT) :: equationType !<On exit, the type of equation the solver DOF is mapped to. \see SolverMappingRoutines_EquationsTypes
    INTEGER(INTG), INTENT(OUT) :: equationIdx !<On exit, the index of either the equations set or interface condition the solver DOF is mapped to in the solver mapping.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMappingSDOFToVDOFsMap_EquationsInfoGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(solverDOFToVariableDOFsMap)) &
      & CALL FlagError("The solver DOF to variable DOFs map is not associated.",err,error,*999)
    IF(equationDOFIdx<1.OR.equationDOFIdx>solverDOFToVariableDOFsMap%numberOfEquationDOFs) THEN
      localError="The specified equation DOF index of "//TRIM(NumberToVString(equationDOFIdx,"*",err,error))// &
        & " is invalid. The equation DOF index should be >= 1 and <= "// &
        & TRIM(NumberToVString(solverDOFToVariableDOFsMap%numberOfEquationDOFs,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(solverDOFToVariableDOFsMap%equationTypes)) &
      & CALL FlagError("The equations types array is not allocated for the solver DOF to variable DOFs map.",err,error,*999)
    IF(.NOT.ALLOCATED(solverDOFToVariableDOFsMap%equationIndices)) &
      & CALL FlagError("The equations indices array is not allocated for the solver DOF to variable DOFs map.",err,error,*999)    
#endif    

    equationType=solverDOFToVariableDOFsMap%equationTypes(equationDOFIdx)
    equationIdx=solverDOFToVariableDOFsMap%equationIndices(equationDOFIdx)

!!TODO: check validity?
      
    EXITS("SolverMappingSDOFToVDOFsMap_EquationsInfoGet")
    RETURN
999 equationType=0
    equationIdx=0
    ERRORS("SolverMappingSDOFToVDOFsMap_EquationsInfoGet",err,error)
    EXITS("SolverMappingSDOFToVDOFsMap_EquationsInfoGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingSDOFToVDOFsMap_EquationsInfoGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of equations DOFs that the solver DOF is mapped to for the solver DOF to variable DOFs map.
  SUBROUTINE SolverMappingSDOFToVDOFsMap_NumberOfEquationDOFsGet(solverDOFToVariableDOFsMap,numberOfEquationDOFs,err,error,*)

    !Argument variables
    TYPE(SolverDOFToVariableDOFsMapType), POINTER :: solverDOFToVariableDOFsMap !<A pointer to the solver DOF to variable DOFs map to get the variable DOF coupling information for
    INTEGER(INTG), INTENT(OUT) :: numberOfEquationDOFs !<On exit, the number of equation DOFs that the solver DOF is mapped to.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingSDOFToVDOFsMap_NumberOfEquationDOFsGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(solverDOFToVariableDOFsMap)) &
      & CALL FlagError("The solver DOF to variable DOFs map is not associated.",err,error,*999)
#endif    

    numberOfEquationDOFs=solverDOFToVariableDOFsMap%numberOfEquationDOFs
      
    EXITS("SolverMappingSDOFToVDOFsMap_NumberOfEquationDOFsGet")
    RETURN
999 numberOfEquationDOFs=0
    ERRORS("SolverMappingSDOFToVDOFsMap_NumberOfEquationDOFsGet",err,error)
    EXITS("SolverMappingSDOFToVDOFsMap_NumberOfEquationDOFsGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingSDOFToVDOFsMap_NumberOfEquationDOFsGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns variable DOF coupling information in the map for the solver DOF mapped to variable DOFs.
  SUBROUTINE SolverMappingSDOFToVDOFsMap_VariableDOFCouplingGet(solverDOFToVariableDOFsMap,equationDOFIdx,variableDOF, &
    & variableCoefficient,additiveConstant,err,error,*)

    !Argument variables
    TYPE(SolverDOFToVariableDOFsMapType), POINTER :: solverDOFToVariableDOFsMap !<A pointer to the solver DOF to variable DOFs map to get the variable DOF coupling information for
    INTEGER(INTG), INTENT(IN) :: equationDOFIdx !<The equation DOF index to get the variable DOF coupling for.
    INTEGER(INTG), INTENT(OUT) :: variableDOF !<On exit, the variable DOF that the solver DOF is mapped to.
    REAL(DP), INTENT(OUT) :: variableCoefficient !<On exit, the coefficient multiplying the solver DOF for the variable DOF coupling.
    REAL(DP), INTENT(OUT) :: additiveConstant !<On exit, the constant added to the solver DOF for the variable DOF coupling.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMappingSDOFToVDOFsMap_VariableDOFCouplingGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(solverDOFToVariableDOFsMap)) &
      & CALL FlagError("The solver DOF to variable DOFs map is not associated.",err,error,*999)
    IF(equationDOFIdx<1.OR.equationDOFIdx>solverDOFToVariableDOFsMap%numberOfEquationDOFs) THEN
      localError="The specified equations DOF index of "//TRIM(NumberToVString(equationDOFIdx,"*",err,error))// &
        & " is invalid. The equations DOF index should be >= 1 and <= "// &
        & TRIM(NumberToVString(solverDOFToVariableDOFsMap%numberOfEquationDOFs,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(solverDOFToVariableDOFsMap%variableDOF)) &
      & CALL FlagError("The variable DOF array is not allocated for the solver DOF to variable DOFs map.",err,error,*999)
    IF(.NOT.ALLOCATED(solverDOFToVariableDOFsMap%variableCoefficient)) &
      & CALL FlagError("The variable coefficient array is not allocated for the solver DOF to variable DOFs map.",err,error,*999)
    IF(.NOT.ALLOCATED(solverDOFToVariableDOFsMap%additiveConstant)) &
      & CALL FlagError("The additive constant array is not allocated for the solver DOF to variable DOFs map.",err,error,*999)
#endif    

    variableDOF=solverDOFToVariableDOFsMap%variableDOF(equationDOFIdx)
    variableCoefficient=solverDOFToVariableDOFsMap%variableCoefficient(equationDOFIdx)
    additiveConstant=solverDOFToVariableDOFsMap%additiveConstant(equationDOFIdx)

!!TODO: check validity?
      
    EXITS("SolverMappingSDOFToVDOFsMap_VariableDOFCouplingGet")
    RETURN
999 variableDOF=0
    variableCoefficient=0.0_DP
    additiveConstant=0.0_DP
    ERRORS("SolverMappingSDOFToVDOFsMap_VariableDOFCouplingGet",err,error)
    EXITS("SolverMappingSDOFToVDOFsMap_VariableDOFCouplingGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingSDOFToVDOFsMap_VariableDOFCouplingGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the field variable containing the variable DOF for solver DOF mapped to variable DOFs.
  SUBROUTINE SolverMappingSDOFToVDOFsMap_VariableGet(solverDOFToVariableDOFsMap,equationDOFIdx,fieldVariable,err,error,*)

    !Argument variables
    TYPE(SolverDOFToVariableDOFsMapType), POINTER :: solverDOFToVariableDOFsMap !<A pointer to the solver DOF to variable DOFs map to get the equations information for
    INTEGER(INTG), INTENT(IN) :: equationDOFIdx !<The equation DOF index to get the equations information for.
    TYPE(FieldVariableType), POINTER :: fieldVariable !<On exit, a pointer to the field variable containing the DOF the solver DOF is mapped to. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMappingSDOFToVDOFsMap_VariableGet",err,error,*998)

#ifdef WITH_PRECHECKS
    IF(ASSOCIATED(fieldVariable)) CALL FlagError("The field variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverDOFToVariableDOFsMap)) &
      & CALL FlagError("The solver DOF to variable DOFs map is not associated.",err,error,*999)
    IF(equationDOFIdx<1.OR.equationDOFIdx>solverDOFToVariableDOFsMap%numberOfEquationDOFs) THEN
      localError="The specified equations DOF index of "//TRIM(NumberToVString(equationDOFIdx,"*",err,error))// &
        & " is invalid. The equations DOF index should be >= 1 and <= "// &
        & TRIM(NumberToVString(solverDOFToVariableDOFsMap%numberOfEquationDOFs,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(solverDOFToVariableDOFsMap%variable)) &
      & CALL FlagError("The variable array is not allocated for the solver DOF to variable DOFs map.",err,error,*999)
#endif    

    fieldVariable=>solverDOFToVariableDOFsMap%variable(equationDOFIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(fieldVariable)) THEN
      localError="The field variable for equations DOF index "//TRIM(NumberToVString(equationDOFIdx,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("SolverMappingSDOFToVDOFsMap_VariableGet")
    RETURN
999 NULLIFY(fieldVariable)
998 ERRORS("SolverMappingSDOFToVDOFsMap_VariableGet",err,error)
    EXITS("SolverMappingSDOFToVDOFsMap_VariableGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingSDOFToVDOFsMap_VariableGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the column DOFs mapping for a solver matrix to equations map.
  SUBROUTINE SolverMappingSMToEQSMap_ColumnDOFsMappingGet(solverMatrixToEquationsMap,columnDOFsMapping,err,error,*)

    !Argument variables
    TYPE(SolverMatrixToEquationsMapType), POINTER :: solverMatrixToEquationsMap !<A pointer to the solver matrix to equations map to get the column DOFs mapping for
    TYPE(DomainMappingType), POINTER :: columnDOFsMapping !<On exit, a pointer to the specified column DOFs mapping. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingSMToEQSMap_ColumnDOFsMappingGet",err,error,*998)

#ifdef WITH_PRECHECKS
    IF(ASSOCIATED(columnDOFsMapping)) CALL FlagError("Column DOFs mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMatrixToEquationsMap)) &
      & CALL FlagError("Solver matrix to equations map is not associated.",err,error,*999)
#endif    

    columnDOFsMapping=>solverMatrixToEquationsMap%columnDOFsMapping

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(columnDOFsMapping)) &
      & CALL FlagError("The column DOFs mapping for the specified solver matrix to equations map is not associated. ", &
      & err,error,*999)
#endif    
      
    EXITS("SolverMappingSMToEQSMap_ColumnDOFsMappingGet")
    RETURN
999 NULLIFY(columnDOFsMapping)
998 ERRORS("SolverMappingSMToEQSMap_ColumnDOFsMappingGet",err,error)
    EXITS("SolverMappingSMToEQSMap_ColumnDOFsMappingGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingSMToEQSMap_ColumnDOFsMappingGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of columns for a solver matrix to equations map.
  SUBROUTINE SolverMappingSMToEQSMap_NumberOfColumnsGet(solverMatrixToEquationsMap,numberOfColumns,err,error,*)

    !Argument variables
    TYPE(SolverMatrixToEquationsMapType), POINTER :: solverMatrixToEquationsMap !<A pointer to the solver matrix to equations map to get the number of columns for
    INTEGER(INTG), INTENT(OUT) :: numberOfColumns !<On exit, the number of columns in the solver matrix to equations map
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingSMToEQSMap_NumberOfColumnsGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(solverMatrixToEquationsMap)) &
      & CALL FlagError("Solver matrix to equations map is not associated.",err,error,*999)
#endif    

    numberOfColumns=solverMatrixToEquationsMap%numberOfColumns
     
    EXITS("SolverMappingSMToEQSMap_NumberOfColumnsGet")
    RETURN
999 ERRORS("SolverMappingSMToEQSMap_NumberOfColumnsGet",err,error)
    EXITS("SolverMappingSMToEQSMap_NumberOfColumnsGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingSMToEQSMap_NumberOfColumnsGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of DOFs for a solver matrix to equations map.
  SUBROUTINE SolverMappingSMToEQSMap_NumberOfDOFsGet(solverMatrixToEquationsMap,numberOfDOFs,err,error,*)

    !Argument variables
    TYPE(SolverMatrixToEquationsMapType), POINTER :: solverMatrixToEquationsMap !<A pointer to the solver matrix to equations map to get the number of DOFs for
    INTEGER(INTG), INTENT(OUT) :: numberOfDOFs !<On exit, the number of DOFs in the solver matrix to equations map
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingSMToEQSMap_NumberOfDOFsGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(solverMatrixToEquationsMap)) &
      & CALL FlagError("Solver matrix to equations map is not associated.",err,error,*999)
#endif    

    numberOfDOFs=solverMatrixToEquationsMap%numberOfDOFs
     
    EXITS("SolverMappingSMToEQSMap_NumberOfDOFsGet")
    RETURN
999 ERRORS("SolverMappingSMToEQSMap_NumberOfDOFsGet",err,error)
    EXITS("SolverMappingSMToEQSMap_NumberOfDOFsGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingSMToEQSMap_NumberOfDOFsGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of equations sets for a solver matrix to equations map.
  SUBROUTINE SolverMappingSMToEQSMap_NumberOfEquationsSetsGet(solverMatrixToEquationsMap,numberOfEquationsSets,err,error,*)

    !Argument variables
    TYPE(SolverMatrixToEquationsMapType), POINTER :: solverMatrixToEquationsMap !<A pointer to the solver matrix to equations map to get the number of equations sets for
    INTEGER(INTG), INTENT(OUT) :: numberOfEquationsSets !<On exit, the number of equations sets in the solver matrix to equations map
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingSMToEQSMap_NumberOfEquationsSetsGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(solverMatrixToEquationsMap)) &
      & CALL FlagError("Solver matrix to equations map is not associated.",err,error,*999)
#endif    

    numberOfEquationsSets=solverMatrixToEquationsMap%numberOfEquationsSets
     
    EXITS("SolverMappingSMToEQSMap_NumberOfEquationsSetsGet")
    RETURN
999 ERRORS("SolverMappingSMToEQSMap_NumberOfEquationsSetsGet",err,error)
    EXITS("SolverMappingSMToEQSMap_NumberOfEquationsSetsGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingSMToEQSMap_NumberOfEquationsSetsGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of global DOFs for a solver matrix to equations map.
  SUBROUTINE SolverMappingSMToEQSMap_NumberOfGlobalDOFsGet(solverMatrixToEquationsMap,numberOfGlobalDOFs,err,error,*)

    !Argument variables
    TYPE(SolverMatrixToEquationsMapType), POINTER :: solverMatrixToEquationsMap !<A pointer to the solver matrix to equations map to get the number of global DOFs for
    INTEGER(INTG), INTENT(OUT) :: numberOfGlobalDOFs !<On exit, the number of global DOFs in the solver matrix to equations map
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingSMToEQSMap_NumberOfGlobalDOFsGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(solverMatrixToEquationsMap)) &
      & CALL FlagError("Solver matrix to equations map is not associated.",err,error,*999)
#endif    

    numberOfGlobalDOFs=solverMatrixToEquationsMap%numberOfGlobalDOFs
     
    EXITS("SolverMappingSMToEQSMap_NumberOfGlobalDOFsGet")
    RETURN
999 ERRORS("SolverMappingSMToEQSMap_NumberOfGlobalDOFsGet",err,error)
    EXITS("SolverMappingSMToEQSMap_NumberOfGlobalDOFsGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingSMToEQSMap_NumberOfGlobalDOFsGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of interface conditions for a solver matrix to equations map.
  SUBROUTINE SolverMappingSMToEQSMap_NumberOfInterfaceConditionsGet(solverMatrixToEquationsMap,numberOfInterfaceConditions, &
    & err,error,*)

    !Argument variables
    TYPE(SolverMatrixToEquationsMapType), POINTER :: solverMatrixToEquationsMap !<A pointer to the solver matrix to equations map to get the number of interface conditions for
    INTEGER(INTG), INTENT(OUT) :: numberOfInterfaceConditions !<On exit, the number of interface conditions in the solver matrix to equations map
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingSMToEQSMap_NumberOfInterfaceConditionsGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(solverMatrixToEquationsMap)) &
      & CALL FlagError("Solver matrix to equations map is not associated.",err,error,*999)
#endif    

    numberOfInterfaceConditions=solverMatrixToEquationsMap%numberOfInterfaceConditions
     
    EXITS("SolverMappingSMToEQSMap_NumberOfInterfaceConditionsGet")
    RETURN
999 ERRORS("SolverMappingSMToEQSMap_NumberOfInterfaceConditionsGet",err,error)
    EXITS("SolverMappingSMToEQSMap_NumberOfInterfaceConditionsGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingSMToEQSMap_NumberOfInterfaceConditionsGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the solver DOF to variable DOFs map for a solver matrix to equations map.
  SUBROUTINE SolverMappingSMToEQSMap_SolverDOFToVariableDOFsMapGet(solverMatrixToEquationsMap,solverDOFIdx, &
    & solverDOFToVariableDOFsMap,err,error,*)

    !Argument variables
    TYPE(SolverMatrixToEquationsMapType), POINTER :: solverMatrixToEquationsMap !<A pointer to the solver matrix to equations map to get the solver DOF to variable DOFs map for
    INTEGER(INTG), INTENT(IN) :: solverDOFIdx !<The solver DOF index to get the solver DOF to variables DOF map for
    TYPE(SolverDOFToVariableDOFsMapType), POINTER :: solverDOFToVariableDOFsMap !<On exit, a pointer to the solver DOF to varible DOFs map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMappingSMToEQSMap_SolverDOFToVariableDOFsMapGet",err,error,*998)

#ifdef WITH_PRECHECKS
    IF(ASSOCIATED(solverDOFToVariableDOFsMap)) &
      & CALL FlagError("Solver DOF to variable DOFs map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMatrixToEquationsMap)) &
      & CALL FlagError("Solver matrix to equations map is not associated.",err,error,*999)
    IF(solverDOFIdx<1.OR.solverDOFIdx>solverMatrixToEquationsMap%numberOfDOFs) THEN
      localError="The specified solver DOF index of "//TRIM(NumberToVString(solverDOFIdx,"*",err,error))// &
        & " is invalid. The solver DOF index should be >= 1 and <= "// &
        & TRIM(NumberToVString(solverMatrixToEquationsMap%numberOfDOFs,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(solverMatrixToEquationsMap%solverDOFToVariableDOFsMaps)) &
      & CALL FlagError("The solver DOF to variable DOFs maps is not allocated for the solver matrix to equations map.", &
      & err,error,*999)    
#endif    

    solverDOFToVariableDOFsMap=>solverMatrixToEquationsMap%solverDOFToVariableDOFsMaps(solverDOFIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(solverDOFToVariableDOFsMap)) THEN
      localError="The solver DOF to variable DOFs map is not associated for solver DOF index "// &
        & TRIM(NumberToVString(solverDOFIdx,"*",err,error))//" of the solver matrix to equations map."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("SolverMappingSMToEQSMap_SolverDOFToVariableDOFsMapGet")
    RETURN
999 NULLIFY(solverDOFToVariableDOFsMap)
998 ERRORS("SolverMappingSMToEQSMap_SolverDOFToVariableDOFsMapGet",err,error)
    EXITS("SolverMappingSMToEQSMap_SolverDOFToVariableDOFsMapGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingSMToEQSMap_SolverDOFToVariableDOFsMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the solver mapping for a solver matrix to equations map.
  SUBROUTINE SolverMappingSMToEQSMap_SolverMappingGet(solverMatrixToEquationsMap,solverMapping,err,error,*)

    !Argument variables
    TYPE(SolverMatrixToEquationsMapType), POINTER :: solverMatrixToEquationsMap !<A pointer to the solver matrix to equations map to get the solver mapping for
    TYPE(SolverMappingType), POINTER :: solverMapping !<On exit, a pointer to the specified solver mapping. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingSMToEQSMap_SolverMappingGet",err,error,*998)

#ifdef WITH_PRECHECKS
    IF(ASSOCIATED(solverMapping)) CALL FlagError("Solver mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMatrixToEquationsMap)) &
      & CALL FlagError("Solver matrix to equations map is not associated.",err,error,*999)
#endif    

    solverMapping=>solverMatrixToEquationsMap%solverMapping

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(solverMapping)) &
      & CALL FlagError("The solver mapping for the specified solver matrix to equations map is not associated. ", &
      & err,error,*999)
#endif    
      
    EXITS("SolverMappingSMToEQSMap_SolverMappingGet")
    RETURN
999 NULLIFY(solverMapping)
998 ERRORS("SolverMappingSMToEQSMap_SolverMappingGet",err,error)
    EXITS("SolverMappingSMToEQSMap_SolverMappingGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingSMToEQSMap_SolverMappingGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the solver matrix for a solver matrix to equations map.
  SUBROUTINE SolverMappingSMToEQSMap_SolverMatrixGet(solverMatrixToEquationsMap,solverMatrix,err,error,*)

    !Argument variables
    TYPE(SolverMatrixToEquationsMapType), POINTER :: solverMatrixToEquationsMap !<A pointer to the solver matrix to equations map to get the solver matrix for
    TYPE(SolverMatrixType), POINTER :: solverMatrix !<On exit, a pointer to the specified solver matrix. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingSMToEQSMap_SolverMatrixGet",err,error,*998)

#ifdef WITH_PRECHECKS
    IF(ASSOCIATED(solverMatrix)) CALL FlagError("Solver matrix is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMatrixToEquationsMap)) &
      & CALL FlagError("Solver matrix to equations map is not associated.",err,error,*999)
#endif    

    solverMatrix=>solverMatrixToEquationsMap%solverMatrix

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(solverMatrix)) &
      & CALL FlagError("The solver matrix for the specified solver matrix to equations map is not associated. ", &
      & err,error,*999)
#endif    
      
    EXITS("SolverMappingSMToEQSMap_SolverMatrixGet")
    RETURN
999 NULLIFY(solverMatrix)
998 ERRORS("SolverMappingSMToEQSMap_SolverMatrixGet",err,error)
    EXITS("SolverMappingSMToEQSMap_SolverMatrixGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingSMToEQSMap_SolverMatrixGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the solver matrix to equations set map for an equations set index in a solver matrix to equations map.
  SUBROUTINE SolverMappingSMToEQSMap_SolverMatrixToEquationsSetMapGet(solverMatrixToEquationsMap,equationsSetIdx, &
    & solverMatrixToEquationsSetMap,err,error,*)

    !Argument variables
    TYPE(SolverMatrixToEquationsMapType), POINTER :: solverMatrixToEquationsMap !<A pointer to the solver matrix to equations map to get the solver matrix to equations set map for
    INTEGER(INTG), INTENT(IN) :: equationsSetIdx !<The equations set index to get the solver matrix to equations set map for
    TYPE(SolverMatrixToEquationsSetMapType), POINTER :: solverMatrixToEquationsSetMap !<On exit, a pointer to the specified solver matrix to equations set map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMappingSMToEQSMap_SolverMatrixToEquationsSetMapGet",err,error,*998)

#ifdef WITH_PRECHECKS
    IF(ASSOCIATED(solverMatrixToEquationsSetMap)) &
      & CALL FlagError("Solver matrix to equations set map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMatrixToEquationsMap)) &
      & CALL FlagError("Solver matrix to equations map is not associated.",err,error,*999)
    IF(equationsSetIdx<1.OR.equationsSetIdx>solverMatrixToEquationsMap%numberOfEquationsSets) THEN
      localError="The specified equations set index of "//TRIM(NumberToVString(equationsSetIdx,"*",err,error))// &
        & " is invalid. The equations set index should be >= 1 and <= "// &
        & TRIM(NumberToVString(solverMatrixToEquationsMap%numberOfEquationsSets,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(solverMatrixToEquationsMap%solverMatrixToEquationsSetMaps)) &
      & CALL FlagError("The solver matrix to equations sets maps is not allocated for the solver matrix to equations map.", &
      & err,error,*999)
#endif    
    
    solverMatrixToEquationsSetMap=>solverMatrixToEquationsMap%solverMatrixToEquationsSetMaps(equationsSetIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(solverMatrixToEquationsSetMap)) THEN
      localError="The solver matrix to equation set map is not associated for equations set index "// &
        & TRIM(NumberToVString(equationsSetIdx,"*",err,error))//" of the solver matrix to equations map."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("SolverMappingSMToEQSMap_SolverMatrixToEquationsSetMapGet")
    RETURN
999 NULLIFY(solverMatrixToEquationsSetMap)
998 ERRORS("SolverMappingSMToEQSMap_SolverMatrixToEquationsSetMapGet",err,error)
    EXITS("SolverMappingSMToEQSMap_SolverMatrixToEquationsSetMapGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingSMToEQSMap_SolverMatrixToEquationsSetMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the solver matrix to interface condition map for an equations set index in a solver matrix to equations map.
  SUBROUTINE SolverMappingSMToEQSMap_SolverMatrixToInterfaceConditionMapGet(solverMatrixToEquationsMap,interfaceConditionIdx, &
    & solverMatrixToInterfaceConditionMap,err,error,*)

    !Argument variables
    TYPE(SolverMatrixToEquationsMapType), POINTER :: solverMatrixToEquationsMap !<A pointer to the solver matrix to equations map to get the solver matrix to interface condition map for
    INTEGER(INTG), INTENT(IN) :: interfaceConditionIdx !<The interface condition index to get the solver matrix to interface condition map for
    TYPE(SolverMatrixToInterfaceConditionMapType), POINTER :: solverMatrixToInterfaceConditionMap !<On exit, a pointer to the specified solver matrix to interface condition map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMappingSMToEQSMap_SolverMatrixToInterfaceConditionMapGet",err,error,*998)

#ifdef WITH_PRECHECKS
    IF(ASSOCIATED(solverMatrixToInterfaceConditionMap)) &
      & CALL FlagError("Solver matrix to equations set map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMatrixToEquationsMap)) &
      & CALL FlagError("Solver matrix to equations map is not associated.",err,error,*999)
    IF(interfaceConditionIdx<1.OR.interfaceConditionIdx>solverMatrixToEquationsMap%numberOfInterfaceConditions) THEN
      localError="The specified interface condition index of "//TRIM(NumberToVString(interfaceConditionIdx,"*",err,error))// &
        & " is invalid. The interface condition index should be >= 1 and <= "// &
        & TRIM(NumberToVString(solverMatrixToEquationsMap%numberOfInterfaceConditions,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(solverMatrixToEquationsMap%solverMatrixToInterfaceConditionMaps)) &
      & CALL FlagError("The solver matrix to interface conditions maps is not allocated for the solver matrix to equations map.", &
      & err,error,*999)
#endif    
    
    solverMatrixToInterfaceConditionMap=>solverMatrixToEquationsMap%solverMatrixToInterfaceConditionMaps(interfaceConditionIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(solverMatrixToInterfaceConditionMap)) THEN
      localError="The solver matrix to interface condition map is not associated for interface condition index "// &
        & TRIM(NumberToVString(interfaceConditionIdx,"*",err,error))//" of the solver matrix to equations map."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("SolverMappingSMToEQSMap_SolverMatrixToInterfaceConditionMapGet")
    RETURN
999 NULLIFY(solverMatrixToInterfaceConditionMap)
998 ERRORS("SolverMappingSMToEQSMap_SolverMatrixToInterfaceConditionMapGet",err,error)
    EXITS("SolverMappingSMToEQSMap_SolverMatrixToInterfaceConditionMapGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingSMToEQSMap_SolverMatrixToInterfaceConditionMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the total number of DOFs for a solver matrix to equations map.
  SUBROUTINE SolverMappingSMToEQSMap_TotalNumberOfDOFsGet(solverMatrixToEquationsMap,totalNumberOfDOFs,err,error,*)

    !Argument variables
    TYPE(SolverMatrixToEquationsMapType), POINTER :: solverMatrixToEquationsMap !<A pointer to the solver matrix to equations map to get the number of DOFs for
    INTEGER(INTG), INTENT(OUT) :: totalNumberOfDOFs !<On exit, the total number of DOFs in the solver matrix to equations map
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingSMToEQSMap_TotalNumberOfDOFsGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(solverMatrixToEquationsMap)) &
      & CALL FlagError("Solver matrix to equations map is not associated.",err,error,*999)
#endif    

    totalNumberOfDOFs=solverMatrixToEquationsMap%totalNumberOfDOFs
     
    EXITS("SolverMappingSMToEQSMap_TotalNumberOfDOFsGet")
    RETURN
999 ERRORS("SolverMappingSMToEQSMap_TotalNumberOfDOFsGet",err,error)
    EXITS("SolverMappingSMToEQSMap_TotalNumberOfDOFsGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingSMToEQSMap_TotalNumberOfDOFsGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to equations variables list for a solver matrix to equations map.
  SUBROUTINE SolverMappingSMToEQSMap_VariablesListGet(solverMatrixToEquationsMap,variablesList,err,error,*)

    !Argument variables
    TYPE(SolverMatrixToEquationsMapType), POINTER :: solverMatrixToEquationsMap !<A pointer to the solver matrix to equations map to get the variables list in the solver matrix to equations map for
    TYPE(SolverMappingVariablesType), POINTER :: variablesList !<On exit, a pointer to variables in the solver matrix to equations map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingSMToEQSMap_VariablesListGet",err,error,*998)

#ifdef WITH_PRECHECKS
    IF(ASSOCIATED(variablesList)) CALL FlagError("Solver variables list is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMatrixToEquationsMap)) &
      & CALL FlagError("Solver matrix to equations map is not associated.",err,error,*999)
#endif    
    
    variablesList=>solverMatrixToEquationsMap%variablesList

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(variablesList)) &
      & CALL FlagError("The solver variables list is not associated for the solver matrix to equations map.",err,error,*999)
#endif    
      
    EXITS("SolverMappingSMToEQSMap_VariablesListGet")
    RETURN
999 NULLIFY(variablesList)
998 ERRORS("SolverMappingSMToEQSMap_VariablesListGet",err,error)
    EXITS("SolverMappingSMToEQSMap_VariablesListGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingSMToEQSMap_VariablesListGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the equations set row coupling information for a solver row to equations set row.
  SUBROUTINE SolverMappingSRowToEQSMap_EquationsSetCouplingInfoGet(solverRowToEquationsMap,coupledRowIdx,equationsSetIdx, &
    & equationsSetRow,couplingCoefficient,err,error,*)

    !Argument variables
    TYPE(SolverRowToEquationsMapType), POINTER :: solverRowToEquationsMap !<A pointer to the solver row to equations map to get the equations set row coupling information for
    INTEGER(INTG), INTENT(IN) :: coupledRowIdx !<The coupled row index to get the coupling information for
    INTEGER(INTG), INTENT(OUT) :: equationsSetIdx !<On exit, the equations set index of the coupled solver row
    INTEGER(INTG), INTENT(OUT) :: equationsSetRow !<On exit, the equations set row of the coupled solver row
    REAL(DP), INTENT(OUT) :: couplingCoefficient !<On exit, the equations set row coupling coefficient of the coupled solver row
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMappingSRowToEQSMap_EquationsSetCouplingInfoGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(solverRowToEquationsMap)) &
      & CALL FlagError("Solver row to equations map is not associated.",err,error,*999)
    !Check that we have an equations set coupling
    IF(solverRowToEquationsMap%interfaceConditionIndex/=0) &
      & CALL FlagError("The solver row to equations map is coupled to interface equations and not equations set equations.", &
      & err,error,*999)
    IF(coupledRowIdx<1.OR.coupledRowIdx>solverRowToEquationsMap%numberOfEquationsSetRows) THEN
      localError="The specified equations set row index of "//TRIM(NumberToVString(coupledRowIdx,"*",err,error))// &
        & " is invalid. The equations set row index should be >= 1 and <= "// &
        & TRIM(NumberToVString(solverRowToEquationsMap%numberOfEquationsSetRows,"*",err,error))// &
        & " for the solver row to equations map."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(solverRowToEquationsMap%equationsIndex)) &
      & CALL FlagError("The equations index array is not allocated for the solver row to equations map.", &
      & err,error,*999)
    IF(.NOT.ALLOCATED(solverRowToEquationsMap%rowColNumber)) &
      & CALL FlagError("The row column number array is not allocated for the solver row to equations map.", &
      & err,error,*999)
    IF(.NOT.ALLOCATED(solverRowToEquationsMap%couplingCoefficients)) &
      & CALL FlagError("The coupling coefficients array is not allocated for the solver row to equations map.", &
      & err,error,*999)
#endif    

    equationsSetIdx=solverRowToEquationsMap%equationsIndex(coupledRowIdx)
    equationsSetRow=solverRowToEquationsMap%rowColNumber(coupledRowIdx)
    couplingCoefficient=solverRowToEquationsMap%couplingCoefficients(coupledRowIdx)
     
    EXITS("SolverMappingSRowToEQSMap_EquationsSetCouplingInfoGet")
    RETURN
999 equationsSetIdx=0
    equationsSetRow=0
    couplingCoefficient=0.0_DP
    ERRORS("SolverMappingSRowToEQSMap_EquationsSetCouplingInfoGet",err,error)
    EXITS("SolverMappingSRowToEQSMap_EquationsSetCouplngInfoGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingSRowToEQSMap_EquationsSetCouplingInfoGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the interface condition column coupling information for a solver row to interface condition column.
  SUBROUTINE SolverMappingSRowToEQSMap_InterfaceConditionCouplingInfoGet(solverRowToEquationsMap,interfaceConditionColumn, &
    & couplingCoefficient,err,error,*)

    !Argument variables
    TYPE(SolverRowToEquationsMapType), POINTER :: solverRowToEquationsMap !<A pointer to the solver row to equation map to get the interface condition column coupling information for
    INTEGER(INTG), INTENT(OUT) :: interfaceConditionColumn !<On exit, the interface condition column of the coupled solver row
    REAL(DP), INTENT(OUT) :: couplingCoefficient !<On exit, the interface condition coupling coefficient of the coupled solver row
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingSRowToEQSMap_InterfaceConditionCouplingInfoGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(solverRowToEquationsMap)) &
      & CALL FlagError("Solver row to equations map is not associated.",err,error,*999)
    !Check that we have an interface condition coupling
    IF(solverRowToEquationsMap%numberOfEquationsSetRows/=0) &
      & CALL FlagError("The solver row to equations map is coupled to equations set equations and not "// &
      & "interface conditions equations.",err,error,*999)
    IF(.NOT.ALLOCATED(solverRowToEquationsMap%rowColNumber)) &
      & CALL FlagError("The row column number array is not allocated for the solver row to equations map.", &
      & err,error,*999)
    IF(.NOT.ALLOCATED(solverRowToEquationsMap%couplingCoefficients)) &
      & CALL FlagError("The coupling coefficients array is not allocated for the solver row to equations map.", &
      & err,error,*999)
#endif    

    interfaceConditionColumn=solverRowToEquationsMap%rowColNumber(1)
    couplingCoefficient=solverRowToEquationsMap%couplingCoefficients(1)
     
    EXITS("SolverMappingSRowToEQSMap_InterfaceConditionCouplingInfoGet")
    RETURN
999 interfaceConditionColumn=0
    couplingCoefficient=0.0_DP
    ERRORS("SolverMappingSRowToEQSMap_InterfaceConditionCouplingInfoGet",err,error)
    EXITS("SolverMappingSRowToEQSMap_InterfaceConditionCouplingInfoGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingSRowToEQSMap_InterfaceConditionCouplingInfoGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the interface condition index for a solver row to interface condition column.
  SUBROUTINE SolverMappingSRowToEQSMap_InterfaceConditionIndexGet(solverRowToEquationsMap,interfaceConditionIndex,err,error,*)

    !Argument variables
    TYPE(SolverRowToEquationsMapType), POINTER :: solverRowToEquationsMap !<A pointer to the solver row to equation map to get the interface condition column coupling information for
    INTEGER(INTG), INTENT(OUT) :: interfaceConditionIndex !<On exit, the interface condition index of the coupled solver row
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingSRowToEQSMap_InterfaceConditionIndexGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(solverRowToEquationsMap)) &
      & CALL FlagError("Solver row to equations map is not associated.",err,error,*999)
#endif    

    interfaceConditionIndex=solverRowToEquationsMap%interfaceConditionIndex
     
    EXITS("SolverMappingSRowToEQSsMap_InterfaceConditionIndexGet")
    RETURN
999 interfaceConditionIndex=0
    ERRORS("SolverMappingSRowToEQSMap_InterfaceConditionIndexGet",err,error)
    EXITS("SolverMappingSRowToEQSMap_InterfaceConditionIndexGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingSRowToEQSMap_InterfaceConditionIndexGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of equations set rows for a solver row to equations map.
  SUBROUTINE SolverMappingSRowToEQSMap_NumberOfEquationsSetRowsGet(solverRowToEquationsMap,numberOfEquationsSetRows,err,error,*)

    !Argument variables
    TYPE(SolverRowToEquationsMapType), POINTER :: solverRowToEquationsMap !<A pointer to the solver row to equation map to get the number of equaitons set rows for
    INTEGER(INTG), INTENT(OUT) :: numberOfEquationsSetRows !<On exit, the number of equations set rows of the coupled solver row
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingSRowToEQSMap_NumberOfEquationsSetRowsGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(solverRowToEquationsMap)) &
      & CALL FlagError("Solver row to equations map is not associated.",err,error,*999)
#endif    

    numberOfEquationsSetRows=solverRowToEquationsMap%numberOfEquationsSetRows
     
    EXITS("SolverMappingSRowToEQSMap_NumberOfEquationsSetRowsGet")
    RETURN
999 numberOfEquationsSetRows=0
    ERRORS("SolverMappingSRowToEQSMap_NumberOfEquationsSetRowsGet",err,error)
    EXITS("SolverMappingSRowToEQSMap_NumberOfEquationsSetRowsGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingSRowToEQSMap_NumberOfEquationsSetRowsGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the equations for a solver matrix to equations set map.
  SUBROUTINE SolverMappingSMToESMap_EquationsGet(solverMatrixToEquationsSetMap,equations,err,error,*)

    !Argument variables
    TYPE(SolverMatrixToEquationsSetMapType), POINTER :: solverMatrixToEquationsSetMap !<A pointer to the solver matrix to equations set map to get the equations for
    TYPE(EquationsType), POINTER :: equations !<On exit, a pointer to the equations fo the solver matrix to equations set map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingSMToESMap_EquationsGet",err,error,*998)

#ifdef WITH_PRECHECKS
    IF(ASSOCIATED(equations)) CALL FlagError("Equations is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMatrixToEquationsSetMap)) &
      & CALL FlagError("The solver matrix to equations set map is not associated.",err,error,*999)
#endif    
    
    equations=>solverMatrixToEquationsSetMap%equations

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(equations)) &
      & CALL FlagError("The equations is not associated for the solver matrix to equations set map.",err,error,*999)
#endif    
      
    EXITS("SolverMappingSMToESMap_EquationsGet")
    RETURN
999 NULLIFY(equations)
998 ERRORS("SolverMappingSMToESMap_EquationsGet",err,error)
    EXITS("SolverMappingSMToESMap_EquationsGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingSMToESMap_EquationsGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the have equations flags for a solver matrix to equations set map.
  SUBROUTINE SolverMappingSMToESMap_HaveEquationsGet(solverMatrixToEquationsSetMap,haveDynamic,haveLinear,haveNonlinear, &
    & err,error,*)

    !Argument variables
    TYPE(SolverMatrixToEquationsSetMapType), POINTER :: solverMatrixToEquationsSetMap !<A pointer to the solver matrix to equations set map to get the have equations for
    LOGICAL, INTENT(OUT) :: haveDynamic !<On exit, the have dynamic equations flag in the solver matrix to equations set map.
    LOGICAL, INTENT(OUT) :: haveLinear !<On exit, the have linear equations flag in the solver matrix to equations set map.
    LOGICAL, INTENT(OUT) :: haveNonlinear !<On exit, the have nonlinear equations flag in the solver matrix to equations set map.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingSMToESMap_HaveEquationsGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(solverMatrixToEquationsSetMap)) &
      & CALL FlagError("The solver matrix to equations set map is not associated.",err,error,*999)
#endif    
    
    haveDynamic=solverMatrixToEquationsSetMap%haveDynamic
    haveLinear=solverMatrixToEquationsSetMap%haveLinear
    haveNonlinear=solverMatrixToEquationsSetMap%haveNonlinear

    EXITS("SolverMappingSMToESMap_HaveEquationsGet")
    RETURN
999 ERRORS("SolverMappingSMToESMap_HaveEquationsGet",err,error)
    EXITS("SolverMappingSMToESMap_HaveEquationsGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingSMToESMap_HaveEquationsGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the solver column to dynamic equations map for a column of a solver matrix to equations set map.
  SUBROUTINE SolverMappingSMToESMap_SolverColToDynamicEquationsMapGet(solverMatrixToEquationsSetMap,columnIdx, &
    & solverColToDynamicEquationsMap,err,error,*)

    !Argument variables
    TYPE(SolverMatrixToEquationsSetMapType), POINTER :: solverMatrixToEquationsSetMap !<A pointer to the solver matrix to equations set map to get the solver column to dynamic equations map for
    INTEGER(INTG), INTENT(IN) :: columnIdx !<The column index to get the solver column to dynamic equations map for
    TYPE(SolverColToDynamicEquationsMapType), POINTER :: solverColToDynamicEquationsMap !<On exit, a pointer to the solver column to dynamic equations map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMappingSMToESMap_SolverColToDynamicEquationsMapGet",err,error,*998)

#ifdef WITH_PRECHECKS
    IF(ASSOCIATED(solverColToDynamicEquationsMap)) &
      & CALL FlagError("The solver column to dynamic equations map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMatrixToEquationsSetMap)) &
      & CALL FlagError("The solver matrix to equations set map is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(solverMatrixToEquationsSetMap%solverColToDynamicEquationsMaps)) &
      & CALL FlagError("The solver column to dynamic equations map is not allocated for the solver matrix to equations set map.", &
      & err,error,*999)
    IF(columnIdx<1.OR.columnIdx>SIZE(solverMatrixToEquationsSetMap%solverColToDynamicEquationsMaps,1)) THEN
      localError="The specified column index of "//TRIM(NumberToVString(columnIdx,"*",err,error))// &
        & " is invalid. The column index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(solverMatrixToEquationsSetMap%solverColToDynamicEquationsMaps,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF    
#endif    

    solverColToDynamicEquationsMap=>solverMatrixToEquationsSetMap%solverColToDynamicEquationsMaps(columnIdx)%ptr

#ifdef WITH_POSTCHECKS
    IF(.NOT.ASSOCIATED(solverColToDynamicEquationsMap)) THEN
      localError="The solver column to dynamic equations map is not associated for column index "// &
        & TRIM(NumberToVString(columnIdx,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif

    EXITS("SolverMappingSMToESMap_SolverColToDynamicEquationsMapGet")
    RETURN
999 NULLIFY(solverColToDynamicEquationsMap)
998 ERRORS("SolverMappingSMToESMap_SolverColToDynamicEquationsMapGet",err,error)
    EXITS("SolverMappingSMToESMap_SolverColToDynamicEquationsMapGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingSMToESMap_SolverColToDynamicEquationsMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the solver column to linear equations map for a column of a solver matrix to equations set map.
  SUBROUTINE SolverMappingSMToESMap_SolverColToLinearEquationsMapGet(solverMatrixToEquationsSetMap,columnIdx, &
    & solverColToLinearEquationsMap,err,error,*)

    !Argument variables
    TYPE(SolverMatrixToEquationsSetMapType), POINTER :: solverMatrixToEquationsSetMap !<A pointer to the solver matrix to equations set map to get the solver column to linear equations map for
    INTEGER(INTG), INTENT(IN) :: columnIdx !<The column index to get the solver column to linear equations map for
    TYPE(SolverColToLinearEquationsMapType), POINTER :: solverColToLinearEquationsMap !<On exit, a pointer to the solver column to linear equations map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMappingSMToESMap_SolverColToLinearEquationsMapGet",err,error,*998)

#ifdef WITH_PRECHECKS
    IF(ASSOCIATED(solverColToLinearEquationsMap)) &
      & CALL FlagError("The solver column to linear equations map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMatrixToEquationsSetMap)) &
      & CALL FlagError("The solver matrix to equations set map is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(solverMatrixToEquationsSetMap%solverColToLinearEquationsMaps)) &
      & CALL FlagError("The solver column to linear equations map is not allocated for the solver matrix to equations set map.", &
      & err,error,*999)
    IF(columnIdx<1.OR.columnIdx>SIZE(solverMatrixToEquationsSetMap%solverColToLinearEquationsMaps,1)) THEN
      localError="The specified column index of "//TRIM(NumberToVString(columnIdx,"*",err,error))// &
        & " is invalid. The column index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(solverMatrixToEquationsSetMap%solverColToLinearEquationsMaps,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF    
#endif    

    solverColToLinearEquationsMap=>solverMatrixToEquationsSetMap%solverColToLinearEquationsMaps(columnIdx)%ptr

#ifdef WITH_POSTCHECKS
    IF(.NOT.ASSOCIATED(solverColToLinearEquationsMap)) THEN
      localError="The solver column to linear equations map is not associated for column index "// &
        & TRIM(NumberToVString(columnIdx,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif

    EXITS("SolverMappingSMToESMap_SolverColToLinearEquationsMapGet")
    RETURN
999 NULLIFY(solverColToLinearEquationsMap)
998 ERRORS("SolverMappingSMToESMap_SolverColToLinearEquationsMapGet",err,error)
    EXITS("SolverMappingSMToESMap_SolverColToLinearEquationsMapGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingSMToESMap_SolverColToLinearEquationsMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the solver column to nonlinear equations map for a column of a solver matrix to equations set map.
  SUBROUTINE SolverMappingSMToESMap_SolverColToNonlinearEquationsMapGet(solverMatrixToEquationsSetMap,columnIdx, &
    & solverColToNonlinearEquationsMap,err,error,*)

    !Argument variables
    TYPE(SolverMatrixToEquationsSetMapType), POINTER :: solverMatrixToEquationsSetMap !<A pointer to the solver matrix to equations set map to get the solver column to nonlinear equations map for
    INTEGER(INTG), INTENT(IN) :: columnIdx !<The column index to get the solver column to nonlinear equations map for
    TYPE(SolverColToNonlinearEquationsMapType), POINTER :: solverColToNonlinearEquationsMap !<On exit, a pointer to the solver column to nonlinear equations map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMappingSMToESMap_SolverColToNonlinearEquationsMapGet",err,error,*998)

#ifdef WITH_PRECHECKS
    IF(ASSOCIATED(solverColToNonlinearEquationsMap)) &
      & CALL FlagError("The solver column to nonlinear equations map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMatrixToEquationsSetMap)) &
      & CALL FlagError("The solver matrix to equations set map is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(solverMatrixToEquationsSetMap%solverColToNonlinearEquationsMaps)) &
      & CALL FlagError("The solver column to nonlinear equations map is not allocated for the solver matrix to "// &
      & "equations set map.",err,error,*999)
    IF(columnIdx<1.OR.columnIdx>SIZE(solverMatrixToEquationsSetMap%solverColToNonlinearEquationsMaps,1)) THEN
      localError="The specified column index of "//TRIM(NumberToVString(columnIdx,"*",err,error))// &
        & " is invalid. The column index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(solverMatrixToEquationsSetMap%solverColToNonlinearEquationsMaps,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF    
#endif    

    solverColToNonlinearEquationsMap=>solverMatrixToEquationsSetMap%solverColToNonlinearEquationsMaps(columnIdx)%ptr

#ifdef WITH_POSTCHECKS
    IF(.NOT.ASSOCIATED(solverColToNonlinearEquationsMap)) THEN
      localError="The solver column to nonlinear equations map is not associated for column index "// &
        & TRIM(NumberToVString(columnIdx,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif

    EXITS("SolverMappingSMToESMap_SolverColToNonlinearEquationsMapGet")
    RETURN
999 NULLIFY(solverColToNonlinearEquationsMap)
998 ERRORS("SolverMappingSMToESMap_SolverColToNonlinearEquationsMapGet",err,error)
    EXITS("SolverMappingSMToESMap_SolverColToNonlinearEquationsMapGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingSMToESMap_SolverColToNonlinearEquationsMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the interface equations for a solver matrix to interface condiion map.
  SUBROUTINE SolverMappingSMToICMap_InterfaceEquationsGet(solverMatrixToInterfaceConditionMap,interfaceEquations,err,error,*)

    !Argument variables
    TYPE(SolverMatrixToInterfaceConditionMapType), POINTER :: solverMatrixToInterfaceConditionMap !<A pointer to the solver matrix to interface condition map to get the interface equations for
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations !<On exit, a pointer to the interface equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverMappingSMToICMap_InterfaceEquationsGet",err,error,*998)

#ifdef WITH_PRECHECKS
    IF(ASSOCIATED(interfaceEquations)) CALL FlagError("The interface equations is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMatrixToInterfaceConditionMap)) &
      & CALL FlagError("The solver matrix to interface condition map is not associated.",err,error,*999)
#endif    

    interfaceEquations=>solverMatrixToInterfaceConditionMap%interfaceEquations

#ifdef WITH_POSTCHECKS
    IF(.NOT.ASSOCIATED(interfaceEquations)) &
      CALL FlagError("The interface equations map is not associated for solver matrix to interface condition map.",err,error,*999)
#endif

    EXITS("SolverMappingSMToICMap_InterfaceEquationsGet")
    RETURN
999 NULLIFY(interfaceEquations)
998 ERRORS("SolverMappingSMToICMap_InterfaceEquationsGet",err,error)
    EXITS("SolverMappingSMToICMap_InterfaceEquationsGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingSMToICMap_InterfaceEquationsGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the solver column to interface equations map for a column of a solver matrix to interface condiion map.
  SUBROUTINE SolverMappingSMToICMap_SolverColToInterfaceEquationsMapGet(solverMatrixToInterfaceConditionMap,columnIdx, &
    & solverColToInterfaceEquationsMap,err,error,*)

    !Argument variables
    TYPE(SolverMatrixToInterfaceConditionMapType), POINTER :: solverMatrixToInterfaceConditionMap !<A pointer to the solver matrix to interface condition map to get the solver column to interface equations map for
    INTEGER(INTG), INTENT(IN) :: columnIdx !<The column index to get the solver column to interface equations map for
    TYPE(SolverColToInterfaceEquationsMapType), POINTER :: solverColToInterfaceEquationsMap !<On exit, a pointer to the solver column to interface equations map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMappingSMToICMap_SolverColToInterfaceEquationsMapGet",err,error,*998)

#ifdef WITH_PRECHECKS
    IF(ASSOCIATED(solverColToInterfaceEquationsMap)) &
      & CALL FlagError("The solver column to interface equations map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMatrixToInterfaceConditionMap)) &
      & CALL FlagError("The solver matrix to interface condition map is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(solverMatrixToInterfaceConditionMap%solverColToInterfaceEquationsMaps)) &
      & CALL FlagError("The solver column to interface equations map is not allocated for the solver matrix to "// &
      & "interface condition map.",err,error,*999)
    IF(columnIdx<1.OR.columnIdx>SIZE(solverMatrixToInterfaceConditionMap%solverColToInterfaceEquationsMaps,1)) THEN
      localError="The specified column index of "//TRIM(NumberToVString(columnIdx,"*",err,error))// &
        & " is invalid. The column index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(solverMatrixToInterfaceConditionMap%solverColToInterfaceEquationsMaps,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF    
#endif    

    solverColToInterfaceEquationsMap=>solverMatrixToInterfaceConditionMap%solverColToInterfaceEquationsMaps(columnIdx)%ptr

#ifdef WITH_POSTCHECKS
    IF(.NOT.ASSOCIATED(solverColToInterfaceEquationsMap)) THEN
      localError="The solver column to interface equations map is not associated for column index "// &
        & TRIM(NumberToVString(columnIdx,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif

    EXITS("SolverMappingSMToICMap_SolverColToInterfaceEquationsMapGet")
    RETURN
999 NULLIFY(solverColToInterfaceEquationsMap)
998 ERRORS("SolverMappingSMToICMap_SolverColToInterfaceEquationsMapGet",err,error)
    EXITS("SolverMappingSMToICMap_SolverColToInterfaceEquationsMapGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingSMToICMap_SolverColToInterfaceEquationsMapGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the equation information for a variable in the solver mapping variables list.
  SUBROUTINE SolverMappingVariable_EquationInfoGet(solverMappingVariable,equationIdx,equationType,equationIndex,err,error,*)

    !Argument variables
    TYPE(SolverMappingVariableType), POINTER :: solverMappingVariable !<A pointer to the solver mapping variable to get the equation information for
    INTEGER(INTG), INTENT(IN) :: equationIdx !<The equation index of the solver mapping variable to get the equations information for.
    INTEGER(INTG), INTENT(OUT) :: equationType !<On exit, the equation type (equations set or interface condition) that the solver mapping variable is mapped to. \see SolverMappingRoutines_EquationsTypes
    INTEGER(INTG), INTENT(OUT) :: equationIndex !<On exit, the equation set of interface condition equation index that the solver mapping variable is mapped to.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMappingVariable_EquationInfoGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(solverMappingVariable)) CALL FlagError("Solver mapping variable is not associated.",err,error,*999)
    IF(equationIdx<1.OR.equationIdx>solverMappingVariable%numberOfEquations) THEN
      localError="The specified equations index of "//TRIM(NumberToVString(equationIdx,"*",err,error))// &
        & " is invalid. The equations index should be >= 1 and <= "// &
        & TRIM(NumberToVString(solverMappingVariable%numberOfEquations,"*",err,error))// &
        & " for the solver mapping variable with index "// &
        & TRIM(NumberToVString(solverMappingVariable%variableIdx,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(solverMappingVariable%equationTypes)) THEN
      localError="The equations types array is not allocated for the solver mapping variable with index "// &
        & TRIM(NumberToVString(solverMappingVariable%variableIdx,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(solverMappingVariable%equationIndices)) THEN
      localError="The equations indices array is not allocated for the solver mapping variable with index "// &
        & TRIM(NumberToVString(solverMappingVariable%variableIdx,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    equationType=solverMappingVariable%equationTypes(equationIdx)
    equationIndex=solverMappingVariable%equationIndices(equationIdx)
            
    EXITS("SolverMappingVariable_EquationInfoGet")
    RETURN
999 equationType=0
    equationIndex=0
    ERRORS("SolverMappingVariable_EquationInfoGet",err,error)
    EXITS("SolverMappingVariable_EquationInfoGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingVariable_EquationInfoGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the matrix number for a matrix in an equation in a solver mapping variable.
  SUBROUTINE SolverMappingVariable_EquationMatrixNumberGet(solverMappingVariable,equationsIdx,matrixIdx,matrixNumber,err,error,*)

    !Argument variables
    TYPE(SolverMappingVariableType), POINTER :: solverMappingVariable !<A pointer to the solver mapping variable to get the matrix number for
    INTEGER(INTG), INTENT(IN) :: equationsIdx !<The equations index to get the matrix number in the equationIdx'th equation that the solver mapping variable is mapped to.
    INTEGER(INTG), INTENT(IN) :: matrixIdx !<The matrix index to get the matrix number in the equationsIdx'th equation that the solver mapping variable is mapped to.
    INTEGER(INTG), INTENT(OUT) :: matrixNumber !<On exit, the matrix number for the matrixIdx'th matrix in the equationsIdx'th equation that the solver mapping variable is mapped to.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMappingVariable_EquationMatrixNumberGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(solverMappingVariable)) CALL FlagError("Solver mapping variable is not associated.",err,error,*999)
    IF(equationsIdx<1.OR.equationsIdx>solverMappingVariable%numberOfEquations) THEN
      localError="The specified equations index of "//TRIM(NumberToVString(equationsIdx,"*",err,error))// &
        & " is invalid. The equations index should be >= 1 and <= "// &
        & TRIM(NumberToVString(solverMappingVariable%numberOfEquations,"*",err,error))// &
        & " for the solver mapping variable with variable index "// &
        & TRIM(NumberToVString(solverMappingVariable%variableIdx,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(solverMappingVariable%numberOfMatrices)) THEN
      localError="The number of matrices array is not allocated for the solver mapping variable with variable index "// &
        & TRIM(NumberToVString(solverMappingVariable%variableIdx,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(matrixIdx<1.OR.matrixIdx>solverMappingVariable%numberOfMatrices(equationsIdx)) THEN
      localError="The specified matrix index of "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
        & " is invalid. The matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(solverMappingVariable%numberOfMatrices(equationsIdx),"*",err,error))// &
        & " for equations index "//TRIM(NumberToVString(equationsIdx,"*",err,error))// &
        & " of the solver mapping variable with variable index "// &
        & TRIM(NumberToVString(solverMappingVariable%variableIdx,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(solverMappingVariable%matrixNumbers)) THEN
      localError="The matrix numbers array is not allocated for the solver matrix variable with variable index "// &
        & TRIM(NumberToVString(solverMappingVariable%variableIdx,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    matrixNumber=solverMappingVariable%matrixNumbers(matrixIdx,equationsIdx)
            
    EXITS("SolverMappingVariable_EquationMatrixNumberGet")
    RETURN
999 matrixNumber=0
    ERRORS("SolverMappingVariable_EquationMatrixNumberGet",err,error)
    EXITS("SolverMappingVariable_EquationMatrixNumberGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingVariable_EquationMatrixNumberGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of matrices for an equation in a solver mapping variable.
  SUBROUTINE SolverMappingVariable_EquationNumberOfMatricesGet(solverMappingVariable,equationsIdx,numberOfMatrices,err,error,*)

    !Argument variables
    TYPE(SolverMappingVariableType), POINTER :: solverMappingVariable !<A pointer to the solver mapping variable to get the number of matrices for
    INTEGER(INTG), INTENT(IN) :: equationsIdx !<The equations index to get the number of matrices in the equation that the solver mapping variable is mapped to.
    INTEGER(INTG), INTENT(OUT) :: numberOfMatrices !<On exit, the number of matrices that the equationsIdx'th equation that the solver mapping variable is mapped to.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMappingVariable_EquationNumberOfMatricesGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(solverMappingVariable)) CALL FlagError("Solver mapping variable is not associated.",err,error,*999)
    IF(equationsIdx<1.OR.equationsIdx>solverMappingVariable%numberOfEquations) THEN
      localError="The specified equations index of "//TRIM(NumberToVString(equationsIdx,"*",err,error))// &
        & " is invalid. The equations index should be >= 1 and <= "// &
        & TRIM(NumberToVString(solverMappingVariable%numberOfEquations,"*",err,error))// &
        & " for the solver mapping variable with variable index "// &
        & TRIM(NumberToVString(solverMappingVariable%variableIdx,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(solverMappingVariable%numberOfMatrices)) THEN
      localError="The number of matrices array is not allocated for the solver mapping variable with variable index "// &
        & TRIM(NumberToVString(solverMappingVariable%variableIdx,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    numberOfMatrices=solverMappingVariable%numberOfMatrices(equationsIdx)
            
    EXITS("SolverMappingVariable_EquationNumberOfMatricesGet")
    RETURN
999 numberOfMatrices=0
    ERRORS("SolverMappingVariable_EquationNumberOfMatricesGet",err,error)
    EXITS("SolverMappingVariable_EquationNumberOfMatricesGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingVariable_EquationNumberOfMatricesGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the equations field variable for a solver mapping variable.
  SUBROUTINE SolverMappingVariable_FieldVariableGet(solverMappingVariable,fieldVariable,err,error,*)

    !Argument variables
    TYPE(SolverMappingVariableType), POINTER :: solverMappingVariable !<A pointer to the solver mapping variable to get the field variables for
    TYPE(FieldVariableType), POINTER :: fieldVariable !<On exit, a pointer to the field variable for the solver mapping variable. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMappingVariable_FieldVariableGet",err,error,*998)

#ifdef WITH_PRECHECKS
    IF(ASSOCIATED(fieldVariable)) CALL FlagError("The field variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMappingVariable)) CALL FlagError("Solver mapping variable is not associated.",err,error,*999)
#endif    
    
    fieldVariable=>solverMappingVariable%variable

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(fieldVariable)) THEN
      localError="The field variable is not associated for the solver mapping variable with index "// &
        & TRIM(NumberToVString(solverMappingVariable%variableIdx,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    EXITS("SolverMappingVariable_FieldVariableGet")
    RETURN
999 NULLIFY(fieldVariable)
998 ERRORS("SolverMappingVariable_FieldVariableGet",err,error)
    EXITS("SolverMappingVariable_FieldVariableGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingVariable_FieldVariableGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of equations for a solver mapping variable.
  SUBROUTINE SolverMappingVariable_NumberOfEquationsGet(solverMappingVariable,numberOfEquations,err,error,*)

    !Argument variables
    TYPE(SolverMappingVariableType), POINTER :: solverMappingVariable !<A pointer to the solver mapping variable to get the number of equations for
    INTEGER(INTG), INTENT(OUT) :: numberOfEquations !<On exit, the number of equations that the solver mapping variable is mapped to.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("SolverMappingVariable_NumberOfEquationsGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(solverMappingVariable)) CALL FlagError("Solver mapping variable is not associated.",err,error,*999)
#endif    
    
    numberOfEquations=solverMappingVariable%numberOfEquations
            
    EXITS("SolverMappingVariable_NumberOfEquationsGet")
    RETURN
999 numberOfEquations=0
    ERRORS("SolverMappingVariable_NumberOfEquationsGet",err,error)
    EXITS("SolverMappingVariable_NumberOfEquationsGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingVariable_NumberOfEquationsGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the number of variables in the solver mapping variables list.
  SUBROUTINE SolverMappingVariables_NumberOfVariablesGet(solverMappingVariables,numberOfVariables,err,error,*)

    !Argument variables
    TYPE(SolverMappingVariablesType), POINTER :: solverMappingVariables !<A pointer to the solver variables list to get the number of variables for
    INTEGER(INTG), INTENT(OUT) :: numberOfVariables !<On exit, the number of varaibles int he solver variables list.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMappingVariables_NumberOfVariablesGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(solverMappingVariables)) CALL FlagError("Solver mapping variables is not associated.",err,error,*999)
#endif    
    
    numberOfVariables=solverMappingVariables%numberOfVariables
      
    EXITS("SolverMappingVariables_NumberOfVariablesGet")
    RETURN
999 ERRORS("SolverMappingVariables_NumberOfVariablesGet",err,error)
    EXITS("SolverMappingVariables_NumberOfVariablesGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingVariables_NumberOfVariablesGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns the specified solver mapping variable corresponding to the specified variable index in a list of solver mapping variables.
  SUBROUTINE SolverMappingVariables_VariableGet(solverMappingVariables,variableIdx,solverMappingVariable,err,error,*)

    !Argument variables
    TYPE(SolverMappingVariablesType), POINTER :: solverMappingVariables !<A pointer to the solver variables list to get the variable for
    INTEGER(INTG), INTENT(IN) :: variableIdx !<The variable index in the solver mapping variables list to get.
    TYPE(SolverMappingVariableType), POINTER :: solverMappingVariable !<On return, the specified solver mapping variable in the solver mapping variables list. Must not be associated on entry. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
    
    ENTERS("SolverMappingVariables_VariableGet",err,error,*998)

#ifdef WITH_PRECHECKS
    IF(ASSOCIATED(solverMappingVariable)) CALL FlagError("Solver mapping variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMappingVariables)) CALL FlagError("Solver mapping variables is not associated.",err,error,*999)
    IF(variableIdx<1.OR.variableIdx>solverMappingVariables%numberOfVariables) THEN
      localError="The specified variable index of "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
        & " is invalid. The variable index should be >= 1 and <= "// &
        & TRIM(NumberToVString(solverMappingVariables%numberOfVariables,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(solverMappingVariables%variables)) &
      & CALL FlagError("The solver mapping variables list variables is not allocated.",err,error,*999)
#endif    

    solverMappingVariable=>solverMappingVariables%variables(variableIdx)%ptr

#ifdef WITH_POSTCHECKS
    IF(.NOT.ASSOCIATED(solverMappingVariable)) THEN
      localError="The solver mapping field variable for variable index "// &
        & TRIM(NumberToVString(variableIdx,"*",err,error))//" is not associated."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    EXITS("SolverMappingVariables_VariableGet")
    RETURN
999 NULLIFY(solverMappingVariable)
998 ERRORS("SolverMappingVariables_VariableGet",err,error)
    EXITS("SolverMappingVariables_VariableGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingVariables_VariableGet
  
  !
  !================================================================================================================================
  !
  
  !>Checks to see if a field variable is in the solver mapping variables. If the field variable is in the list then solver mapping variable will point to the required solver mapping variable that corresponds to the field variable. If the field variable is not in the list then the returned solver mapping variable will be null. 
  SUBROUTINE SolverMappingVariables_VariableInListCheck(solverMappingVariables,fieldVariable,solverMappingVariable,err,error,*)

    !Argument variables
    TYPE(SolverMappingVariablesType), POINTER :: solverMappingVariables !<A pointer to the solver variables list to get the variables for
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to variable to check if it is the solver variables list.
    TYPE(SolverMappingVariableType), POINTER :: solverMappingVariable !<On return, the solver mapping variable will point to the solver mapping variable in the list of solver mapping variables that corresponds to the field variable. If the field variable is not in the list then the solver mapping variable will be null. The pointer should not be associated on entry. 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: variableIdx
    TYPE(FieldVariableType), POINTER :: solverVariable
    TYPE(SolverMappingVariableType), POINTER :: solverMappingTestVariable
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("SolverMappingVariables_VariableInListCheck",err,error,*998)

#ifdef WITH_PRECHECKS
    IF(ASSOCIATED(solverMappingVariable)) CALL FlagError("Solver mapping variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverMappingVariables)) CALL FlagError("Solver mapping variables is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("The field variable is not associated.",err,error,*999)
#endif    

    NULLIFY(solverMappingVariable)
    DO variableIdx=1,solverMappingVariables%numberOfVariables
      solverMappingTestVariable=>solverMappingVariables%variables(variableIdx)%ptr
      IF(ASSOCIATED(solverMappingTestVariable)) THEN
        solverVariable=>solverMappingTestVariable%variable
        IF(ASSOCIATED(fieldVariable,solverVariable)) THEN
          solverMappingVariable=>solverMappingTestVariable
          EXIT
        ENDIF
      ELSE
        localError="The solver mapping variable is not associated for variable index "// &
          & TRIM(NumberToVString(variableIdx,"*",err,error))//" is not associated."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDDO !variableIdx
    
    EXITS("SolverMappingVariables_VariableInListCheck")
    RETURN
999 NULLIFY(solverMappingVariable)
998 ERRORS("SolverMappingVariables_VariableInListCheck",err,error)
    EXITS("SolverMappingVariables_VariableInListCheck")
    RETURN 1
    
  END SUBROUTINE SolverMappingVariables_VariableInListCheck
  
  !
  !================================================================================================================================
  !
  
  !>Returns the solver DOF coupling information from a variable DOF in a variable DOF to SolverDOFs map
  SUBROUTINE SolverMappingVDOFToSDOFsMap_SolverDOFCouplingGet(variableDOFToSolverDOFsMap,variableDOFIdx,solverDOFIdx, &
    & couplingCoefficient,additiveConstant,err,error,*)    

    !Argument variables
    TYPE(VariableDOFToSolverDOFsMapType), POINTER :: variableDOFToSolverDOFsMap !<A pointer to the variable DOF to solver DOFs map to get the 
    INTEGER(INTG), INTENT(IN) :: variableDOFIdx !<The variable DOF index in the variable DOF to sovler DOFs map to get the solver DOF coupling information.
    INTEGER(INTG), INTENT(OUT) :: solverDOFIdx !<On exit, the solver DOF that the variable DOF is coupled to
    REAL(DP), INTENT(OUT) :: couplingCoefficient !<On exit, the coupling coefficient between the variable DOF and the solver DOF
    REAL(DP), INTENT(OUT) :: additiveConstant !<On exit, the additive constant (offset) between the variable DOF and the solver DOF
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("SolverMappingVDOFToSDOFsMap_SolverDOFCouplingGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(variableDOFToSolverDOFsMap)) &
      & CALL FlagError("The variable DOF to solver DOFs map is not associated.",err,error,*999)    
    IF(.NOT.ALLOCATED(variableDOFToSolverDOFsMap%dofNumbers)) &
      & CALL FlagError("The dof numbers array is not allocated for the variable DOF to solver DOFs map.",err,error,*999)
    IF(variableDOFIdx<1.OR.variableDOFIdx>SIZE(variableDOFToSolverDOFsMap%dofNumbers,1)) THEN
      localError="The specified variable DOF index of "//TRIM(NumberToVString(variableDOFIdx,"*",err,error))// &
        & " is invalid. The variable index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(variableDOFToSolverDOFsMap%dofNumbers,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(variableDOFToSolverDOFsMap%couplingCoefficients)) &
      & CALL FlagError("The coupling coefficients array is not allocated for the variable DOF to solver DOFs map.",err,error,*999)
    IF(.NOT.ALLOCATED(variableDOFToSolverDOFsMap%additiveConstants)) &
      & CALL FlagError("The additive constants array is not allocated for the variable DOF to solver DOFs map.",err,error,*999)
#endif    
    
    solverDOFIdx=variableDOFToSolverDOFsMap%dofNumbers(variableDOFIdx)
    couplingCoefficient=variableDOFToSolverDOFsMap%couplingCoefficients(variableDOFIdx)
    additiveConstant=variableDOFToSolverDOFsMap%additiveConstants(variableDOFIdx)
            
    EXITS("SolverMappingVDOFToSDOFsMap_SolverDOFCouplingGet")
    RETURN
999 ERRORS("SolverMappingVDOFToSDOFsMap_SolverDOFCouplingGet",err,error)
    EXITS("SolverMappingVDOFToSDOFsMap_SolverDOFCouplingGet")
    RETURN 1
    
  END SUBROUTINE SolverMappingVDOFToSDOFsMap_SolverDOFCouplingGet
  
  !
  !================================================================================================================================
  !
  
END MODULE SolverMappingAccessRoutines
