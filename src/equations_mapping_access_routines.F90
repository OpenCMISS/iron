!> \file
!> \author Chris Bradley
!> \brief This module contains all equations mapping access method routines.
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

!> This module contains all equations mapping access method routines.
MODULE EquationsMappingAccessRoutines
  
  USE BaseRoutines
  USE ISO_VARYING_STRING
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

  PUBLIC EquationsMappingDynamic_DynamicVariableGet

  PUBLIC EquationsMappingDynamic_DynamicVariableTypeGet

  PUBLIC EquationsMappingDynamic_EquationsMatrixToVariableMapGet

  PUBLIC EquationsMappingDynamic_NumberOfDynamicMatricesGet

  PUBLIC EquationsMappingDynamic_VariableNumberOfMatricesGet

  PUBLIC EquationsMappingDynamic_VariableToEquationsMatricesMapGet

  PUBLIC EquationsMappingLHS_EquationsRowToLHSDOFMapGet

  PUBLIC EquationsMappingLHS_LHSDOFToEquationsRowMapGet

  PUBLIC EquationsMappingLHS_LHSVariableGet

  PUBLIC EquationsMappingLHS_LHSVariableTypeGet
  
  PUBLIC EquationsMappingLHS_RowDOFsMappingGet

  PUBLIC EquationsMappingLHS_VectorMappingGet

  PUBLIC EquationsMappingLinear_EquationsMatrixToVariableMapGet

  PUBLIC EquationsMappingLinear_LinearMatrixVariableGet

  PUBLIC EquationsMappingLinear_LinearMatrixVariableTypeGet

  PUBLIC EquationsMappingLinear_LinearVariableGet

  PUBLIC EquationsMappingLinear_LinearVariableIndexGet

  PUBLIC EquaitonsMappingLinear_NumberOfLinearMatricesGet

  PUBLIC EquaitonsMappingLinear_NumberOfLinearVariablesGet

  PUBLIC EquationsMappingLinear_VariableNumberOfMatricesGet

  PUBLIC EquationsMappingLinear_VariableToEquationsMatricesMapGet

  PUBLIC EquationsMappingNonlinear_NumberOfResidualsGet

  PUBLIC EquationsMappingNonlinear_ResidualMappingGet

  PUBLIC EquationsMappingResidual_JacobianMatrixToVariableGet

  PUBLIC EquationsMappingResidual_JacobianMatrixVariableGet

  PUBLIC EquationsMappingResidual_JacobianMatrixVariableTypeGet

  PUBLIC EquationsMappingResidual_NumberOfJacobianMatricesGet
  
  PUBLIC EquationsMappingResidual_NumberOfResidualVariablesGet
  
  PUBLIC EquationsMappingResidual_VariableGet

  PUBLIC EquationsMappingResidual_VariableIndexGet
  
  PUBLIC EquaitonsMappingResidual_VariableToJacobianMatrixMapGet
  
  PUBLIC EquationsMappingResidual_VariableTypeGet

  PUBLIC EquationsMappingRHS_RHSVariableGet
  
  PUBLIC EquationsMappingRHS_RHSVariableTypeGet
  
  PUBLIC EquationsMappingScalar_AssertIsFinish,EquationsMappingScalar_AssertNotFinished
  
  PUBLIC EquationsMappingScalar_CreateValuesCacheGet
  
  PUBLIC EquationsMappingScalar_ScalarEquationsGet

  PUBLIC EquationsMappingSource_SourceVariableGet

  PUBLIC EquationsMappingSources_SourceMappingGet

  PUBLIC EquationsMappingVector_AssertIsFinish,EquationsMappingVector_AssertNotFinished
  
  PUBLIC EquationsMappingVector_CreateValuesCacheGet

  PUBLIC EquationsMappingVector_DynamicMappingExists
  
  PUBLIC EquationsMappingVector_DynamicMappingGet
  
  PUBLIC EquationsMappingVector_LHSMappingGet
  
  PUBLIC EquationsMappingVector_LinearMappingExists
  
  PUBLIC EquationsMappingVector_LinearMappingGet
  
  PUBLIC EquationsMappingVector_NonlinearMappingExists
  
  PUBLIC EquationsMappingVector_NonlinearMappingGet
  
  PUBLIC EquationsMappingVector_RHSMappingExists
  
  PUBLIC EquationsMappingVector_RHSMappingGet
  
  PUBLIC EquationsMappingVector_SourcesMappingExists

  PUBLIC EquationsMappingVector_SourcesMappingGet

  PUBLIC EquationsMappingVector_VectorEquationsGet
  
  PUBLIC EquationsMappingVector_VectorMatricesGet
 
  PUBLIC EquationsMappingVectorCVC_DynamicMatrixCoefficientGet
  
  PUBLIC EquationsMappingVectorCVC_DynamicVariableTypeGet
  
  PUBLIC EquationsMappingVectorCVC_LinearMatrixCoefficientGet
  
  PUBLIC EquationsMappingVectorCVC_LinearMatrixVariableTypeGet

  PUBLIC EquationsMappingVectorCVC_NumberOfResidualVariablesGet
  
  PUBLIC EquationsMappingVectorCVC_ResidualCoefficientGet

  PUBLIC EquationsMappingVectorCVC_ResidualVariableTypeGet

  PUBLIC EquationsMappingVectorCVC_SourceVectorCoefficientGet
  
  PUBLIC EquationsMappingVectorCVC_SourceVectorVariableTypeGet

  PUBLIC EquationsMappingVectorEMToVMap_EquationsMatrixGet

  PUBLIC EquationsMappingVectorEMToVMap_VariableGet

  PUBLIC EquationsMappingVectorJMToVMap_VariableGet

  PUBLIC EquationsMappingVectorJMToVMap_JacobianMatrixGet

  PUBLIC EquationsMappingVectorVToEMSMap_VariableGet

  PUBLIC EquationsMappingVectorVToJMMap_VariableGet

CONTAINS

  !
  !================================================================================================================================
  !

  !>Gets the specified dynamic variable for a dynamic mapping.
  SUBROUTINE EquationsMappingDynamic_DynamicVariableGet(dynamicMapping,dynamicVariable,err,error,*)

    !Argument variables
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping !<A pointer to the dynamic mapping to get the dynamic variable for
    TYPE(FieldVariableType), POINTER :: dynamicVariable !<On exit, a pointer to the requested dynamic variable. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingDynamic_DynamicVariableGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dynamicVariable)) CALL FlagError("Dynamic variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(dynamicMapping)) CALL FlagError("Dynamic mapping is not associated.",err,error,*999)
#endif    
    
    dynamicVariable=>dynamicMapping%dynamicVariable

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(dynamicVariable)) CALL FlagError("The dynamic variable is not associated.",err,error,*999)
#endif    
    
    EXITS("EquationsMappingDynamic_DynamicVariableGet")
    RETURN
999 NULLIFY(dynamicVariable)
998 ERRORS("EquationsMappingDynamic_DynamicVariableGet",err,error)
    EXITS("EquationsMappingDynamic_DynamicVariableGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingDynamic_DynamicVariableGet

  !
  !================================================================================================================================
  !

  !>Gets the specified dynamic variable for a dynamic mapping.
  SUBROUTINE EquationsMappingDynamic_DynamicVariableGet(dynamicMapping,dynamicVariable,err,error,*)

    !Argument variables
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping !<A pointer to the dynamic mapping to get the dynamic variable for
    TYPE(FieldVariableType), POINTER :: dynamicVariable !<On exit, a pointer to the requested dynamic variable. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingDynamic_DynamicVariableGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dynamicVariable)) CALL FlagError("Dynamic variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(dynamicMapping)) CALL FlagError("Dynamic mapping is not associated.",err,error,*999)
#endif    
    
    dynamicVariable=>dynamicMapping%dynamicVariable

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(dynamicVariable)) CALL FlagError("The dynamic variable is not associated.",err,error,*999)
#endif    
    
    EXITS("EquationsMappingDynamic_DynamicVariableGet")
    RETURN
999 NULLIFY(dynamicVariable)
998 ERRORS("EquationsMappingDynamic_DynamicVariableGet",err,error)
    EXITS("EquationsMappingDynamic_DynamicVariableGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingDynamic_DynamicVariableGet

  !
  !================================================================================================================================
  !

  !>Gets the specified dynamic variable type for a dynamic mapping.
  SUBROUTINE EquationsMappingDynamic_DynamicVariableTypeGet(dynamicMapping,dynamicVariableType,err,error,*)

    !Argument variables
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping !<A pointer to the dynamic mapping to get the dynamic variable typefor
    INTEGER(INTG), INTENT(OUT) :: dynamicVariableType !<On exit, the requested dynamic variable type.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingDynamic_DynamicVariableTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(dynamicMapping)) CALL FlagError("Dynamic mapping is not associated.",err,error,*999)
#endif    
    
    dynamicVariableType=dynamicMapping%dynamicVariableType 
    
    EXITS("EquationsMappingDynamic_DynamicVariableTypeGet")
    RETURN
999 ERRORS("EquationsMappingDynamic_DynamicVariableTypeGet",err,error)
    EXITS("EquationsMappingDynamic_DynamicVariableTypeGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingDynamic_DynamicVariableTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the equations matrix to variable mapping a dynamic mapping.
  SUBROUTINE EquationsMappingDynamic_EquationsMatrixToVariableMapGet(dynamicMapping,dynamicMatrixIdx,equationsMatrixToVariableMap, &
    & err,error,*)

    !Argument variables
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping !<A pointer to the dynamic mapping to get the equations matrix to variable map for
    INTEGER(INTG), INTENT(IN) :: dynamicMatrixIdx !<The dynamic matrix index to get the equations matrix to variable map for
    TYPE(EquationsMatrixToVarMapType), POINTER :: equationsMatrixToVariableMap !<On exit, a pointer to the requested equations matrix to variable map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("EquationsMappingDynamic_EquationsMatrixToVariableMap",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(equationsMatrixToVariableMap)) &
      & CALL FlagError("Equations matrix to variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(dynamicMapping)) CALL FlagError("Dynamic mapping is not associated.",err,error,*999)
    IF(dynamicMatrixIdx<1.OR.dynamicMatrixIdx>dynamicMapping%numberOfDynamicMatrices) THEN
      localError="The specified dynamic matrix index of "//TRIM(NumberToVString(dynamicMatrixIdx,"*",err,error))// &
        & " is invalid. The dynamic matrix index should be >=1 and <= "// &
        & TRIM(NumberToVString(dynamicMapping%numberOfDynamicMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(dynamicMapping%equationsMatrixToVarMaps)) &
      & CALL FlagError("The equations matrix to variable maps is not allocated for the dynamic mapping.",err,error,*999)
#endif    
    
    equationsMatrixToVariableMap=>dynamicMapping%equationsMatrixToVarMaps(dynamicMatrixIdx)%ptr

#ifdef    
    IF(.NOT.ASSOCIATED(equationsMatrixToVariableMap)) THEN
      localError="The equations matrix to variable map is not associated for dynamic matrix index "// &
        & TRIM(NumberToVString(dynamicMatrixIdx,"*",err,error))//" of the dynamic mapping."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    EXITS("EquationsMappingDynamic_EquationsMatrixToVariableMapGet")
    RETURN
999 NULLIFY(equationsMatrixToVariableMap)
998 ERRORS("EquationsMappingDynamic_EquationsMatrixToVariableMapGet",err,error)
    EXITS("EquationsMappingDynamic_EquationsMatrixToVariableMapGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingDynamic_EquationsMatrixToVariableMapGet

  !
  !================================================================================================================================
  !

  !>Gets the number of dynamic matrices for a dynamic mapping.
  SUBROUTINE EquationsMappingDynamic_NumberOfDynamicMatricesGet(dynamicMapping,numberOfDynamicMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping !<A pointer to the dynamic mapping to get the number of dynamic matrices for
    INTEGER(INTG), INTENT(OUT) :: numberOfDynamicMatrices !<On exit, the number of dynamic matrices.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingDynamic_NumberOfDynamicMatricesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(dynamicMapping)) CALL FlagError("Dynamic mapping is not associated.",err,error,*999)
#endif    
    
    numberOfDynamicMatrices=dynamicMapping%numberOfDynamicMatrices 
    
    EXITS("EquationsMappingDynamic_NumberOfDynamicMatricesGet")
    RETURN
999 ERRORS("EquationsMappingDynamic_NumberOfDynamicMatricesGet",err,error)
    EXITS("EquationsMappingDynamic_NumberOfDynamicMatricesGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingDynamic_NumberOfDynamicMatricesGet

  !
  !================================================================================================================================
  !

  !>Gets the number of equations matrices a dynamic variable is mapped to in a dynamic mapping.
  SUBROUTINE EquationsMappingDynamic_VariableNumberOfMatricesGet(dynamicMapping,numberOfEquationsMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping !<A pointer to the dynamic mapping to get the number of equations matrices for
    INTEGER(INTG), INTENT(OUT) :: numberOfEquationsMatrices !<On exit, the number of equations matrices that this dynamic variable is mapped to in this linear mapping.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(varToEquationsMatricesMapType), POINTER :: variableToEquationsMatricesMap
#endif    
 
    ENTERS("EquationsMappingDynamic_VariableNumberOfMatricesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    NULLIFY(variableToEquationsMatricesMap)
    CALL EquationsMappingDynamic_VariableToEquationsMatricesMapGet(dynamicMapping,variableToEquationsMatricesMap,err,error,*999)
#endif    
    
    numberOfEquationsMatrices=dynamicMapping%varToEquationsMatricesMap%numberOfEquationsMatrices
   
    EXITS("EquationsMappingDynamic_VariableNumberOfMatricesGet")
    RETURN
999 ERRORS("EquationsMappingDynamic_VariableNumberOfMatricesGet",err,error)
    EXITS("EquationsMappingDynamic_VariableNumberOfMatricesGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingDynamic_VariableNumberOfMatricesGet

  !
  !================================================================================================================================
  !

  !>Gets the variable to equations matrices map for the dynamic variable in a dynamic mapping.
  SUBROUTINE EquationsMappingDynamic_VariableToEquationsMatricesMapGet(dynamicMapping,variableToEquationsMatricesMap,err,error,*)

    !Argument variables
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping !<A pointer to the dynamic mapping to get the variable to equations matrices map for
    TYPE(varToEquationsMatricesMapType), POINTER :: variableToEquationsMatricesMap !<On exit, a pointer to the requested dynamic variable to equations matrices map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingDynamic_VariableToEquationsMatricesMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(variableToEquationsMatricesMap)) &
      & CALL FlagError("Variable to equations matrices map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(dynamicMapping)) CALL FlagError("Dynamic mapping is not associated.",err,error,*999)
#endif    
    
    variableToEquationsMatricesMap=>dynamicMapping%varToEquationsMatricesMap

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(variableToEquationsMatricesMap)) &
      & CALL FlagError("The dynamic variable to equations matrices map is not associated.",err,error,*999)
#endif    
    
    EXITS("EquationsMappingDynamic_VariableToEquationsMatricesMapGet")
    RETURN
999 NULLIFY(variableToEquationsMatricesMap)
998 ERRORS("EquationsMappingDynamic_VariableToEquationsMatricesMapGet",err,error)
    EXITS("EquationsMappingDynamic_VariableToEquationsMatricesMapGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingDynamic_VariableToEquationsMatricesMapGet

  !
  !================================================================================================================================
  !

  !>Gets the equations row to LHS DOF map for a LHS mapping.
  SUBROUTINE EquationsMappingLHS_EquationsRowToLHSDOFMapGet(lhsMapping,equationsRowToLHSDOFMap,err,error,*)

    !Argument variables
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping !<A pointer to the LHS mapping to get the equations row to LHS DOF map for
    INTEGER(INTG), POINTER :: equationsRowToLHSDOMap(:) !<On exit, a pointer to the equations row to LHS DOF map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingLHS_EquationsRowToLHSDOFGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(equationsRowToLHSDOFMap)) CALL FlagError("Equations row to LHS DOF map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(lhsMapping)) CALL FlagError("LHS mapping is not associated.",err,error,*999)
#endif    
   
    equationsRowToLHSDOFMap=>lhsMapping%equationsRowToLHSDOFMap

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(equationsRowToLHSDOFMap)) &
      & CALL FlagError("The equations row to LHS DOF map is not associated for the LHS mapping.",err,error,*999)
#endif    
    
    EXITS("EquationsMappingLHS_EquationRowToLHSDOFMapGet")
    RETURN
999 NULLIFY(equationsRowToLHSDOFMap)
998 ERRORS("EquationsMappingLHS_EquationsRowToLHSDOFMapGet",err,error)
    EXITS("EquationsMappingLHS_EquationsRowToLHSDOFMapGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingLHS_EquationsRowToLHSDOFMapGet

  !
  !================================================================================================================================
  !

  !>Gets the LHS DOF to equations row map for a LHS mapping.
  SUBROUTINE EquationsMappingLHS_LHSDOFToEquationsRowMapGet(lhsMapping,lhsDOFToEquationsRowMap,err,error,*)

    !Argument variables
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping !<A pointer to the LHS mapping to get the LHS DOF to equations row map for
    INTEGER(INTG), POINTER :: lhsDOFToEquationsRowMap(:) !<On exit, a pointer to the LHS DOF to equations row map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingLHS_LHSDOFToEquationsRowGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(lhsDOFToEquationsRowMap)) CALL FlagError("LHS DOF to equations row map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(lhsMapping)) CALL FlagError("LHS mapping is not associated.",err,error,*999)
#endif    
   
    lhsDOFToEquationsRowMap=>lhsMapping%lhsDOFToEquationsRowMap

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(lhsDOFToEquationsRowMap)) &
      & CALL FlagError("The LHS DOF to equations row map is not associated for the LHS mapping.",err,error,*999)
#endif    
    
    EXITS("EquationsMappingLHS_LHSDOFToEquationRowMapGet")
    RETURN
999 NULLIFY(lhsDOFToEquationsRowMap)
998 ERRORS("EquationsMappingLHS_LHSDOFToEquationsRowMapGet",err,error)
    EXITS("EquationsMappingLHS_LHSDOFToEquationsRowMapGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingLHS_LHSDOFToEquationsRowMapGet

  !
  !================================================================================================================================
  !

  !>Gets the specified LHS variable type for a LHS mapping.
  SUBROUTINE EquationsMappingLHS_LHSVariableTypeGet(lhsMapping,lhsVariableType,err,error,*)

    !Argument variables
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping !<A pointer to the LHS mapping to get the LHS variable type for
    INTEGER(INTG), INTENT(OUT) :: lhsVariableType !<On exit, the LHS variable type.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingLHS_LHSVariableTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(lhsMapping)) CALL FlagError("LHS mapping is not associated.",err,error,*999)
#endif    
   
    lhsVariableType=lhsMapping%lhsVariableType
    
    EXITS("EquationsMappingLHS_LHSVariableTypeGet")
    RETURN
999 ERRORS("EquationsMappingLHS_LHSVariableTypeGet",err,error)
    EXITS("EquationsMappingLHS_LHSVariableTypeGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingLHS_LHSVariableTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the row DOFs mapping for a LHS mapping.
  SUBROUTINE EquationsMappingLHS_RowDOFsMappingGet(lhsMapping,rowDOFsMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping !<A pointer to the LHS mapping to get the row DOFs mapping for
    TYPE(DomainMappingType), POINTER :: rowDOFsMapping !<On exit, a pointer to the row DOFs mapping for the LHS mapping. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingLHS_RowDOFsMappingGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(rowDOFsMapping)) CALL FlagError("Row DOFs mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(lhsMapping)) CALL FlagError("LHS mapping is not associated.",err,error,*999)
#endif    
   
    rowDOFsMapping=>lhsMapping%rowDOFsMapping

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(rowDOFsMapping)) CALL FlagError("The row DOFs mapping is not associated for the LHS mapping.",err,error,*999)
#endif    
    
    EXITS("EquationsMappingLHS_RowDOFsMappingGet")
    RETURN
999 NULLIFY(rowDOFsMapping)
998 ERRORS("EquationsMappingLHS_RowDOFsMappingGet",err,error)
    EXITS("EquationsMappingLHS_RowDOFsMappingGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingLHS_RowDOFsMappingGet

  !
  !================================================================================================================================
  !

  !>Gets the vector mapping for a LHS mapping.
  SUBROUTINE EquationsMappingLHS_VectorMappingGet(lhsMapping,vectorMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping !<A pointer to the LHS mapping to get the vector mapping for
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<On exit, a pointer to the vector mapping for the LHS mapping. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingLHS_VectorMappingGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(vectorMapping)) CALL FlagError("Vector mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(lhsMapping)) CALL FlagError("LHS mapping is not associated.",err,error,*999)
#endif    
   
    vectorMapping=>lhsMapping%vectorMapping

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("The vector mapping is not associated for the LHS mapping.",err,error,*999)
#endif    
    
    EXITS("EquationsMappingLHS_VectorMappingGet")
    RETURN
999 NULLIFY(vectorMapping)
998 ERRORS("EquationsMappingLHS_VectorMappingGet",err,error)
    EXITS("EquationsMappingLHS_VectorMappingGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingLHS_VectorMappingGet

  !
  !================================================================================================================================
  !

  !>Gets the equations matrix to variable mapping a linear mapping.
  SUBROUTINE EquationsMappingLinear_EquationsMatrixToVariableMapGet(linearMapping,linearMatrixIdx,equationsMatrixToVariableMap, &
    & err,error,*)

    !Argument variables
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping !<A pointer to the linear mapping to get the equations matrix to variable map for
    INTEGER(INTG), INTENT(IN) :: linearMatrixIdx !<The linear matrix index to get the equations matrix to variable map for
    TYPE(EquationsMatrixToVarMapType), POINTER :: equationsMatrixToVariableMap !<On exit, a pointer to the requested equations matrix to variable map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("EquationsMappingLinear_EquationsMatrixToVariableMap",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(equationsMatrixToVariableMap)) &
      & CALL FlagError("Equations matrix to variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(linearMapping)) CALL FlagError("Linear mapping is not associated.",err,error,*999)
    IF(linearMatrixIdx<1.OR.linearMatrixIdx>linearMapping%numberOfLinearMatrices) THEN
      localError="The specified linear matrix index of "//TRIM(NumberToVString(linearMatrixIdx,"*",err,error))// &
        & " is invalid. The linear matrix index should be >=1 and <= "// &
        & TRIM(NumberToVString(linearMapping%numberOfLinearMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(linearMapping%equationsMatrixToVarMaps)) &
      & CALL FlagError("The equations matrix to variable maps is not allocated for the linear mapping.",err,error,*999)
#endif    
    
    equationsMatrixToVariableMap=>linearMapping%equationsMatrixToVarMaps(linearMatrixIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(equationsMatrixToVariableMap)) THEN
      localError="The equations matrix to variable map is not associated for linear matrix index "// &
        & TRIM(NumberToVString(linearMatrixIdx,"*",err,error))//" of the linear mapping."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    EXITS("EquationsMappingLinear_EquationsMatrixToVariableMapGet")
    RETURN
999 NULLIFY(equationsMatrixToVariableMap)
998 ERRORS("EquationsMappingLinear_EquationsMatrixToVariableMapGet",err,error)
    EXITS("EquationsMappingLinear_EquationsMatrixToVariableMapGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingLinear_EquationsMatrixToVariableMapGet
  
  !
  !================================================================================================================================
  !

  !>Gets the specified linear matrix variable for a linear matrix mapping.
  SUBROUTINE EquationsMappingLinear_LinearMatrixVariableGet(linearMapping,linearMatrixIdx,linearMatrixVariable,err,error,*)

    !Argument variables
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping !<A pointer to the linear mapping to get the linear matrix variable for
    INTEGER(INTG), INTENT(IN) :: linearMatrixIdx !<The index of the linear matrix to get the linear variable for. 
    TYPE(FieldVariableType), POINTER :: linearMatrixVariable !<On exit, a pointer to the requested linear matrix variable. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(EquationsMatrixToVarMapType), POINTER :: equationsMatrixToVariableMap
#endif
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("EquationsMappingLinear_LinearMatrixVariableGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(linearMatrixVariable)) CALL FlagError("Linear matrix variable is already associated.",err,error,*998)
    NULLIFY(equationsMatrixToVariableMap)
    CALL EquationsMappingLinear_EquationsMatrixToVarMapGet(linearMapping,linearMatrixIdx,equationsMatrixToVarMap,err,error,*999)
#endif    
   
    linearMatrixVariable=>linearMapping%equationsMatrixToVarMaps(linearMatrixIdx)%ptr%variable

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(linearMatrixVariable)) THEN
      localError="The linear matrix variable is not associated for linear matrix index "// &
        & TRIM(NumberToVString(linearMatrixIdx,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    EXITS("EquationsMappingLinear_LinearMatrixVariableGet")
    RETURN
999 NULLIFY(linearMatrixVariable)
998 ERRORS("EquationsMappingLinear_LinearMatrixVariableGet",err,error)
    EXITS("EquationsMappingLinear_LinearMatrixVariableGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingLinear_LinearMatrixVariableGet

  !
  !================================================================================================================================
  !

  !>Gets the specified linear matrix variable type for a linear matrix mapping.
  SUBROUTINE EquationsMappingLinear_LinearMatrixVariableTypeGet(linearMapping,linearMatrixIdx,linearMatrixVariableType,err,error,*)

    !Argument variables
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping !<A pointer to the linear mapping to get the linear matrix variable type for
    INTEGER(INTG), INTENT(IN) :: linearMatrixIdx !<The index of the linear matrix to get the linear variable for. 
    INTEGER(INTG), INTENT(OUT) :: linearMatrixVariableType !<On exit, the requested linear matrix variable type.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(EquationsMatrixToVarMapType), POINTER :: equationsMatrixToVariableMap
#endif
    
    ENTERS("EquationsMappingLinear_LinearMatrixVariableTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(linearMapping)) CALL FlagError("Linear mapping is not associated.",err,error,*999)
    NULLIFY(equationsMatrixToVariableMap)
    NULLIFY(equationsMatrixToVariableMap)
    CALL EquationsMappingLinear_EquationsMatrixToVarMapGet(linearMapping,linearMatrixIdx,equationsMatrixToVarMap,err,error,*999)
#endif    
    
    linearMatrixVariableType=linearMapping%equationsMatrixToVarMaps(linearMatrixIdx)%ptr%variableType
    
    EXITS("EquationsMappingLinear_LinearMatrixVariableTypeGet")
    RETURN
999 ERRORS("EquationsMappingLinear_LinearMatrixVariableTypeGet",err,error)
    EXITS("EquationsMappingLinear_LinearMatrixVariableTypeGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingLinear_LinearMatrixVariableTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the specified linear variable for a linear mapping.
  SUBROUTINE EquationsMappingLinear_LinearVariableGet(linearMapping,linearVariableIdx,linearVariable,err,error,*)

    !Argument variables
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping !<A pointer to the linear mapping to get the linear variable for
    INTEGER(INTG), INTENT(IN) :: linearVariableIdx !<The index of the linear variable to get the linear variable for. 
    TYPE(FieldVariableType), POINTER :: linearVariable !<On exit, a pointer to the requested linear variable. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("EquationsMappingLinear_LinearVariableGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(linearVariable)) CALL FlagError("Linear variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(linearMapping)) CALL FlagError("Linear mapping is not associated.",err,error,*999)
    IF(linearVariableIdx<1.OR.linearVariableIdx>linearMapping%numberOfLinearVariables) THEN
      localError="The specified linear variable index of "//TRIM(NumberToVString(linearVariableIdx,"*",err,error))// &
        & " is invalid for a linear mapping which has "// &
        & TRIM(NumberToVString(linearMapping%numberOfLinearVariables,"*",err,error))//" linear variables."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(linearMapping%linearVariables)) &
      & CALL FlagError("Linear mapping linear variables is not allocated.",err,error,*999)
#endif    
    
    linearVariable=>linearMapping%linearVariables(linearVariableIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(linearVariable)) THEN
      localError="The linear variable is not associated for linear variable index "// &
        & TRIM(NumberToVString(linearVariableIdx,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    EXITS("EquationsMappingLinear_LinearVariableGet")
    RETURN
999 NULLIFY(linearVariable)
998 ERRORS("EquationsMappingLinear_LinearVariableGet",err,error)
    EXITS("EquationsMappingLinear_LinearVariableGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingLinear_LinearVariableGet

  !
  !================================================================================================================================
  !

  !>Gets the specified linear variable index for a variable type in a linear mapping.
  SUBROUTINE EquationsMappingLinear_LinearVariableIndexGet(linearMapping,linearVariableType,linearVariableIdx,err,error,*)

    !Argument variables
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping !<A pointer to the linear mapping to get the linear variable index for
    INTEGER(INTG), INTENT(IN) :: linearVariableType !<The type of the linear variable to get the linear variable index for. 
    INTEGER(INTG), INTENT(OUT) :: linearVariableIdx !<On exit, the requested linear variable index corresponding to the variable type. If there is no variable index the returned value will be zero.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("EquationsMappingLinear_LinearVariableIndexGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(linearMapping)) CALL FlagError("Linear mapping is not associated.",err,error,*999)    
    IF(linearVariableType<1.OR.linearVariableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified linear variable type of "//TRIM(NumberToVString(linearVariabletType,"*",err,error))// &
        & " is invalid. The linear variable type should be >=1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(linearMapping%linearVariableTypesMap)) &
      & CALL FlagError("Linear mapping linear variables type map is not allocated.",err,error,*999)
#endif    
    
    linearVariableIdx=linearMapping%linearVariableTypeMap(linearVariableType)
   
    EXITS("EquationsMappingLinear_LinearVariableIndexGet")
    RETURN
999 ERRORS("EquationsMappingLinear_LinearVariableIndexGet",err,error)
    EXITS("EquationsMappingLinear_LinearVariableIndexGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingLinear_LinearVariableIndexGet

  !
  !================================================================================================================================
  !

  !>Gets the number of linear matrices for a linear mapping.
  SUBROUTINE EquationsMappingLinear_NumberOfLinearMatricesGet(linearMapping,numberOfLinearMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMappingLinearType), POINTER :: LinearMapping !<A pointer to the linear mapping to get the number of linear matrices for
    INTEGER(INTG), INTENT(OUT) :: numberOfLinearMatrices !<On exit, the number of linear matrices.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingLinear_NumberOfLinearMatricesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(linearMapping)) CALL FlagError("Linear mapping is not associated.",err,error,*999)
#endif    
    
    numberOfLinearMatrices=linearMapping%numberOfLinearMatrices 
    
    EXITS("EquationsMappingLinear_NumberOfLinearMatricesGet")
    RETURN
999 ERRORS("EquationsMappingLinear_NumberOfLinearMatricesGet",err,error)
    EXITS("EquationsMappingLinear_NumberOfLinearMatricesGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingLinear_NumberOfLinearMatricesGet

  !
  !================================================================================================================================
  !

  !>Gets the number of linear variables for a linear mapping.
  SUBROUTINE EquationsMappingLinear_NumberOfLinearVariablesGet(linearMapping,numberOfLinearVariables,err,error,*)

    !Argument variables
    TYPE(EquationsMappingLinearType), POINTER :: LinearMapping !<A pointer to the linear mapping to get the number of linear variables for
    INTEGER(INTG), INTENT(OUT) :: numberOfLinearVariables !<On exit, the number of linear variables.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingLinear_NumberOfLinearVariablesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(linearMapping)) CALL FlagError("Linear mapping is not associated.",err,error,*999)
#endif    
    
    numberOfLinearVariables=linearMapping%numberOfLinearVariables 
    
    EXITS("EquationsMappingLinear_NumberOfLinearVariablesGet")
    RETURN
999 ERRORS("EquationsMappingLinear_NumberOfLinearVariblesGet",err,error)
    EXITS("EquationsMappingLinear_NumberOfLinearVariablesGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingLinear_NumberOfLinearVariablesGet

  !
  !================================================================================================================================
  !

  !>Gets the number of equations matrices a linear variable index is mapped to in a linear mapping.
  SUBROUTINE EquationsMappingLinear_VariableNumberOfMatricesGet(linearMapping,linearVariableIdx,numberOfEquationsMatrices, &
    & err,error,*)

    !Argument variables
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping !<A pointer to the linear mapping to get the number of matrices for
    INTEGER(INTG), INTENT(IN) :: linearVariableIdx !<The linear variable index to get the number of equations matrices for
    INTEGER(INTG), INTENT(OUT) :: numberOfEquationsMatrices !<On exit, the number of equations matrices that this linear variable is mapped to in this linear mapping.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(varToEquationsMatricesMapType), POINTER :: variableToEquationsMatricesMap
#endif    
 
    ENTERS("EquationsMappingLinear_VariableNumberOfMatricesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    NULLIFY(variableToEquationsMatricesMap)    
    CALL EquationsMappingLinear_VariableToEquationsMatricesMapGet(linearMapping,linearVariableIdx,variableToEquationsMatricesMap, &
      & err,error,*999)
#endif    
    
    numberOfEquationsMatrices=linearMapping%varToEquationsMatricesMaps(linearVariableIdx)%ptr%numberOfEquationsMatrices
   
    EXITS("EquationsMappingLinear_VariableNumberOfMatricesGet")
    RETURN
999 ERRORS("EquationsMappingLinear_VariableNumberOfMatricesGet",err,error)
    EXITS("EquationsMappingLinear_VariableNumberOfMatricesGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingLinear_VariableNumberOfMatricesGet

  !
  !================================================================================================================================
  !

  !>Gets the variable to equations matrices map for the linear variable in a linear mapping.
  SUBROUTINE EquationsMappingLinear_VariableToEquationsMatricesMapGet(linearMapping,linearVariableIdx, &
    & variableToEquationsMatricesMap,err,error,*)

    !Argument variables
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping !<A pointer to the linear mapping to get the variable to equations matrices map for
    INTEGER(INTG), INTENT(IN) :: linearVariableIdx !<The linear variable index to get the linear variable to equations matrices map for
    TYPE(varToEquationsMatricesMapType), POINTER :: variableToEquationsMatricesMap !<On exit, a pointer to the requested linear variable to equations matrices map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("EquationsMappingLinear_VariableToEquationsMatricesMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(variableToEquationsMatricesMap)) &
      & CALL FlagError("Variable to equations matrices map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(linearMapping)) CALL FlagError("Linear mapping is not associated.",err,error,*999)
    IF(linearVariableIdx<1.OR.linearVariableIdx>linearMapping%numberOfLinearVariables) THEN
      localError="The specified linear variable index of "//TRIM(NumberToVString(linearVariableIdx,"*",err,error))// &
        & " is invalid. The linear variable index should be >= 1 and <= " &
        & TRIM(NumberToVString(linearMapping%numberOfLinearVariables,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(linearMapping%varToEquationsMatricesMaps)) &
      & CALL FlagError("The variable to equations matrices maps is not allocated for the linear mapping.",err,error,*999)
#endif    
    
    variableToEquationsMatricesMap=>linearMapping%varToEquationsMatricesMaps(linearVariableIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(variableToEquationsMatricesMap)) THEN
      localError="The variable to equations matrices map is not associated for linear variable index "// &
        & TRIM(NumberToVString(linearVariableIdx,"*",err,error)//" of the linear mapping."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    EXITS("EquationsMappingLinear_VariableToEquationsMatricesMapGet")
    RETURN
999 NULLIFY(variableToEquationsMatricesMap)
998 ERRORS("EquationsMappingLinear_VariableToEquationsMatricesMapGet",err,error)
    EXITS("EquationsMappingLinear_VariableToEquationsMatricesMapGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingLinear_VariableToEquationsMatricesMapGet

  !
  !================================================================================================================================
  !

  !>Gets the number of residuals in a nonlinear mapping.
  SUBROUTINE EquationsMappingNonlinear_NumberOfResidualsGet(nonlinearMapping,numberOfResiduals,err,error,*)

    !Argument variables
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping !<A pointer to the nonlinear mapping to get the number of residuals for
    INTEGER(INTG), INTENT(OUT) :: numberOfResiduals !<On exit, the number of residuals in this nonlinear mapping.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingNonlinear_NumberOfResidualsGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(nonlinearMapping)) CALL FlagError("Nonlinear mapping is not associated.",err,error,*999)
#endif    
    
    numberOfResiduals=nonlinearMapping%numberOfResiduals
   
    EXITS("EquationsMappingNonlinear_NumberOfResidualsGet")
    RETURN
999 ERRORS("EquationsMappingNonlinear_NumberOfResidualsGet",err,error)
    EXITS("EquationsMappingNonlinear_NumberOfResidualsGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingNonlinear_NumberOfResidualsGet

  !
  !================================================================================================================================
  !

  !>Gets the residual mapping for a nonlinear mapping.
  SUBROUTINE EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,residualIdx,residualMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping !<A pointer to the nonlinear mapping to get the residual mapping for
    INTEGER(INTG), INTENT(IN) :: residualIdx !<The resididual index to get the residual mapping for
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping !<On exit, a pointer to the requested residual mapping. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("EquationsMappingNonlinear_ResidualMappingGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(residualMapping)) &
      & CALL FlagError("Residual mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(nonlinearMapping)) CALL FlagError("Nonlinear mapping is not associated.",err,error,*999)
    IF(residualIdx<1.OR.residualIdx>nonlinearMapping%numberOfResiduals) THEN
      localError="The specified residual index of "//TRIM(NumberToVString(residualIdx,"*",err,error))// &
        & " is invalid. The residual index should be >= 1 and <= " &
        & TRIM(NumberToVString(nonlinearMapping%numberOfResiduals,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(nonlinearMapping%residuals)) &
      & CALL FlagError("The residuals is not allocated for the nonlinear mapping.",err,error,*999)
#endif    
    
    residualMapping=>nonlinearMapping%residuals(residualIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(residualMapping)) THEN
      localError="The residual mapping is not associated for residual index "// &
        & TRIM(NumberToVString(residualIdx,"*",err,error)//" of the nonlinear mapping."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    EXITS("EquationsMappingNonlinear_ResidualMappingGet")
    RETURN
999 NULLIFY(residualMapping)
998 ERRORS("EquationsMappingNonlinear_ResidualMappingGet",err,error)
    EXITS("EquationsMappingNonlinear_ResidualMappingGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingNonlinear_ResidualMappingGet

  !
  !================================================================================================================================
  !

  !>Gets the Jacobian matrix to variable mapping in a residual mapping.
  SUBROUTINE EquationsMappingResidual_JacobianMatrixToVariableMapGet(residualMapping,jacobianMatrixIdx, &
    & jacobianMatrixToVariableMap,err,error,*)

    !Argument variables
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping !<A pointer to the residual mapping to get the Jacobian matrix to variable map for
    INTEGER(INTG), INTENT(IN) :: jacobianMatrixIdx !<The Jacobian matrix index of the residual to get the Jacobian matrix to variable map for
    TYPE(JacobianMatrixToVarMapType), POINTER :: jacobianMatrixToVariableMap !<On exit, a pointer to the requested Jacobian matrix to variable map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif
    
    ENTERS("EquationsMappingResidual_JacobianMatrixToVariableMap",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(jacobianMatrixToVariableMap)) &
      & CALL FlagError("Jacobian matrix to variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(residualMapping)) CALL FlagError("Residual mapping is not associated.",err,error,*999)
    IF(jacobianMatrixIdx<1.OR.jacobianMatrixIdx>linearMapping%numberOfJacobianMatrices) THEN
      localError="The specified Jacobian matrix index of "//TRIM(NumberToVString(jacobianMatrixIdx,"*",err,error))// &
        & " of residual number "//TRIM(NumberToVString(residualMapping%residualNumber,"*",err,error))// &
        & " is invalid. The Jacobian matrix index should be >=1 and <= "// &
        & TRIM(NumberToVString(residualMapping%numberOfJacobianMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(residualMapping%jacobianMatrixToVarMaps)) THEN
      localError="The Jacobian matrix to variable maps is not allocated for the residual mapping for residual number "// &
        & TRIM(NumberToVString(residualMapping%residualNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    jacobianMatrixToVariableMap=>residualMapping%jacobianMatrixToVarMaps(jacobianMatrixIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(jacobianMatrixToVariableMap)) THEN
      localError="The Jacobian matrix to variable map is not associated for Jacobian matrix index "// &
        & TRIM(NumberToVString(linearMatrixIdx,"*",err,error))//" of the mapping for residual number "// &
        & TRIM(NumberToVString(residualMapping,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    EXITS("EquationsMappingResidual_JacobianMatrixToVariableMapGet")
    RETURN
999 NULLIFY(jacobianMatrixToVariableMap)
998 ERRORS("EquationsMappingResidual_JacobianMatrixToVariableMapGet",err,error)
    EXITS("EquationsMappingResidual_JacobianMatrixToVariableMapGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingResidual_JacobianMatrixToVariableMapGet
  
  !
  !================================================================================================================================
  !

  !>Gets the specified Jacobian matrix variable for a residual mapping.
  SUBROUTINE EquationsMappingResidual_JacobianMatrixVariableGet(residualMapping,jacobianMatrixIdx,jacobianMatrixVariable, &
    & err,error,*)

    !Argument variables
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping !<A pointer to the residual mapping to get the Jacobian matrix variable for
    INTEGER(INTG), INTENT(IN) :: jacobianMatrixIdx !<The index of the Jacobian matrix to get the Jacobian variable for. 
    TYPE(FieldVariableType), POINTER :: jacobianMatrixVariable !<On exit, a pointer to the requested Jacobian matrix variable. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(JacobianMatrixToVarMapType), POINTER :: jacobianMatrixToVariableMap
#endif
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("EquationsMappingResidual_JacobianMatrixVariableGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(jacobianMatrixVariable)) CALL FlagError("Jacobian matrix variable is already associated.",err,error,*998)
    NULLIFY(jacobianMatrixToVariableMap)
    CALL EquationsMappingResidual_JacobianMatrixToVarMapGet(residualMapping,jacobianMatrixIdx,jacobianMatrixToVarMap,err,error,*999)
#endif    
   
    jacobianMatrixVariable=>residualMapping%jacobianMatrixToVarMaps(jacobianMatrixIdx)%ptr%variable

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(jacobianMatrixVariable)) THEN
      localError="The Jacobian matrix variable is not associated for Jacobian matrix index "// &
        & TRIM(NumberToVString(jacobianMatrixIdx,"*",err,error))//" of residual number "// &
        & TRIM(NumberToVString(residualMapping%residualNumber,"*",err,error)//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    EXITS("EquationsMappingResidual_JacobianMatrixVariableGet")
    RETURN
999 NULLIFY(jacobianMatrixVariable)
998 ERRORS("EquationsMappingResidual_JacobianMatrixVariableGet",err,error)
    EXITS("EquationsMappingResidual_JacoibanMatrixVariableGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingResidual_JacobianMatrixVariableGet

  !
  !================================================================================================================================
  !

  !>Gets the specified Jacobian matrix variable type for a residual mapping.
  SUBROUTINE EquationsMappingResidual_JacobianMatrixVariableTypeGet(residualMapping,jacobianMatrixIdx,jacobianMatrixVariableType, &
    & err,error,*)

    !Argument variables
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping !<A pointer to the residual mapping to get the Jacobian matrix variable type for
    INTEGER(INTG), INTENT(IN) :: jacobianMatrixIdx !<The index of the Jacobian matrix to get the Jacobian variable type for. 
    INTEGER(INTG), INTENT(OUT) :: jacobianMatrixVariableType !<On exit, the requested Jacobian matrix variable type.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(JacobianMatrixToVarMapType), POINTER :: jacobianMatrixToVariableMap
#endif
 
    ENTERS("EquationsMappingResidual_JacobianMatrixVariableTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    NULLIFY(jacobianMatrixToVariableMap)
    CALL EquationsMappingResidual_JacobianMatrixToVarMapGet(residualMapping,jacobianMatrixIdx,jacobianMatrixToVarMap,err,error,*999)
#endif    
   
    jacobianMatrixVariableType=residualMapping%jacobianMatrixToVarMaps(jacobianMatrixIdx)%ptr%variableType

    EXITS("EquationsMappingResidual_JacobianMatrixVariableTypeGet")
    RETURN
999 ERRORS("EquationsMappingResidual_JacobianMatrixVariableTypeGet",err,error)
    EXITS("EquationsMappingResidual_JacoibanMatrixVariableGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingResidual_JacobianMatrixVariableTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the number of Jacobian matrices for a residual mapping.
  SUBROUTINE EquationsMappingResidual_NumberOfJacobianMatricesGet(residualMapping,numberOfJacobianMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping !<A pointer to the residual mapping to get the number of Jacobian matrices for
    INTEGER(INTG), INTENT(OUT) :: numberOfJacobianMatrices !<On exit, the number of Jacobians matrices for the residual mapping.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingResidual_NumberOfJacobianMatricesGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(residualMapping)) CALL FlagError("Residual mapping is not associated.",err,error,*999)
#endif    
      
    numberOfJacobianMatrices=residualMapping%numberOfJacobianMatrices
   
    EXITS("EquationsMappingResidual_NumberOfJacobianMatricesGet")
    RETURN
999 ERRORS("EquationsMappingResidual_NumberOfJacobianMatricesGet",err,error)
    EXITS("EquationsMappingResidual_NumberOfJacobianMatricesGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingResidual_NumberOfJacobianMatricesGet

  !
  !================================================================================================================================
  !

  !>Gets the number of residual variables for a residual mapping.
  SUBROUTINE EquationsMappingResidual_NumberOfResidualVariablesGet(residualMapping,numberOfResidualVariables,err,error,*)

    !Argument variables
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping !<A pointer to the residual mapping to get the number of residual variables for
    INTEGER(INTG), INTENT(OUT) :: numberOfResidualVariables !<On exit, the number of residual variables for the residual mapping.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMappingResidual_NumberOfResidualVariablesGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(residualMapping)) CALL FlagError("Residual mapping is not associated.",err,error,*999)
#endif    
      
    numberOfResidualVariables=residualMapping%numberOfResidualVariables
   
    EXITS("EquationsMappingResidual_NumberOfResidualVariablesGet")
    RETURN
999 ERRORS("EquationsMappingResidual_NumberOfResidualVariablesGet",err,error)
    EXITS("EquationsMappingResidual_NumberOfResidualVariablesGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingResidual_NumberOfResidualVariablesGet

  !
  !================================================================================================================================
  !

  !>Gets the specified residual variable for a residual mapping.
  SUBROUTINE EquationsMappingResidual_VariableGet(residualMapping,variableIdx,fieldVariable,err,error,*)

    !Argument variables
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping !<A pointer to the residual mapping to get the residual variable for
    INTEGER(INTG), INTENT(IN) :: variableIdx !<The index of the variable in the residual to get the field variable for. 
    TYPE(FieldVariableType), POINTER :: fieldVariable !<On exit, a pointer to the requested field variable. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("EquationsMappingNonlinear_VariableGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(residualMapping)) CALL FlagError("Residual mapping is not associated.",err,error,*999)
    IF(variableIdx<1.OR.variableIdx>residualMapping%numberOfResidualVariables) THEN
      localError="The specified residual variable index of "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
        & " is invalid for residual number "//TRIM(NumberToVString(residualMapping%residualNumber,"*",err,error))// &
        & " which has "//TRIM(NumberToVString(residualMapping%numberOfVariables,"*",err,error))//" variables."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(residualMapping%variables)) THEN
      localError="Variables is not allocated for the residual mapping for residual number "// &
        & TRIM(NumberToVString(residualMapping%residualNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    fieldVariable=>residualMapping%variables(variableIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(fieldVariable)) THEN
      localError="The field variable is not associated for variable index "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
        & " of residual number "//TRIM(NumberToVString(residualMapping%residualNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    EXITS("EquationsMappingResidual_VariableGet")
    RETURN
999 NULLIFY(fieldVariable)
998 ERRORS("EquationsMappingResidual_VariableGet",err,error)
    EXITS("EquationsMappingResidual_VariableGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingResidual_VariableGet

  !
  !================================================================================================================================
  !

  !>Gets the specified variable index for a variable type in a residual mapping.
  SUBROUTINE EquationsMappingResidual_VariableIndexGet(residualMapping,variableType,variableIdx,err,error,*)

    !Argument variables
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping !<A pointer to the residual mapping to get the variable index for
    INTEGER(INTG), INTENT(IN) :: variableType !<The type of the variable to get the variable index for. 
    INTEGER(INTG), INTENT(OUT) :: variableIdx !<On exit, the requested variable index corresponding to the variable type for the residual mapping. If there is no variable type the returned value will be zero.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("EquationsMappingResidual_VariableIndexGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(residualMapping)) CALL FlagError("Residual mapping is not associated.",err,error,*999)    
    IF(variableType<1.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified variable type of "//TRIM(NumberToVString(variabletType,"*",err,error))// &
        & " is invalid. The variable type should be >=1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(residualMapping%variableTypesMap)) THEN
      localError="Residual mapping variables type map is not allocated for residual number "// &
        & TRIM(NumberToVString(residualMapping%residualNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    variableIdx=residualMapping%variableTypeMap(variableType)
   
    EXITS("EquationsMappingResidual_VariableIndexGet")
    RETURN
999 ERRORS("EquationsMappingResidual_VariableIndexGet",err,error)
    EXITS("EquationsMappingResidual_VariableIndexGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingLinear_LinearVariableIndexGet

  !
  !================================================================================================================================
  !

  !>Gets the variable to Jacobian matrix map for the residual variable in a residual mapping.
  SUBROUTINE EquationsMappingResidual_VariableToJacobianMatrixMapGet(residualMapping,variableIdx, &
    & variableToJacobianMatrixMap,err,error,*)
    
    !Argument variables
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping !<A pointer to the residual mapping to get the variable to Jacobian matrix map for
    INTEGER(INTG), INTENT(IN) :: variableIdx !<The residual variable index to get the variable to Jacobian matrices map for
    TYPE(varToJacobianMatrixMapType), POINTER :: variableToJacobianMatrixMap !<On exit, a pointer to the requested variable to Jacobain matrix map. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("EquationsMappingResidual_VariableToJacobianMatrixMapGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(variableToJacobianMatrixMap)) &
      & CALL FlagError("Variable to Jacobian matrix map is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(residualMapping)) CALL FlagError("Residual mapping is not associated.",err,error,*999)
    IF(variableIdx<1.OR.variableIdx>residualMapping%numberOfVariables) THEN
      localError="The specified variable index of "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
        & " is invalid for residual number "//TRIM(NumberToVString(residualMapping%residualNumber))// &
        & ". The variable index should be >= 1 and <= "// &
        & TRIM(NumberToVString(residualMapping%numberOfVariables,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(residualMapping%varToJacobianMatrixsMaps)) THEN
      localError="The variable to Jacobian matrix maps is not allocated for the residual mapping for residual number "// &
        & TRIM(NumberToVString(residualMapping%residualNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    variableToJacobianMatrixMap=>residualMapping%varToJacobianMatrixMaps(variableIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(variableToJacobianMatrixMap)) THEN
      localError="The variable to Jacobian matrix map is not associated for variable index "// &
        & TRIM(NumberToVString(variableIdx,"*",err,error)//" of the mapping for residual number "// &
        & TRIM(NumberToVString(residualMapping%residualNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    EXITS("EquationsMappingResidual_VariableToJacobianMatrixMapGet")
    RETURN
999 NULLIFY(variableToJacobianMatrixMap)
998 ERRORS("EquationsMappingResidual_VariableToJacobianMatrixMapGet",err,error)
    EXITS("EquationsMappingResidual_VariableToJacobianMatrixMapGet")
    RETURN 1
e    
  END SUBROUTINE EquationsMappingResidual_VariableToJacobianMatrixMapGet
  
  !
  !================================================================================================================================
  !

  !>Gets the specified variable type for a residual mapping.
  SUBROUTINE EquationsMappingResidual_VariableTypeGet(residualMapping,variableIdx,variableType,err,error,*)

    !Argument variables
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping !<A pointer to the residual mapping to get the variable type for
    INTEGER(INTG), INTENT(IN) :: variableIdx !<The index of the variable in the residual to get the field variable for. 
    INTEGER(INTG), INTENT(OUT) :: variableType !<On exit, a pointer to the requested residual variable type.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("EquationsMappingResidual_VariableTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(residualMapping)) CALL FlagError("Residual mapping is not associated.",err,error,*999)
    IF(variableIdx<1.OR.variableIdx>residualMapping%numberOfVariables) THEN
      localError="The specified variable index of "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
        & " is invalid for residual number "//TRIM(NumberToVString(residualMapping%residualNumber,"*",err,error))// &
        & " which has "//TRIM(NumberToVString(residualMapping%numberOfVariables,"*",err,error))//" variables."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(residualMapping%variableTypes)) THEN
      localError="Residual variable types is not allocated for residual number "// &
        & TRIM(NumberToVString(residualMappping%residualNumber,err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
      
    variableType=residualMapping%variableTypes(variableIdx)
   
    EXITS("EquationsMappingResidual_VariableTypeGet")
    RETURN
999 ERRORS("EquationsMappingResidual_VariableTypeGet",err,error)
    EXITS("EquationsMappingResidual_VariableTypeGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingResidual_VariableTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the specified RHS variable for a rhs mapping.
  SUBROUTINE EquationsMappingRHS_RHSVariableGet(rhsMapping,rhsVariable,err,error,*)

    !Argument variables
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping !<A pointer to the RHS mapping to get the RHS variable for
    TYPE(FieldVariableType), POINTER :: rhsVariable !<On exit, a pointer to the requested RHS field variable. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingRHS_RHSVariableGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(rhsVariable)) CALL FlagError("RHS Field variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(rhsMapping)) CALL FlagError("RHS mapping is not associated.",err,error,*999)
#endif    
    
    rhsVariable=>rhsMapping%rhsVariable

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(rhsVariable)) CALL FlagError("The RHS field variable is not associated in the RHS mapping.",err,error,*999)
#endif    
    
    EXITS("EquationsMappingRHS_RHSVariableGet")
    RETURN
999 NULLIFY(rhsVariable)
998 ERRORS("EquationsMappingRHS_RHSVariableGet",err,error)
    EXITS("EquationsMappingRHS_RHSVariableGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingRHS_RHSVariableGet

  !
  !================================================================================================================================
  !

  !>Gets the specified RHS variable type for a rhs mapping.
  SUBROUTINE EquationsMappingRHS_RHSVariableTypeGet(rhsMapping,rhsVariableType,err,error,*)

    !Argument variables
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping !<A pointer to the RHS mapping to get the RHS variable for
    INTEGER(INTG), INTENT(OUT) :: rhsVariableType !<On exit, the requested RHS variable type.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingRHS_RHSVariableTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(rhsMapping)) CALL FlagError("RHS mapping is not associated.",err,error,*999)
#endif    
    
    rhsVariableType=rhsMapping%rhsVariableType
    
    EXITS("EquationsMappingRHS_RHSVariableTypeGet")
    RETURN
999 ERRORS("EquationsMappingRHS_RHSVariableTypeGet",err,error)
    EXITS("EquationsMappingRHS_RHSVariableTypeGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingRHS_RHSVariableTypeGet

  !
  !=================================================================================================================================
  !

  !>Assert that a scalar equations mapping has been finished
  SUBROUTINE EquationsMappingScalar_AssertIsFinished(scalarMapping,err,error,*)

    !Argument Variables
    TYPE(EquationsMappingScalarType), POINTER, INTENT(IN) :: scalarMapping !<The scalar mapping to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingScalar_AssertIsFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(scalarMapping)) CALL FlagError("Scalar mapping is not associated.",err,error,*999)
#endif    

    IF(.NOT.scalarMapping%scalarMappingFinished) &
      & CALL FlagError("Scalar equations mapping has not been finished."
    
    EXITS("EquationsMappingScalar_AssertIsFinished")
    RETURN
999 ERRORSEXITS("EquationsMappingScalar_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingScalar_AssertIsFinished

  !
  !=================================================================================================================================
  !

  !>Assert that a scalar equations mapping has not been finished
  SUBROUTINE EquationsMappingScalar_AssertNotFinished(scalarMapping,err,error,*)

    !Argument Variables
    TYPE(EquationsMappingScalarType), POINTER, INTENT(IN) :: scalarMapping !<The scalar mapping to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingScalar_AssertNotFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(scalarMapping)) CALL FlagError("Scalar mapping is not associated.",err,error,*999)
#endif    

    IF(scalarMapping%scalarMappingFinished) &
      & CALL FlagError("Scalar equations mapping has already been finished.",err,error,*999)
    
    EXITS("EquationsMappingScalar_AssertNotFinished")
    RETURN
999 ERRORSEXITS("EquationsMappingScalar_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingScalar_AssertNotFinished
  
  !
  !================================================================================================================================
  !

  !>Gets the create values cache for an scalar equations mapping.
  SUBROUTINE EquationsMappingScalar_CreateValuesCacheGet(scalarMapping,createValuesCache,err,error,*)

    !Argument variables
    TYPE(EquationsMappingScalarType), POINTER :: scalarMapping !<A pointer to the equations scalar mapping to get the create values cache for
    TYPE(EquationsMappingScalarCreateValuesCacheType), POINTER :: createValuesCache !<On exit, a pointer to the create values cache in the specified scalar equations mapping. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingScalar_CreateValuesCacheGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(createValuesCache)) CALL FlagError("Create values cache is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(scalarMapping)) CALL FlagError("Scalar mapping is not associated.",err,error,*999)
#endif    

    createValuesCache=>scalarMapping%createValuesCache

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(createValuesCache)) &
      & CALL FlagError("Create values cache is not associated for the scalar mapping.",err,error,*999)
#endif    
       
    EXITS("EquationsMappingScalar_CreateValuesCacheGet")
    RETURN
999 NULLIFY(createValuesCache)
998 ERRORSEXITS("EquationsMappingScalar_CreateValuesCacheGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingScalar_CreateValuesCacheGet

  !
  !================================================================================================================================
  !

  !>Gets the scalar equations for an scalar equations mapping.
  SUBROUTINE EquationsMappingScalar_ScalarEquationsGet(scalarMapping,scalarEquations,err,error,*)

    !Argument variables
    TYPE(EquationsMappingScalarType), POINTER :: scalarMapping !<A pointer to the equations scalar mapping to get the scalar equations for
    TYPE(EquationsScalarType), POINTER :: scalarEquations !<On exit, a pointer to the scalar equations in the specified scalar equations mapping. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingScalar_ScalarEquationsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(scalarEquations)) CALL FlagError("Scalar equations is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(scalarMapping)) CALL FlagError("Scalar mapping is not associated.",err,error,*999)
#endif    

    scalarEquations=>scalarMapping%scalarEquations

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(scalarEquations)) CALL FlagError("Scalar equations is not associated for the scalar mapping.",err,error,*999)
#endif    
       
    EXITS("EquationsMappingScalar_ScalarEquationsGet")
    RETURN
999 NULLIFY(scalarEquations)
998 ERRORSEXITS("EquationsMappingScalar_ScalarEquationsGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingScalar_ScalarEquationsGet

  !
  !================================================================================================================================
  !

  !>Gets the specified source variable for a source mapping.
  SUBROUTINE EquationsMappingSource_SourceVariableGet(sourceMapping,sourceVariable,err,error,*)

    !Argument variables
    TYPE(EquationsMappingSourceType), POINTER :: sourceMapping !<A pointer to the source mapping to get the source variable for
    TYPE(FieldVariableType), POINTER :: sourceVariable !<On exit, a pointer to the requested source variable. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingSource_SourceVariableGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(sourceVariable)) CALL FlagError("Source variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(sourceMapping)) CALL FlagError("Source mapping is not associated.",err,error,*999)
#endif    
    
    sourceVariable=>sourceMapping%sourceVariable

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(sourceVariable)) CALL FlagError("The source variable is not associated.",err,error,*999)
#endif    
    
    EXITS("EquationsMappingSource_SourceVariableGet")
    RETURN
999 NULLIFY(sourceVariable)
998 ERRORS("EquationsMappingSource_SourceVariableGet",err,error)
    EXITS("EquationsMappingSource_SourceVariableGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingSource_SourceVariableGet

  !
  !================================================================================================================================
  !

  !>Gets the specified source mapping for a sources mapping.
  SUBROUTINE EquationsMappingSources_SourceMappingGet(sourcesMapping,sourceIdx,sourceMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingSourcesType), POINTER :: sourcesMapping !<A pointer to the sources mapping to get the source mapping for
    INTEGER(INTG) :: sourceIdx !<The source index to get the source mapping for
    TYPE(EquationsMappingSourceType), POINTER :: sourceMapping !<On exit, a pointer to the requested source mapping. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif 
 
    ENTERS("EquationsMappingSources_SourceMappingGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(sourceMapping)) CALL FlagError("Source mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(sourcesMapping)) CALL FlagError("Sources mapping is not associated.",err,error,*999)
    IF(sourceIdx<1.OR.sourceIdx>sourcesMapping%numberOfSources) THEN
      localError="The specified source index of "//TRIM(NumberToVString(sourceIdx,"*",err,error))// &
        & " is invalid. The source index should be >= 1 and <= "// &
        & TRIM(NumberToVString(sourcesMapping%numberOfSources,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(sourcesMapping%sources)) &
      & CALL FlagError("The sources array is not allocated for the sources mapping.",err,error,*999)
#endif    
    
    sourceMapping=>sourcesMapping%sources(sourceIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(sourceMapping)) THEN
      localError="The source mapping is not associated for source index "//TRIM(NumberToVString(sourceIdx,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
    
    EXITS("EquationsMappingSources_SourceMappingGet")
    RETURN
999 NULLIFY(sourceMapping)
998 ERRORS("EquationsMappingSources_SourceMappingGet",err,error)
    EXITS("EquationsMappingSources_SourceMappingGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingSources_SourceMappingGet

  !
  !=================================================================================================================================
  !

  !>Assert that a vector equations mapping has been finished
  SUBROUTINE EquationsMappingVector_AssertIsFinished(vectorMapping,err,error,*)

    !Argument Variables
    TYPE(EquationsMappingVectorType), POINTER, INTENT(IN) :: vectorMapping !<The vector mapping to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingVector_AssertIsFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector mapping is not associated.",err,error,*999)
#endif    

    IF(.NOT.vectorMapping%vectorMappingFinished) &
      & CALL FlagError("Vector equations mapping has not been finished."
    
    EXITS("EquationsMappingVector_AssertIsFinished")
    RETURN
999 ERRORSEXITS("EquationsMappingVector_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_AssertIsFinished

  !
  !=================================================================================================================================
  !

  !>Assert that a vector equations mapping has not been finished
  SUBROUTINE EquationsMappingVector_AssertNotFinished(vectorMapping,err,error,*)

    !Argument Variables
    TYPE(EquationsMappingVectorType), POINTER, INTENT(IN) :: vectorMapping !<The vector mapping to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingVector_AssertNotFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector mapping is not associated.",err,error,*999)
#endif    

    IF(vectorMapping%vectorMappingFinished) &
      & CALL FlagError("Vector equations mapping has already been finished.",err,error,*999)
    
    EXITS("EquationsMappingVector_AssertNotFinished")
    RETURN
999 ERRORSEXITS("EquationsMappingVector_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_AssertNotFinished
  
  !
  !================================================================================================================================
  !

  !>Gets the create values cache for an vector equations mapping.
  SUBROUTINE EquationsMappingVector_CreateValuesCacheGet(vectorMapping,createValuesCache,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations vector mapping to get the create values cache for
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache !<On exit, a pointer to the create values cache in the specified vector equations mapping. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingVector_CreateValuesCacheGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(createValuesCache)) CALL FlagError("Create values cache is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector mapping is not associated.",err,error,*999)
#endif    

    createValuesCache=>vectorMapping%createValuesCache

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(createValuesCache)) &
      & CALL FlagError("Create values cache is not associated for the vector mapping.",err,error,*999)
#endif    
       
    EXITS("EquationsMappingVector_CreateValuesCacheGet")
    RETURN
999 NULLIFY(createValuesCache)
998 ERRORSEXITS("EquationsMappingVector_CreateValuesCacheGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_CreateValuesCacheGet

  !
  !================================================================================================================================
  !

  !>Checks if the dynamic vector mapping for an vector mapping exists.
  SUBROUTINE EquationsMappingVector_DynamicMappingExists(vectorMapping,dynamicMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations vector mapping to get the dynamic mapping for
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping !<On exit, a pointer to the dynamic Mapping in the specified vector equations mapping. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingVector_DynamicMappingExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dynamicMapping)) CALL FlagError("Dynamic mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector mapping is not associated.",err,error,*999)
#endif    

    dynamicMapping=>vectorMapping%dynamicMapping
       
    EXITS("EquationsMappingVector_DynamicMappingExists")
    RETURN
999 NULLIFY(dynamicMapping)
998 ERRORSEXITS("EquationsMappingVector_DynamicMappingExists",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_DynamicMappingExists

  !
  !================================================================================================================================
  !

  !>Gets the dynamic vector mapping for an vector mapping.
  SUBROUTINE EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations vector mapping to get the dynamic mapping for
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping !<On exit, a pointer to the dynamic Mapping in the specified vector equations mapping. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingVector_DynamicMappingGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dynamicMapping)) CALL FlagError("Dynamic mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector mapping is not associated.",err,error,*999)
#endif    

    dynamicMapping=>vectorMapping%dynamicMapping

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(dynamicMapping)) CALL FlagError("Dynamic mapping is not associated for the vector mapping.",err,error,*999)
#endif    
       
    EXITS("EquationsMappingVector_DynamicMappingGet")
    RETURN
999 NULLIFY(dynamicMapping)
998 ERRORSEXITS("EquationsMappingVector_DynamicMappingGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_DynamicMappingGet

  !
  !================================================================================================================================
  !

  !>Gets the LHS vector mapping for an vector mapping.
  SUBROUTINE EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations vector mapping to get the LHS mapping for
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping !<On exit, a pointer to the LHS Mapping in the specified vector equations mapping. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingVector_LHSMappingGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(lhsMapping)) CALL FlagError("LHS mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector mapping is not associated.",err,error,*999)
#endif    

    lhsMapping=>vectorMapping%lhsMapping

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(lhsMapping)) CALL FlagError("LHS mapping is not associated for the vector mapping.",err,error,*999)
#endif    
       
    EXITS("EquationsMappingVector_LHSMappingGet")
    RETURN
999 NULLIFY(lhsMapping)
998 ERRORSEXITS("EquationsMappingVector_LHSMappingGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_LHSMappingGet

  !
  !================================================================================================================================
  !

  !>Checks if the linear vector mapping for an vector mapping exists.
  SUBROUTINE EquationsMappingVector_LinearMappingExists(vectorMapping,linearMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations vector mapping to get the linear mapping for
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping !<On exit, a pointer to the linear Mapping in the specified vector equations mapping. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingVector_LinearMappingExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(linearMapping)) CALL FlagError("Linear mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector mapping is not associated.",err,error,*999)
#endif    

    linearMapping=>vectorMapping%linearMapping
       
    EXITS("EquationsMappingVector_LinearMappingExists")
    RETURN
999 NULLIFY(linearMapping)
998 ERRORSEXITS("EquationsMappingVector_LinearMappingExists",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_LinearMappingExists

  !
  !================================================================================================================================
  !

  !>Gets the linear vector mapping for an vector mapping.
  SUBROUTINE EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations vector mapping to get the linear mapping for
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping !<On exit, a pointer to the linear Mapping in the specified vector equations mapping. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingVector_LinearMappingGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(linearMapping)) CALL FlagError("Linear mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector mapping is not associated.",err,error,*999)
#endif    

    linearMapping=>vectorMapping%linearMapping

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(linearMapping)) CALL FlagError("Linear mapping is not associated for the vector mapping.",err,error,*999)
#endif    
       
    EXITS("EquationsMappingVector_LinearMappingGet")
    RETURN
999 NULLIFY(linearMapping)
998 ERRORSEXITS("EquationsMappingVector_LinearMappingGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_LinearMappingGet

  !
  !================================================================================================================================
  !

  !>Checks if the nonlinear vector mapping for an vector mapping exists.
  SUBROUTINE EquationsMappingVector_NonlinearMappingExists(vectorMapping,nonlinearMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations vector mapping to get the nonlinear mapping for
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping !<On exit, a pointer to the nonlinear Mapping in the specified vector equations mapping. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingVector_NonlinearMappingExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(nonlinearMapping)) CALL FlagError("Nonlinear mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector mapping is not associated.",err,error,*999)
#endif    

    nonlinearMapping=>vectorMapping%nonlinearMapping
       
    EXITS("EquationsMappingVector_NonlinearMappingExists")
    RETURN
999 NULLIFY(nonlinearMapping)
998 ERRORS("EquationsMappingVector_NonlinearMappingExists",err,error)
    EXITS("EquationsMappingVector_NonlinearMappingExists")
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_NonlinearMappingExists

  !
  !================================================================================================================================
  !

  !>Gets the nonlinear vector mapping for an vector mapping.
  SUBROUTINE EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations vector mapping to get the nonlinear mapping for
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping !<On exit, a pointer to the nonlinear Mapping in the specified vector equations mapping. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingVector_NonlinearMappingGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(nonlinearMapping)) CALL FlagError("Nonlinear mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector mapping is not associated.",err,error,*999)
#endif    

    nonlinearMapping=>vectorMapping%nonlinearMapping

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(nonlinearMapping))  &
      & CALL FlagError("Nonlinear mapping is not associated for the vector mapping.",err,error,*999)
#endif    
       
    EXITS("EquationsMappingVector_NonlinearMappingGet")
    RETURN
999 NULLIFY(nonlinearMapping)
998 ERRORSEXITS("EquationsMappingVector_NonlinearMappingGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_NonlinearMappingGet

  !
  !================================================================================================================================
  !

  !>Checks if the rhs vector mapping for an vector mapping exists.
  SUBROUTINE EquationsMappingVector_RHSMappingExists(vectorMapping,rhsMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations vector mapping to get the rhs mapping for
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping !<On exit, a pointer to the rhs Mapping in the specified vector equations mapping. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingVector_RHSMappingExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(rhsMapping)) CALL FlagError("RHS mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector mapping is not associated.",err,error,*999)
#endif    

    rhsMapping=>vectorMapping%rhsMapping
       
    EXITS("EquationsMappingVector_RHSMappingExists")
    RETURN
999 NULLIFY(rhsMapping)
998 ERRORSEXITS("EquationsMappingVector_RHSMappingExists",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_RHSMappingExists

  !
  !================================================================================================================================
  !

  !>Gets the rhs vector mapping for an vector mapping.
  SUBROUTINE EquationsMappingVector_RHSMappingGet(vectorMapping,rhsMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations vector mapping to get the rhs mapping for
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping !<On exit, a pointer to the rhs Mapping in the specified vector equations mapping. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingVector_RHSMappingGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(rhsMapping)) CALL FlagError("RHS mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector mapping is not associated.",err,error,*999)
#endif    

    rhsMapping=>vectorMapping%rhsMapping

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(rhsMapping)) CALL FlagError("RHS mapping is not associated for the vector mapping.",err,error,*999)
#endif    
       
    EXITS("EquationsMappingVector_RHSMappingGet")
    RETURN
999 NULLIFY(rhsMapping)
998 ERRORSEXITS("EquationsMappingVector_RHSMappingGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_RHSMappingGet

  !
  !================================================================================================================================
  !

  !>Checks if the sources vector mapping for an vector mapping exists.
  SUBROUTINE EquationsMappingVector_SourcesMappingExists(vectorMapping,sourcesMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations vector mapping to get the source mapping for
    TYPE(EquationsMappingSourcesType), POINTER :: sourcesMapping !<On exit, a pointer to the source Mapping in the specified vector equations mapping. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingVector_SourcesMappingExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(sourcseMapping)) CALL FlagError("Sources mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector mapping is not associated.",err,error,*999)
#endif    

    sourcesMapping=>vectorMapping%sourcesMapping
       
    EXITS("EquationsMappingVector_SourcesMappingExists")
    RETURN
999 NULLIFY(sourcesMapping)
998 ERRORSEXITS("EquationsMappingVector_SourcesMappingExists",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_SourcesMappingExists

  !
  !================================================================================================================================
  !

  !>Gets the sources vector mapping for an vector mapping.
  SUBROUTINE EquationsMappingVector_SourcesMappingGet(vectorMapping,sourcesMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations vector mapping to get the source mapping for
    TYPE(EquationsMappingSourcesType), POINTER :: sourcesMapping !<On exit, a pointer to the sources mapping in the specified vector equations mapping. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingVector_SourcesMappingGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(sourcesMapping)) CALL FlagError("Sources mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector mapping is not associated.",err,error,*999)
#endif    

    sourcesMapping=>vectorMapping%sourcesMapping

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(sourcesMapping)) CALL FlagError("Sources mapping is not associated for the vector mapping.",err,error,*999)
#endif    
       
    EXITS("EquationsMappingVector_SourcesMappingGet")
    RETURN
999 NULLIFY(sourceMapping)
998 ERRORSEXITS("EquationsMappingVector_SourcesMappingGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_SourcesMappingGet

  !
  !================================================================================================================================
  !

  !>Gets the vector equations for an vector equations mapping.
  SUBROUTINE EquationsMappingVector_VectorEquationsGet(vectorMapping,vectorEquations,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations vector mapping to get the vector equations for
    TYPE(EquationsVectorType), POINTER :: vectorEquations !<On exit, a pointer to the vector equations in the specified vector equations mapping. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingVector_VectorEquationsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(vectorEquations)) CALL FlagError("Vector equations is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector mapping is not associated.",err,error,*999)
#endif    

    vectorEquations=>vectorMapping%vectorEquations

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(vectorEquations)) CALL FlagError("Vector equations is not associated for the vector mapping.",err,error,*999)
#endif    
       
    EXITS("EquationsMappingVector_VectorEquationsGet")
    RETURN
999 NULLIFY(vectorEquations)
998 ERRORSEXITS("EquationsMappingVector_VectorEquationsGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_VectorEquationsGet

  !
  !================================================================================================================================
  !

  !>Gets the vector matrices for vector mapping.
  SUBROUTINE EquationsMappingVector_VectorMatricesGet(vectorMapping,vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the vector equations mapping to get the vector matrices for
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<On exit, a pointer to the vector matrices for the specified vector equations mapping. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingVector_VectorMatricesGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector mapping is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(vectorMapping%vectorEquations)) &
      & CALL FlagError("Vector mapping vector equations is not associated.",err,error,*999)
#endif    
    
    vectorMatrices=>vectorMapping%vectorEquations%vectorMatrices

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is not associated for the vector mapping.",err,error,*999)
#endif    
       
    EXITS("EquationsMappingVector_VectorMatricesGet")
    RETURN
999 NULLIFY(vectorMatrices)
998 ERRORSEXITS("EquationsMappingVector_VectorMatricesGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_VectorMatricesGet

  !
  !================================================================================================================================
  !

  !>Gets the specified dynamic matrix coefficient for a vector equations mapping create values cache.
  SUBROUTINE EquationsMappingVectorCVC_DynamicMatrixCoefficientGet(createValuesCahce,dynamicMatrixIdx, &
    & dynamicMatrixCoefficient,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache !<A pointer to the vector mapping create values cache to get the dynamic matrix coefficient for
    INTEGER(INTG), INTENT(IN) :: dynamicMatrixIdx !<The dynamic matrix index to get the coefficient for
    REAL(DP), INTENT(OUT) :: dynamicMatrixCoefficient !<On exit, the coefficient for the specified dynamic matrix in the create values cache.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingVectorCVC_DynamicMatrixCoefficientGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(createValuesCache)) &
      & CALL FlagError("Vector equations mapping create values cache is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(createValuesCache%dynamicMatrixCoefficients)) &
      & CALL FlagError("Create values cache dynamic matrix coefficients is not allocated.",err,error,*999)
    IF(dynamicMatrixIdx<1.OR.dynamicMatrixIdx>SIZE(createValuesCache%dynamicMatrixCoefficients,1)) THEN
      localError="The specified dynamic matrix index of "//TRIM(NumberToVString(dynamicMatrixIdx,"*",err,error))// &
        & " is invalid. The dynamic matrix index should be be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(createValuesCache%dynamicMatrixCoefficients,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF    
#endif    
    
    dynamicMatrixCoefficient=createValuesCache%dynamicMatrixCoefficients(dynamicMatrixIdx)
    
    EXITS("EquationsMappingVectorCVC_DynamicMatrixCoefficientGet")
    RETURN
999 ERRORS("EquationsMappingVectorCVC_DynamicMatrixCoefficientGet",err,error)
    EXITS("EquationsMappingVectorCVC_DynamicMatrixCoefficientGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingVectorCVC_DynamicMatrixCoefficientGet

  !
  !================================================================================================================================
  !

  !>Gets the specified dynamic variable type for a vector equations mapping create values cache.
  SUBROUTINE EquationsMappingVectorCVC_DynamicVariableTypeGet(createValuesCahce,dynamicVariableType,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache !<A pointer to the vector mapping create values cache to get the dynamic variable type for
    INTEGER(INTG), INTENT(OUT) :: dynamicVariableType !<On exit, the dynamic variable type in the create values cache.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingVectorCVC_DynamicVariableTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(createValuesCache)) &
      & CALL FlagError("Vector equations mapping create values cache is not associated.",err,error,*999)
#endif    
    
    dynamicVariableType=createValuesCache%dynamicVariableType
    
    EXITS("EquationsMappingVectorCVC_DynamicVariableTypeGet")
    RETURN
999 ERRORS("EquationsMappingVectorCVC_DynamicVariableTypeGet",err,error)
    EXITS("EquationsMappingVectorCVC_DynamicVariableTypeGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingVectorCVC_DynamicVariableTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the specified linear matrix coefficient for a vector equations mapping create values cache.
  SUBROUTINE EquationsMappingVectorCVC_LinearMatrixCoefficientGet(createValuesCahce,linearMatrixIdx, &
    & linearMatrixCoefficient,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache !<A pointer to the vector mapping create values cache to get the linear matrix coefficient for
    INTEGER(INTG), INTENT(IN) :: linearMatrixIdx !<The linear matrix index to get the coefficient for
    REAL(DP), INTENT(OUT) :: linearMatrixCoefficient !<On exit, the coefficient for the specified linear matrix in the create values cache.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingVectorCVC_LinearMatrixCoefficientGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(createValuesCache)) &
      & CALL FlagError("Vector equations mapping create values cache is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(createValuesCache%linearMatrixCoefficients)) &
      & CALL FlagError("Create values cache linear matrix coefficients is not allocated.",err,error,*999)
    IF(linearMatrixIdx<1.OR.linearMatrixIdx>SIZE(createValuesCache%linearMatrixCoefficients,1)) THEN
      localError="The specified linear matrix index of "//TRIM(NumberToVString(linearMatrixIdx,"*",err,error))// &
        & " is invalid. The linear matrix index should be be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(createValuesCache%linearMatrixCoefficients,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF    
#endif    
    
    linearMatrixCoefficient=createValuesCache%linearMatrixCoefficients(linearMatrixIdx)
    
    EXITS("EquationsMappingVectorCVC_LinearMatrixCoefficientGet")
    RETURN
999 ERRORS("EquationsMappingVectorCVC_LinearMatrixCoefficientGet",err,error)
    EXITS("EquationsMappingVectorCVC_LinearMatrixCoefficientGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingVectorCVC_LinearMatrixCoefficientGet

  !
  !================================================================================================================================
  !

  !>Gets the specified linear matrix variable type for a vector equations mapping create values cache.
  SUBROUTINE EquationsMappingVectorCVC_LinearMatrixVariableTypeGet(createValuesCahce,linearMatrixIdx, &
    & linearMatrixVariableType,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache !<A pointer to the vector mapping create values cache to get the linear matrix variable type for
    INTEGER(INTG), INTENT(IN) :: linearMatrixIdx !<The linear matrix index to get the variable type for
    INTEGER(INTG), INTENT(OUT) :: linearMatrixVariableType !<On exit, the variable type for the specified linear matrix in the create values cache.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingVectorCVC_LinearMatrixVariableTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(createValuesCache)) &
      & CALL FlagError("Vector equations mapping create values cache is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(createValuesCache%linearMatrixVariableTypes)) &
      & CALL FlagError("Create values cache linear matrix variable types is not allocated.",err,error,*999)
    IF(linearMatrixIdx<1.OR.linearMatrixIdx>SIZE(createValuesCache%linearMatrixVariableTypes,1)) THEN
      localError="The specified linear matrix index of "//TRIM(NumberToVString(linearMatrixIdx,"*",err,error))// &
        & " is invalid. The linear matrix index should be be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(createValuesCache%linearMatrixVariableTypes,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF    
#endif    
    
    linearMatrixVariableType=createValuesCache%linearMatrixVariableTypes(linearMatrixIdx)
    
    EXITS("EquationsMappingVectorCVC_LinearMatrixVariableTypeGet")
    RETURN
999 ERRORS("EquationsMappingVectorCVC_LinearMatrixVariableTypeGet",err,error)
    EXITS("EquationsMappingVectorCVC_LinearMatrixVariableTypeGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingVectorCVC_LinearMatrixVariableTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the numbr of residual variables for a residual in a vector equations mapping create values cache.
  SUBROUTINE EquationsMappingVectorCVC_NumberOfResidualVariablesGet(createValuesCahce,residualIdx, &
    & numberOfResidualVariables,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache !<A pointer to the vector mapping create values cache to get the number of residual variables for
    INTEGER(INTG), INTENT(IN) :: residualIdx !<The residual index to get the number of residual variables for
    INTEGER(INTG), INTENT(OUT) :: numberOfResidualVariables !<On exit, the number residual variables for the specified residual in the create values cache.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingVectorCVC_NumberOfResidualVariablesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(createValuesCache)) &
      & CALL FlagError("Vector equations mapping create values cache is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(createValuesCache%numberOfResdiualVariables)) &
      & CALL FlagError("Create values cache linear matrix variable types is not allocated.",err,error,*999)
    IF(residualIdx<1.OR.residualIdx>SIZE(createValuesCache%numberOfResidualVariables,1)) THEN
      localError="The specified residual index of "//TRIM(NumberToVString(residualIdx,"*",err,error))// &
        & " is invalid. The residual index should be be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(createValuesCache%numberOfResidualVariables,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF    
#endif    
    
    numberOfResidualVariables=createValuesCache%numberOfResidualVariables(residualIdx)
    
    EXITS("EquationsMappingVectorCVC_NumberOfResidualVariablesGet")
    RETURN
999 ERRORS("EquationsMappingVectorCVC_NumberOfResidualVariablesGet",err,error)
    EXITS("EquationsMappingVectorCVC_NumberOfResidualVariablesGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingVectorCVC_NumberOfResidualVariablesGet

  !
  !================================================================================================================================
  !

  !>Gets the specified residual coefficient for a vector equations mapping create values cache.
  SUBROUTINE EquationsMappingVectorCVC_ResidualCoefficientGet(createValuesCahce,residualIdx,residualCoefficient,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache !<A pointer to the vector mapping create values cache to get the residual coefficient for
    INTEGER(INTG), INTENT(IN) :: residualIdx !<The residual index to get the coefficient for
    REAL(DP), INTENT(OUT) :: residualCoefficient !<On exit, the coefficient for the specified residual in the create values cache.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingVectorCVC_ResidualCoefficientGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(createValuesCache)) &
      & CALL FlagError("Vector equations mapping create values cache is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(createValuesCache%residualCoefficients)) &
      & CALL FlagError("Create values cache residual coefficients is not allocated.",err,error,*999)
    IF(residualIdx<1.OR.residualIdx>SIZE(createValuesCache%residualCoefficients,1)) THEN
      localError="The specified residual index of "//TRIM(NumberToVString(residualIdx,"*",err,error))// &
        & " is invalid. The residual index should be be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(createValuesCache%residualCoefficients,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF    
#endif    
    
    residualCoefficient=createValuesCache%residualCoefficients(residualIdx)
    
    EXITS("EquationsMappingVectorCVC_ResidualCoefficientGet")
    RETURN
999 ERRORS("EquationsMappingVectorCVC_ResidualCoefficientGet",err,error)
    EXITS("EquationsMappingVectorCVC_ResidualCoefficientGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingVectorCVC_ResidualCoefficientGet

  !
  !================================================================================================================================
  !

  !>Gets the specified residual variable type for a vector equations mapping create values cache.
  SUBROUTINE EquationsMappingVectorCVC_ResidualVariableTypeGet(createValuesCahce,variableIdx,residualIdx, &
    & residualVariableType,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache !<A pointer to the vector mapping create values cache to get the residual variable type for
    INTEGER(INTG), INTENT(IN) :: variableIdx !<The variable index to get the variable type for
    INTEGER(INTG), INTENT(IN) :: residualIdx !<The residual index to get the variable type for
    INTEGER(INTG), INTENT(OUT) :: residualVariableType !<On exit, the variable type for the specified variable of the residual in the create values cache.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingVectorCVC_ResidualVariableTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(createValuesCache)) &
      & CALL FlagError("Vector equations mapping create values cache is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(createValuesCache%residualVariableTypes)) &
      & CALL FlagError("Create values cache linear matrix variable types is not allocated.",err,error,*999)
    IF(variableIdx<1.OR.variableIdx>SIZE(createValuesCache%residualVariableTypes,1)) THEN
      localError="The specified variable index of "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
        & " is invalid. The variable index should be be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(createValuesCache%residualVariableTypes,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF    
    IF(residualIdx<1.OR.residualIdx>SIZE(createValuesCache%residualVariableTypes,2)) THEN
      localError="The specified residual index of "//TRIM(NumberToVString(residualIdx,"*",err,error))// &
        & " is invalid. The residual index should be be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(createValuesCache%residualVariableTypes,2),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF    
#endif    
    
    residualVariableType=createValuesCache%residualVariableTypes(variableIdx,residualIdx)
    
    EXITS("EquationsMappingVectorCVC_ResidualVariableTypeGet")
    RETURN
999 ERRORS("EquationsMappingVectorCVC_ResidualVariableTypeGet",err,error)
    EXITS("EquationsMappingVectorCVC_ResidualVariableTypeGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingVectorCVC_ResidualVariableTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the specified source coefficient for a vector equations mapping create values cache.
  SUBROUTINE EquationsMappingVectorCVC_SourceCoefficientGet(createValuesCahce,sourceIdx,sourceCoefficient,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache !<A pointer to the vector mapping create values cache to get the source coefficient for
    INTEGER(INTG), INTENT(IN) :: sourceIdx !<The source index to get the coefficient for
    REAL(DP), INTENT(OUT) :: sourceCoefficient !<On exit, the coefficient for the specified source in the create values cache.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingVectorCVC_SourceCoefficientGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(createValuesCache)) &
      & CALL FlagError("Vector equations mapping create values cache is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(createValuesCache%sourceCoefficients)) &
      & CALL FlagError("Create values cache source coefficients is not allocated.",err,error,*999)
    IF(sourceIdx<1.OR.sourceIdx>SIZE(createValuesCache%sourceCoefficients,1)) THEN
      localError="The specified source index of "//TRIM(NumberToVString(sourceIdx,"*",err,error))// &
        & " is invalid. The source index should be be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(createValuesCache%sourceCoefficients,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF    
#endif    
    
    sourceCoefficient=createValuesCache%sourceCoefficients(sourceIdx)
    
    EXITS("EquationsMappingVectorCVC_SourceCoefficientGet")
    RETURN
999 ERRORS("EquationsMappingVectorCVC_SourceCoefficientGet",err,error)
    EXITS("EquationsMappingVectorCVC_SourceCoefficientGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingVectorCVC_SourceCoefficientGet

  !
  !================================================================================================================================
  !

  !>Gets the specified source variable type for a vector equations mapping create values cache.
  SUBROUTINE EquationsMappingVectorCVC_SourceVariableTypeGet(createValuesCahce,sourceIdx,sourceVariableType,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache !<A pointer to the vector mapping create values cache to get the source variable type for
    INTEGER(INTG), INTENT(IN) :: sourceIdx !<The source index to get the variable type for
    INTEGER(INTG), INTENT(OUT) :: sourceVariableType !<On exit, the variable type for the specified variable of the source in the create values cache.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingVectorCVC_SourceVariableTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(createValuesCache)) &
      & CALL FlagError("Vector equations mapping create values cache is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(createValuesCache%sourceVariableTypes)) &
      & CALL FlagError("Create values cache linear matrix variable types is not allocated.",err,error,*999)
    IF(sourceIdx<1.OR.sourceIdx>SIZE(createValuesCache%sourceVariableTypes,1)) THEN
      localError="The specified source index of "//TRIM(NumberToVString(sourceIdx,"*",err,error))// &
        & " is invalid. The source index should be be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(createValuesCache%sourceVariableTypes,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF    
#endif    
    
    sourceVariableType=createValuesCache%sourceVariableTypes(sourceIdx)
    
    EXITS("EquationsMappingVectorCVC_SourceVariableTypeGet")
    RETURN
999 ERRORS("EquationsMappingVectorCVC_SourceVariableTypeGet",err,error)
    EXITS("EquationsMappingVectorCVC_SourceVariableTypeGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingVectorCVC_SourceVariableTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the equations matrix for a equations matrix to variable map.
  SUBROUTINE EquationsMappingVectorEMToVMap_EquationsMatrixGet(equationsMatrixToVarMap,equationsMatrix,err,error,*)

    !Argument variables
    TYPE(EquationsMatrixToVarMapType), POINTER :: equationsMatrixToVarMap !<A pointer to the equations matrix to variable map to get the equations matrix for
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix !<On exit, a pointer to the equations matrix. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingVectorEMToVMap_EquationsMatrixGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(equationsMatrix)) CALL FlagError("Equations matrix is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsMatrixToVarMap)) &
      & CALL FlagError("Equations matrix to variable map is not associated.",err,error,*999)
#endif    
    
    equationsMatrix=>equationsMatrixToVarMap%equationsMatrix

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(equationsMatrix)) &
      & CALL FlagError("The equations matrix is not associated for the equations matrix to variable map.",err,error,*999)
#endif    
    
    EXITS("EquationsMappingVectorEMToVMap_EquationsMatrixGet")
    RETURN
999 NULLIFY(equationsMatrix)
998 ERRORS("EquationsMappingVectorEMToVMap_EquationsMatrixGet",err,error)
    EXITS("EquationsMappingVectorEMToVMap_EquationsMatrixGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingVectorEMToVMap_EquationsMatrixGet

  !
  !================================================================================================================================
  !

  !>Gets the field variable for a equations matrix to variable map.
  SUBROUTINE EquationsMappingVectorEMToVMap_VariableGet(equationsMatrixToVarMap,variable,err,error,*)

    !Argument variables
    TYPE(EquationsMatrixToVarMapType), POINTER :: equationsMatrixToVarMap !<A pointer to the equations matrix to variable map to get the field variable for
    TYPE(FieldVariableType), POINTER :: variable !<On exit, a pointer to the field variable. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingVectorEMToVMap_VariableGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(variable)) CALL FlagError("Variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsMatrixToVarMap)) &
      & CALL FlagError("Equations matrix to variable map is not associated.",err,error,*999)
#endif    
    
    variable=>equationsMatrixToVarMap%variable

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(variable)) &
      & CALL FlagError("The variable is not associated for the equations matrix to variable map.",err,error,*999)
#endif    
    
    EXITS("EquationsMappingVectorEMToVMap_VariableGet")
    RETURN
999 NULLIFY(variable)
998 ERRORS("EquationsMappingVectorEMToVMap_VariableGet",err,error)
    EXITS("EquationsMappingVectorEMToVMap_VariableGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingVectorEMToVMap_VariableGet

  !
  !================================================================================================================================
  !

  !>Gets the Jacobian matrix for a Jacobian matrix to variable map.
  SUBROUTINE EquationsMappingVectorJMToVMap_JacobianMatrixGet(jacobianMatrixToVarMap,jacobianMatrix,err,error,*)

    !Argument variables
    TYPE(JacobianMatrixToVarMapType), POINTER :: jacobianMatrixToVarMap !<A pointer to the Jacobian matrix to variable map to get the Jacobian matrix for
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix !<On exit, a pointer to the Jacobian matrix. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingVectorJMToVMap_JacobianMatrixGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(jacobianMatrix)) CALL FlagError("Jacobian matrix is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(jacobianMatrixToVarMap)) &
      & CALL FlagError("Jacobian matrix to variable map is not associated.",err,error,*999)
#endif    
    
    jacobianMatrix=>jacobianMatrixToVarMap%jacobianMatrix

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(jacobianMatrix)) &
      & CALL FlagError("The Jacobian matrix is not associated for the Jacobian matrix to variable map.",err,error,*999)
#endif    
    
    EXITS("EquationsMappingVectorJMToVMap_JacobianMatrixGet")
    RETURN
999 NULLIFY(jacobianMatrix)
998 ERRORS("EquationsMappingVectorJMToVMap_JacobianMatrixGet",err,error)
    EXITS("EquationsMappingVectorJMToVMap_JacobianMatrixGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingVectorJMToVMap_JacobianMatrixGet

  !
  !================================================================================================================================
  !

  !>Gets the field variable for a jacobian matrix to variable map.
  SUBROUTINE EquationsMappingVectorJMToVMap_VariableGet(jacobianMatrixToVarMap,variable,err,error,*)

    !Argument variables
    TYPE(JacobianMatrixToVarMapType), POINTER :: jacobianMatrixToVarMap !<A pointer to the jacobian matrix to variable map to get the field variable for
    TYPE(FieldVariableType), POINTER :: variable !<On exit, a pointer to the field variable. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingVectorJMToVMap_VariableGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(variable)) CALL FlagError("Variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(jacobianMatrixToVarMap)) &
      & CALL FlagError("Jacobian matrix to variable map is not associated.",err,error,*999)
#endif    
    
    variable=>jacobianMatrixToVarMap%variable

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(variable)) &
      & CALL FlagError("The variable is not associated for the jacobian matrix to variable map.",err,error,*999)
#endif    
    
    EXITS("EquationsMappingVectorJMToVMap_VariableGet")
    RETURN
999 NULLIFY(variable)
998 ERRORS("EquationsMappingVectorJMToVMap_VariableGet",err,error)
    EXITS("EquationsMappingVectorJMToVMap_VariableGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingVectorJMToVMap_VariableGet

  !
  !================================================================================================================================
  !

  !>Gets the field variable for a variable to equations matrices map.
  SUBROUTINE EquationsMappingVectorVToEMSMap_VariableGet(varToEquationsMatricesMap,variable,err,error,*)

    !Argument variables
    TYPE(VarToEquationsMatricesMapType), POINTER :: varToEquationsMatricesMap !<A pointer to the variable to equations matrices map to get the field variable for
    TYPE(FieldVariableType), POINTER :: variable !<On exit, a pointer to the field variable. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingVectorVToEMSMap_VariableGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(variable)) CALL FlagError("Variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(varToEquationsMatricesMap)) &
      & CALL FlagError("Variable to equations matrices map is not associated.",err,error,*999)
#endif    
    
    variable=>varToEquationsMatricesMap%variable

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(variable)) &
      & CALL FlagError("The variable is not associated for the variable to equations matrices map.",err,error,*999)
#endif    
    
    EXITS("EquationsMappingVectorVToEMSMap_VariableGet")
    RETURN
999 NULLIFY(variable)
998 ERRORS("EquationsMappingVectorVToEMSMap_VariableGet",err,error)
    EXITS("EquationsMappingVectorVToEMSMap_VariableGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingVectorVToEMSMap_VariableGet

  !
  !================================================================================================================================
  !

  !>Gets the field variable for a variable to Jacobian matrix map.
  SUBROUTINE EquationsMappingVectorVToJMMap_VariableGet(varToJacobianMatrixMap,variable,err,error,*)

    !Argument variables
    TYPE(VarToJacobianMatrixMapType), POINTER :: varToJacobianMatrixMap !<A pointer to the variable to Jacobian matrix map to get the field variable for
    TYPE(FieldVariableType), POINTER :: variable !<On exit, a pointer to the field variable. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingVectorVToJMMap_VariableGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(variable)) CALL FlagError("Variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(varToJacobianMatrixMap)) &
      & CALL FlagError("Variable to Jacobian matrix map is not associated.",err,error,*999)
#endif    
    
    variable=>varToJacobianMatrixMap%variable

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(variable)) &
      & CALL FlagError("The variable is not associated for the variable to Jacobian matrix map.",err,error,*999)
#endif    
    
    EXITS("EquationsMappingVectorVToJMMap_VariableGet")
    RETURN
999 NULLIFY(varToJacobianMatrixsMap)
998 ERRORS("EquationsMappingVectorVToJMMap_VariableGet",err,error)
    EXITS("EquationsMappingVectorVToJMMap_VariableGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingVectorVToJMMap_VariableGet

  !
  !================================================================================================================================
  !
  
END MODULE EquationsMappingAccessRoutines
