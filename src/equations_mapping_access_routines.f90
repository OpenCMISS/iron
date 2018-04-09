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

  PUBLIC EquationsMappingLinear_LinearVariableGet

  PUBLIC EquationsMappingNonlinear_ResidualVariableGet
  
  PUBLIC EquationsMappingRHS_RHSVariableGet
  
  PUBLIC EquationsMappingScalar_CreateValuesCacheGet
  
  PUBLIC EquationsMappingScalar_ScalarEquationsGet

  PUBLIC EquationsMappingSource_SourceVariableGet

  PUBLIC EquationsMappingVector_CreateValuesCacheGet

  PUBLIC EquationsMappingVector_DynamicMappingGet
  
  PUBLIC EquationsMappingVector_LHSMappingGet
  
  PUBLIC EquationsMappingVector_LinearMappingGet
  
  PUBLIC EquationsMappingVector_NonlinearMappingGet
  
  PUBLIC EquationsMappingVector_RHSMappingGet
  
  PUBLIC EquationsMappingVector_SourceMappingGet

  PUBLIC EquationsMappingVector_VectorEquationsGet
  
  PUBLIC EquationsMappingVector_VectorMatricesGet
 
CONTAINS

  !
  !================================================================================================================================
  !

  !>Gets the specified dynamic variable for a dynamic mapping.
  SUBROUTINE EquationsMappingDynamic_DynamicVariableGet(dynamicMapping,dynamicVariable,err,error,*)

    !Argument variables
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping !<A pointer to the dynamic mapping to get the dynamic variable for
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: dynamicVariable !<On exit, a pointer to the requested dynamic variable. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingDynamic_DynamicVariableGet",err,error,*998)

    IF(ASSOCIATED(dynamicVariable)) CALL FlagError("Dynamic variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(dynamicMapping)) CALL FlagError("Dynamic mapping is not associated.",err,error,*999)
    
    dynamicVariable=>dynamicMapping%dynamicVariable  
    IF(.NOT.ASSOCIATED(dynamicVariable)) CALL FlagError("The dynamic variable is not associated.",err,error,*999)
    
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

  !>Gets the specified linear variable for a linear mapping.
  SUBROUTINE EquationsMappingLinear_LinearVariableGet(linearMapping,linearMatrixIdx,linearVariable,err,error,*)

    !Argument variables
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping !<A pointer to the linear mapping to get the linear variable for
    INTEGER(INTG), INTENT(IN) :: linearMatrixIdx !<The index of the linear matrix to get the linear variable for. 
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: linearVariable !<On exit, a pointer to the requested linear variable. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsMappingLinear_LinearVariableGet",err,error,*998)

    IF(ASSOCIATED(linearVariable)) CALL FlagError("Linear variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(linearMapping)) CALL FlagError("Linear mapping is not associated.",err,error,*999)
    IF(linearMatrixIdx<1.OR.linearMatrixIdx>linearMapping%numberOfLinearMatrices) THEN
      localError="The specified linear matrix index of "//TRIM(NumberToVString(linearMatrixIdx,"*",err,error))// &
        & " is invalid for a linear mapping which has "// &
        & TRIM(NumberToVString(linearMapping%numberOfLinearMatrices,"*",err,error))//" linear matrices."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(linearMapping%equationsMatrixToVarMaps)) &
      & CALL FlagError("Linear mapping equations matrix to variable maps is not allocated.",err,error,*999)
    
    linearVariable=>linearMapping%equationsMatrixToVarMaps(linearMatrixIdx)%variable
    IF(.NOT.ASSOCIATED(linearVariable)) THEN
      localError="The linear variable is not associated for linear matrix index "// &
        & TRIM(NumberToVString(linearMatrixIdx,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
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

  !>Gets the specified residual variable for a nonlinear mapping.
  SUBROUTINE EquationsMappingNonlinear_ResidualVariableGet(nonlinearMapping,residualIdx,variableIdx,fieldVariable,err,error,*)

    !Argument variables
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping !<A pointer to the nonlinear mapping to get the residual variable for
    INTEGER(INTG), INTENT(IN) :: residualIdx !<The index of the residual to get the field variable for. Currently will just be 1.
    INTEGER(INTG), INTENT(IN) :: variableIdx !<The index of the variable in the residual to get the field variable for. 
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable !<On exit, a pointer to the requested field variable. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsMappingNonlinear_ResidualVariableGet",err,error,*998)

    IF(ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(nonlinearMapping)) CALL FlagError("Nonlinear mapping is not associated.",err,error,*999)
    IF(residualIdx<1.OR.residualIdx>nonlinearMapping%numberOfResiduals) THEN
      localError="The specified residual index of "//TRIM(NumberToVString(residualIdx,"*",err,error))// &
        & " is invalid. The residual index must be >= 1 and <= "// &
        & TRIM(NumberToVString(nonlinearMapping%numberOfResiduals,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(variableIdx<1.OR.variableIdx>nonlinearMapping%numberOfResidualVariables) THEN
      localError="The specified residual variable index of "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
        & " is invalid for residual index "//TRIM(NumberToVString(residualIdx,"*",err,error))//" which has "// &
        & TRIM(NumberToVString(nonlinearMapping%numberOfResidualVariables,"*",err,error))//" variables."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(nonlinearMapping%residualVariables)) &
      & CALL FlagError("Nonlinear mapping residual variables is not allocated.",err,error,*999)
    
    fieldVariable=>nonlinearMapping%residualVariables(variableIdx)%ptr
    IF(.NOT.ASSOCIATED(fieldVariable)) THEN
      localError="The field variable is not associated for variable index "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
        & " of residual index "//TRIM(NumberToVString(residualIdx,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("EquationsMappingNonlinear_ResidualVariableGet")
    RETURN
999 NULLIFY(fieldVariable)
998 ERRORS("EquationsMappingNonlinear_ResidualVariableGet",err,error)
    EXITS("EquationsMappingNonlinear_ResidualVariableGet")
    RETURN 1
    
  END SUBROUTINE EquationsMappingNonlinear_ResidualVariableGet

  !
  !================================================================================================================================
  !

  !>Gets the specified RHS variable for a rhs mapping.
  SUBROUTINE EquationsMappingRHS_RHSVariableGet(rhsMapping,rhsVariable,err,error,*)

    !Argument variables
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping !<A pointer to the RHS mapping to get the RHS variable for
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: rhsVariable !<On exit, a pointer to the requested RHS field variable. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsMappingRHS_RHSVariableGet",err,error,*998)

    IF(ASSOCIATED(rhsVariable)) CALL FlagError("RHS Field variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(rhsMapping)) CALL FlagError("RHS mapping is not associated.",err,error,*999)
    
    rhsVariable=>rhsMapping%rhsVariable
    IF(.NOT.ASSOCIATED(rhsVariable)) CALL FlagError("The RHS field variable is not associated in the RHS mapping.",err,error,*999)
    
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

  !>Gets the create values cache for an scalar equations mapping.
  SUBROUTINE EquationsMappingScalar_CreateValuesCacheGet(scalarMapping,createValuesCache,err,error,*)

    !Argument variables
    TYPE(EquationsMappingScalarType), POINTER :: scalarMapping !<A pointer to the equations scalar mapping to get the create values cache for
    TYPE(EquationsMappingScalarCreateValuesCacheType), POINTER :: createValuesCache !<On exit, a pointer to the create values cache in the specified scalar equations mapping. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingScalar_CreateValuesCacheGet",err,error,*998)

    IF(ASSOCIATED(createValuesCache)) CALL FlagError("Create values cache is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(scalarMapping)) CALL FlagError("Scalar mapping is not associated.",err,error,*999)

    createValuesCache=>scalarMapping%createValuesCache
    IF(.NOT.ASSOCIATED(createValuesCache)) &
      & CALL FlagError("Create values cache is not associated for the scalar mapping.",err,error,*999)
       
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

    IF(ASSOCIATED(scalarEquations)) CALL FlagError("Scalar equations is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(scalarMapping)) CALL FlagError("Scalar mapping is not associated.",err,error,*999)

    scalarEquations=>scalarMapping%scalarEquations
    IF(.NOT.ASSOCIATED(scalarEquations)) CALL FlagError("Scalar equations is not associated for the scalar mapping.",err,error,*999)
       
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
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: sourceVariable !<On exit, a pointer to the requested source variable. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingSoure_SourceVariableGet",err,error,*998)

    IF(ASSOCIATED(sourceVariable)) CALL FlagError("Source variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(sourceMapping)) CALL FlagError("Soure mapping is not associated.",err,error,*999)
    
    sourceVariable=>sourceMapping%sourceVariable  
    IF(.NOT.ASSOCIATED(sourceVariable)) CALL FlagError("The source variable is not associated.",err,error,*999)
    
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

  !>Gets the create values cache for an vector equations mapping.
  SUBROUTINE EquationsMappingVector_CreateValuesCacheGet(vectorMapping,createValuesCache,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations vector mapping to get the create values cache for
    TYPE(EquationsMappingVectorCreateValuesCacheType), POINTER :: createValuesCache !<On exit, a pointer to the create values cache in the specified vector equations mapping. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingVector_CreateValuesCacheGet",err,error,*998)

    IF(ASSOCIATED(createValuesCache)) CALL FlagError("Create values cache is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector mapping is not associated.",err,error,*999)

    createValuesCache=>vectorMapping%createValuesCache
    IF(.NOT.ASSOCIATED(createValuesCache)) &
      & CALL FlagError("Create values cache is not associated for the vector mapping.",err,error,*999)
       
    EXITS("EquationsMappingVector_CreateValuesCacheGet")
    RETURN
999 NULLIFY(createValuesCache)
998 ERRORSEXITS("EquationsMappingVector_CreateValuesCacheGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_CreateValuesCacheGet

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

    IF(ASSOCIATED(dynamicMapping)) CALL FlagError("Dynamic mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector mapping is not associated.",err,error,*999)

    dynamicMapping=>vectorMapping%dynamicMapping
    IF(.NOT.ASSOCIATED(dynamicMapping)) CALL FlagError("Dynamic mapping is not associated for the vector mapping.",err,error,*999)
       
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

    IF(ASSOCIATED(lhsMapping)) CALL FlagError("LHS mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector mapping is not associated.",err,error,*999)

    lhsMapping=>vectorMapping%lhsMapping
    IF(.NOT.ASSOCIATED(lhsMapping)) CALL FlagError("LHS mapping is not associated for the vector mapping.",err,error,*999)
       
    EXITS("EquationsMappingVector_LHSMappingGet")
    RETURN
999 NULLIFY(lhsMapping)
998 ERRORSEXITS("EquationsMappingVector_LHSMappingGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_LHSMappingGet

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

    IF(ASSOCIATED(linearMapping)) CALL FlagError("Linear mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector mapping is not associated.",err,error,*999)

    linearMapping=>vectorMapping%linearMapping
    IF(.NOT.ASSOCIATED(linearMapping)) CALL FlagError("Linear mapping is not associated for the vector mapping.",err,error,*999)
       
    EXITS("EquationsMappingVector_LinearMappingGet")
    RETURN
999 NULLIFY(linearMapping)
998 ERRORSEXITS("EquationsMappingVector_LinearMappingGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_LinearMappingGet

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

    IF(ASSOCIATED(nonlinearMapping)) CALL FlagError("Nonlinear mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector mapping is not associated.",err,error,*999)

    nonlinearMapping=>vectorMapping%nonlinearMapping
    IF(.NOT.ASSOCIATED(nonlinearMapping))  &
      & CALL FlagError("Nonlinear mapping is not associated for the vector mapping.",err,error,*999)
       
    EXITS("EquationsMappingVector_NonlinearMappingGet")
    RETURN
999 NULLIFY(nonlinearMapping)
998 ERRORSEXITS("EquationsMappingVector_NonlinearMappingGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_NonlinearMappingGet

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

    IF(ASSOCIATED(rhsMapping)) CALL FlagError("RHS mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector mapping is not associated.",err,error,*999)

    rhsMapping=>vectorMapping%rhsMapping
    IF(.NOT.ASSOCIATED(rhsMapping)) CALL FlagError("RHS mapping is not associated for the vector mapping.",err,error,*999)
       
    EXITS("EquationsMappingVector_RHSMappingGet")
    RETURN
999 NULLIFY(rhsMapping)
998 ERRORSEXITS("EquationsMappingVector_RHSMappingGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_RHSMappingGet

  !
  !================================================================================================================================
  !

  !>Gets the source vector mapping for an vector mapping.
  SUBROUTINE EquationsMappingVector_SourceMappingGet(vectorMapping,sourceMapping,err,error,*)

    !Argument variables
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<A pointer to the equations vector mapping to get the source mapping for
    TYPE(EquationsMappingSourceType), POINTER :: sourceMapping !<On exit, a pointer to the source Mapping in the specified vector equations mapping. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsMappingVector_SourceMappingGet",err,error,*998)

    IF(ASSOCIATED(sourceMapping)) CALL FlagError("Source mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector mapping is not associated.",err,error,*999)

    sourceMapping=>vectorMapping%sourceMapping
    IF(.NOT.ASSOCIATED(sourceMapping)) CALL FlagError("Source mapping is not associated for the vector mapping.",err,error,*999)
       
    EXITS("EquationsMappingVector_SourceMappingGet")
    RETURN
999 NULLIFY(sourceMapping)
998 ERRORSEXITS("EquationsMappingVector_SourceMappingGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_SourceMappingGet

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

    IF(ASSOCIATED(vectorEquations)) CALL FlagError("Vector equations is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector mapping is not associated.",err,error,*999)

    vectorEquations=>vectorMapping%vectorEquations
    IF(.NOT.ASSOCIATED(vectorEquations)) CALL FlagError("Vector equations is not associated for the vector mapping.",err,error,*999)
       
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

    IF(ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector mapping is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(vectorMapping%vectorEquations)) &
      & CALL FlagError("Vector mapping vector equations is not associated.",err,error,*999)
    
    vectorMatrices=>vectorMapping%vectorEquations%vectorMatrices
    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is not associated for the vector mapping.",err,error,*999)
       
    EXITS("EquationsMappingVector_VectorMatricesGet")
    RETURN
999 NULLIFY(vectorMatrices)
998 ERRORSEXITS("EquationsMappingVector_VectorMatricesGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMappingVector_VectorMatricesGet

  !
  !================================================================================================================================
  !

END MODULE EquationsMappingAccessRoutines
