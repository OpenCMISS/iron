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
  
  USE BaseRoutines
  USE FieldAccessRoutines
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
  
  PUBLIC Equations_InterpolationGet
  
  PUBLIC Equations_ScalarEquationsGet
  
  PUBLIC Equations_VectorEquationsGet
  
  PUBLIC EquationsInterpolation_DependentParametersGet
  
  PUBLIC EquationsInterpolation_DependentPointGet
  
  PUBLIC EquationsInterpolation_DependentPointMetricsGet
  
  PUBLIC EquationsInterpolation_FibreParametersGet
  
  PUBLIC EquationsInterpolation_FibrePointGet
  
  PUBLIC EquationsInterpolation_FibrePointMetricsGet
  
  PUBLIC EquationsInterpolation_GeometricParametersGet
  
  PUBLIC EquationsInterpolation_GeometricPointGet
  
  PUBLIC EquationsInterpolation_GeometricPointMetricsGet
  
  PUBLIC EquationsInterpolation_IndependentParametersGet
  
  PUBLIC EquationsInterpolation_IndependentPointGet
  
  PUBLIC EquationsInterpolation_IndependentPointMetricsGet
  
  PUBLIC EquationsInterpolation_MaterialsParametersGet
  
  PUBLIC EquationsInterpolation_MaterialsPointGet
  
  PUBLIC EquationsInterpolation_PreviousDependentParametersGet
  
  PUBLIC EquationsInterpolation_PreviousDependentPointGet
  
  PUBLIC EquationsInterpolation_PreviousDependentPointMetricsGet
  
  PUBLIC EquationsInterpolation_SourceParametersGet
  
  PUBLIC EquationsInterpolation_SourcePointGet
  
  PUBLIC EquationsScalar_EquationsGet

  PUBLIC EquationsScalar_ScalarMappingGet

  PUBLIC EquationsScalar_ScalarMatricesGet

  PUBLIC EquationsVector_EquationsGet
  
  PUBLIC EquationsVector_VectorMappingGet
  
  PUBLIC EquationsVector_VectorMatricesGet
  
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

  !>Gets the interpolation for an equations.
  SUBROUTINE Equations_InterpolationGet(equations,interpolation,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER :: equations !<A pointer to the equations to get the interpolation for
    TYPE(EquationsInterpolationType), POINTER :: interpolation !<On exit, a pointer to the equations interpolation for the specified equations. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Equations_InterpolationGet",err,error,*998)

    IF(ASSOCIATED(interpolation)) CALL FlagError("Interpolation is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)

    interpolation=>equations%interpolation
    IF(.NOT.ASSOCIATED(interpolation)) CALL FlagError("Interpolation is not associated for the equations.",err,error,*999)
       
    EXITS("Equations_InterpolationGet")
    RETURN
999 NULLIFY(interpolation)
998 ERRORSEXITS("Equations_InterpolationGet",err,error)
    RETURN 1
    
  END SUBROUTINE Equations_InterpolationGet

  !
  !================================================================================================================================
  !

  !>Gets the scalar equations for equations.
  SUBROUTINE Equations_ScalarEquationsGet(equations,scalarEquations,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER :: equations !<A pointer to the equations to get the scalar equations for
    TYPE(EquationsScalarType), POINTER :: scalarEquations !<On exit, a pointer to the scalar equations for the specified equations. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Equations_ScalarEquationsGet",err,error,*998)

    IF(ASSOCIATED(scalarEquations)) CALL FlagError("Scalar equations is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)

    scalarEquations=>equations%scalarEquations
    IF(.NOT.ASSOCIATED(scalarEquations)) CALL FlagError("Scalar equations is not associated for the equations.",err,error,*999)
       
    EXITS("Equations_ScalarEquationsGet")
    RETURN
999 NULLIFY(scalarEquations)
998 ERRORSEXITS("Equations_ScalarEquationsGet",err,error)
    RETURN 1
    
  END SUBROUTINE Equations_ScalarEquationsGet

  !
  !================================================================================================================================
  !

  !>Gets the vector equations for equations.
  SUBROUTINE Equations_VectorEquationsGet(equations,vectorEquations,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER :: equations !<A pointer to the equations to get the vector equations for
    TYPE(EquationsVectorType), POINTER :: vectorEquations !<On exit, a pointer to the vector equations for the specified equations. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Equations_VectorEquationsGet",err,error,*998)

    IF(ASSOCIATED(vectorEquations)) CALL FlagError("Vector equations is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)

    vectorEquations=>equations%vectorEquations
    IF(.NOT.ASSOCIATED(vectorEquations)) CALL FlagError("Vector equations is not associated for the equations.",err,error,*999)
       
    EXITS("Equations_VectorEquationsGet")
    RETURN
999 NULLIFY(vectorEquations)
998 ERRORSEXITS("Equations_VectorEquationsGet",err,error)
    RETURN 1
    
  END SUBROUTINE Equations_VectorEquationsGet

  !
  !================================================================================================================================
  !

  !>Gets the dependent interpolation parameters for an equations interpolation.
  SUBROUTINE EquationsInterpolation_DependentParametersGet(equationsInterpolation,variableType,dependentParameters,err,error,*)

    !Argument variables
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation !<A pointer to the equations interpolation to get the dependent parameters for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to get the dependent parameters for. \see FIELD_ROUTINES_VariableTypes
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: dependentParameters !<On exit, a pointer to the dependent interpolation parameters for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsInterpolation_DependentParametersGet",err,error,*998)

    IF(ASSOCIATED(dependentParameters)) CALL FlagError("Dependent parameters is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
    IF(variableType<1.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    dependentParameters=>equationsInterpolation%dependentInterpParameters(variableType)%ptr
    IF(.NOT.ASSOCIATED(dependentParameters)) THEN
      localError="Equations interpolation dependent parameters is not associated for field variable type "// &
        & TRIM(NumberToVString(variableType,"*",err,error))//" of the equations interpolation."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("EquationsInterpolation_DependentParametersGet")
    RETURN
999 NULLIFY(dependentParameters)
998 ERRORS("EquationsInterpolation_DependentParametersGet",err,error)
    EXITS("EquationsInterpolation_DependentParametersGet")
    RETURN 1
    
  END SUBROUTINE EquationsInterpolation_DependentParametersGet

  !
  !================================================================================================================================
  !

  !>Gets the dependent interpolated point for an equations interpolation.
  SUBROUTINE EquationsInterpolation_DependentPointGet(equationsInterpolation,variableType,dependentPoint,err,error,*)

    !Argument variables
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation !<A pointer to the equations interpolation to get the dependent point for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to get the dependent point for. \see FIELD_ROUTINES_VariableTypes
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: dependentPoint !<On exit, a pointer to the dependent interpolated point for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsInterpolation_DependentPointGet",err,error,*998)

    IF(ASSOCIATED(dependentPoint)) CALL FlagError("Dependent point is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
    IF(variableType<1.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    dependentPoint=>equationsInterpolation%dependentInterpPoint(variableType)%ptr
    IF(.NOT.ASSOCIATED(dependentPoint)) THEN
      localError="Equations interpolated dependent point is not associated for field variable type "// &
        & TRIM(NumberToVString(variableType,"*",err,error))//" of the equations interpolation."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("EquationsInterpolation_DependentPointGet")
    RETURN
999 NULLIFY(dependentPoint)
998 ERRORS("EquationsInterpolation_DependentPointGet",err,error)
    EXITS("EquationsInterpolation_DependentPointGet")
    RETURN 1
    
  END SUBROUTINE EquationsInterpolation_DependentPointGet

  !
  !================================================================================================================================
  !

  !>Gets the dependent interpolated point metrics for an equations interpolation.
  SUBROUTINE EquationsInterpolation_DependentPointMetricsGet(equationsInterpolation,variableType,dependentPointMetrics,err,error,*)

    !Argument variables
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation !<A pointer to the equations interpolation to get the dependent point metrics for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to get the dependent point metrics for. \see FIELD_ROUTINES_VariableTypes
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: dependentPointMetrics !<On exit, a pointer to the dependent interpolated point metrics for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsInterpolation_DependentPointMetricsGet",err,error,*998)

    IF(ASSOCIATED(dependentPointMetrics)) CALL FlagError("Dependent point metrics is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
    IF(variableType<1.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    dependentPointMetrics=>equationsInterpolation%dependentInterpPointMetrics(variableType)%ptr
    IF(.NOT.ASSOCIATED(dependentPointMetrics)) THEN
      localError="Equations interpolated dependent point metrics is not associated for field variable type "// &
        & TRIM(NumberToVString(variableType,"*",err,error))//" of the equations interpolation."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("EquationsInterpolation_DependentPointMetricsGet")
    RETURN
999 NULLIFY(dependentPointMetrics)
998 ERRORS("EquationsInterpolation_DependentPointMetricsGet",err,error)
    EXITS("EquationsInterpolation_DependentPointMetricsGet")
    RETURN 1
    
  END SUBROUTINE EquationsInterpolation_DependentPointMetricsGet

  !
  !================================================================================================================================
  !

  !>Gets the fibre interpolation parameters for an equations interpolation.
  SUBROUTINE EquationsInterpolation_FibreParametersGet(equationsInterpolation,variableType,fibreParameters,err,error,*)

    !Argument variables
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation !<A pointer to the equations interpolation to get the fibre parameters for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to get the fibre parameters for. \see FIELD_ROUTINES_VariableTypes
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: fibreParameters !<On exit, a pointer to the fibre interpolation parameters for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsInterpolation_FibreParametersGet",err,error,*998)

    IF(ASSOCIATED(fibreParameters)) CALL FlagError("Fibre parameters is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
    IF(variableType<1.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    fibreParameters=>equationsInterpolation%fibreInterpParameters(variableType)%ptr
    IF(.NOT.ASSOCIATED(fibreParameters)) THEN
      localError="Equations interpolation fibre parameters is not associated for field variable type "// &
        & TRIM(NumberToVString(variableType,"*",err,error))//" of the equations interpolation."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("EquationsInterpolation_FibreParametersGet")
    RETURN
999 NULLIFY(fibreParameters)
998 ERRORS("EquationsInterpolation_FibreParametersGet",err,error)
    EXITS("EquationsInterpolation_FibreParametersGet")
    RETURN 1
    
  END SUBROUTINE EquationsInterpolation_FibreParametersGet

  !
  !================================================================================================================================
  !

  !>Gets the fibre interpolated point for an equations interpolation.
  SUBROUTINE EquationsInterpolation_FibrePointGet(equationsInterpolation,variableType,fibrePoint,err,error,*)

    !Argument variables
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation !<A pointer to the equations interpolation to get the fibre point for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to get the fibre point for. \see FIELD_ROUTINES_VariableTypes
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: fibrePoint !<On exit, a pointer to the fibre interpolated point for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsInterpolation_FibrePointGet",err,error,*998)

    IF(ASSOCIATED(fibrePoint)) CALL FlagError("Fibre point is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
    IF(variableType<1.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    fibrePoint=>equationsInterpolation%fibreInterpPoint(variableType)%ptr
    IF(.NOT.ASSOCIATED(fibrePoint)) THEN
      localError="Equations interpolated fibre point is not associated for field variable type "// &
        & TRIM(NumberToVString(variableType,"*",err,error))//" of the equations interpolation."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("EquationsInterpolation_FibrePointGet")
    RETURN
999 NULLIFY(fibrePoint)
998 ERRORS("EquationsInterpolation_FibrePointGet",err,error)
    EXITS("EquationsInterpolation_FibrePointGet")
    RETURN 1
    
  END SUBROUTINE EquationsInterpolation_FibrePointGet

  !
  !================================================================================================================================
  !

  !>Gets the fibre interpolated point metrics for an equations interpolation.
  SUBROUTINE EquationsInterpolation_FibrePointMetricsGet(equationsInterpolation,variableType,fibrePointMetrics,err,error,*)

    !Argument variables
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation !<A pointer to the equations interpolation to get the fibre point metrics for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to get the fibre point metrics for. \see FIELD_ROUTINES_VariableTypes
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: fibrePointMetrics !<On exit, a pointer to the fibre interpolated point metrics for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsInterpolation_FibrePointMetricsGet",err,error,*998)

    IF(ASSOCIATED(fibrePointMetrics)) CALL FlagError("Fibre point metrics is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
    IF(variableType<1.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    fibrePointMetrics=>equationsInterpolation%fibreInterpPointMetrics(variableType)%ptr
    IF(.NOT.ASSOCIATED(fibrePointMetrics)) THEN
      localError="Equations interpolated fibre point metrics is not associated for field variable type "// &
        & TRIM(NumberToVString(variableType,"*",err,error))//" of the equations interpolation."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("EquationsInterpolation_FibrePointMetricsGet")
    RETURN
999 NULLIFY(fibrePointMetrics)
998 ERRORS("EquationsInterpolation_FibrePointMetricsGet",err,error)
    EXITS("EquationsInterpolation_FibrePointMetricsGet")
    RETURN 1
    
  END SUBROUTINE EquationsInterpolation_FibrePointMetricsGet

  !
  !================================================================================================================================
  !

  !>Gets the geometric interpolation parameters for an equations interpolation.
  SUBROUTINE EquationsInterpolation_GeometricParametersGet(equationsInterpolation,variableType,geometricParameters,err,error,*)

    !Argument variables
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation !<A pointer to the equations interpolation to get the geometric parameters for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to get the geometric parameters for. \see FIELD_ROUTINES_VariableTypes
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: geometricParameters !<On exit, a pointer to the geometric interpolation parameters for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsInterpolation_GeometricParametersGet",err,error,*998)

    IF(ASSOCIATED(geometricParameters)) CALL FlagError("Geometric parameters is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
    IF(variableType<1.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    geometricParameters=>equationsInterpolation%geometricInterpParameters(variableType)%ptr
    IF(.NOT.ASSOCIATED(geometricParameters)) THEN
      localError="Equations interpolation geometric parameters is not associated for field variable type "// &
        & TRIM(NumberToVString(variableType,"*",err,error))//" of the equations interpolation."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("EquationsInterpolation_GeometricParametersGet")
    RETURN
999 NULLIFY(geometricParameters)
998 ERRORS("EquationsInterpolation_GeometricParametersGet",err,error)
    EXITS("EquationsInterpolation_GeometricParametersGet")
    RETURN 1
    
  END SUBROUTINE EquationsInterpolation_GeometricParametersGet

  !
  !================================================================================================================================
  !

  !>Gets the geometric interpolated point for an equations interpolation.
  SUBROUTINE EquationsInterpolation_GeometricPointGet(equationsInterpolation,variableType,geometricPoint,err,error,*)

    !Argument variables
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation !<A pointer to the equations interpolation to get the geometric point for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to get the geometric point for. \see FIELD_ROUTINES_VariableTypes
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: geometricPoint !<On exit, a pointer to the geometric interpolated point for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsInterpolation_GeometricPointGet",err,error,*998)

    IF(ASSOCIATED(geometricPoint)) CALL FlagError("Geometric point is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
    IF(variableType<1.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    geometricPoint=>equationsInterpolation%geometricInterpPoint(variableType)%ptr
    IF(.NOT.ASSOCIATED(geometricPoint)) THEN
      localError="Equations interpolated geometric point is not associated for field variable type "// &
        & TRIM(NumberToVString(variableType,"*",err,error))//" of the equations interpolation."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("EquationsInterpolation_GeometricPointGet")
    RETURN
999 NULLIFY(geometricPoint)
998 ERRORS("EquationsInterpolation_GeometricPointGet",err,error)
    EXITS("EquationsInterpolation_GeometricPointGet")
    RETURN 1
    
  END SUBROUTINE EquationsInterpolation_GeometricPointGet

  !
  !================================================================================================================================
  !

  !>Gets the geometric interpolated point metrics for an equations interpolation.
  SUBROUTINE EquationsInterpolation_GeometricPointMetricsGet(equationsInterpolation,variableType,geometricPointMetrics,err,error,*)

    !Argument variables
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation !<A pointer to the equations interpolation to get the geometric point metrics for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to get the geometric point metrics for. \see FIELD_ROUTINES_VariableTypes
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: geometricPointMetrics !<On exit, a pointer to the geometric interpolated point metrics for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsInterpolation_GeometricPointMetricsGet",err,error,*998)

    IF(ASSOCIATED(geometricPointMetrics)) CALL FlagError("Geometric point metrics is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
    IF(variableType<1.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    geometricPointMetrics=>equationsInterpolation%geometricInterpPointMetrics(variableType)%ptr
    IF(.NOT.ASSOCIATED(geometricPointMetrics)) THEN
      localError="Equations interpolated geometric point metrics is not associated for field variable type "// &
        & TRIM(NumberToVString(variableType,"*",err,error))//" of the equations interpolation."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("EquationsInterpolation_GeometricPointMetricsGet")
    RETURN
999 NULLIFY(geometricPointMetrics)
998 ERRORS("EquationsInterpolation_GeometricPointMetricsGet",err,error)
    EXITS("EquationsInterpolation_GeometricPointMetricsGet")
    RETURN 1
    
  END SUBROUTINE EquationsInterpolation_GeometricPointMetricsGet

  !
  !================================================================================================================================
  !

  !>Gets the independent interpolation parameters for an equations interpolation.
  SUBROUTINE EquationsInterpolation_IndependentParametersGet(equationsInterpolation,variableType,independentParameters,err,error,*)

    !Argument variables
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation !<A pointer to the equations interpolation to get the independent parameters for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to get the independent parameters for. \see FIELD_ROUTINES_VariableTypes
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: independentParameters !<On exit, a pointer to the independent interpolation parameters for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsInterpolation_IndependentParametersGet",err,error,*998)

    IF(ASSOCIATED(independentParameters)) CALL FlagError("Independent parameters is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
    IF(variableType<1.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    independentParameters=>equationsInterpolation%independentInterpParameters(variableType)%ptr
    IF(.NOT.ASSOCIATED(independentParameters)) THEN
      localError="Equations interpolation independent parameters is not associated for field variable type "// &
        & TRIM(NumberToVString(variableType,"*",err,error))//" of the equations interpolation."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("EquationsInterpolation_IndependentParametersGet")
    RETURN
999 NULLIFY(independentParameters)
998 ERRORS("EquationsInterpolation_IndependentParametersGet",err,error)
    EXITS("EquationsInterpolation_IndependentParametersGet")
    RETURN 1
    
  END SUBROUTINE EquationsInterpolation_IndependentParametersGet

  !
  !================================================================================================================================
  !

  !>Gets the independent interpolated point for an equations interpolation.
  SUBROUTINE EquationsInterpolation_IndependentPointGet(equationsInterpolation,variableType,independentPoint,err,error,*)

    !Argument variables
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation !<A pointer to the equations interpolation to get the independent point for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to get the independent point for. \see FIELD_ROUTINES_VariableTypes
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: independentPoint !<On exit, a pointer to the independent interpolated point for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsInterpolation_IndependentPointGet",err,error,*998)

    IF(ASSOCIATED(independentPoint)) CALL FlagError("Independent point is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
    IF(variableType<1.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    independentPoint=>equationsInterpolation%independentInterpPoint(variableType)%ptr
    IF(.NOT.ASSOCIATED(independentPoint)) THEN
      localError="Equations interpolated independent point is not associated for field variable type "// &
        & TRIM(NumberToVString(variableType,"*",err,error))//" of the equations interpolation."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("EquationsInterpolation_IndependentPointGet")
    RETURN
999 NULLIFY(independentPoint)
998 ERRORS("EquationsInterpolation_IndependentPointGet",err,error)
    EXITS("EquationsInterpolation_IndependentPointGet")
    RETURN 1
    
  END SUBROUTINE EquationsInterpolation_IndependentPointGet

  !
  !================================================================================================================================
  !

  !>Gets the independent interpolated point metrics for an equations interpolation.
  SUBROUTINE EquationsInterpolation_IndependentPointMetricsGet(equationsInterpolation,variableType,independentPointMetrics, &
    & err,error,*)

    !Argument variables
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation !<A pointer to the equations interpolation to get the independent point metrics for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to get the independent point metrics for. \see FIELD_ROUTINES_VariableTypes
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: independentPointMetrics !<On exit, a pointer to the independent interpolated point metrics for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsInterpolation_IndependentPointMetricsGet",err,error,*998)

    IF(ASSOCIATED(independentPointMetrics)) CALL FlagError("Independent point metrics is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
    IF(variableType<1.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    independentPointMetrics=>equationsInterpolation%independentInterpPointMetrics(variableType)%ptr
    IF(.NOT.ASSOCIATED(independentPointMetrics)) THEN
      localError="Equations interpolated independent point metrics is not associated for field variable type "// &
        & TRIM(NumberToVString(variableType,"*",err,error))//" of the equations interpolation."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("EquationsInterpolation_IndependentPointMetricsGet")
    RETURN
999 NULLIFY(independentPointMetrics)
998 ERRORS("EquationsInterpolation_IndependentPointMetricsGet",err,error)
    EXITS("EquationsInterpolation_IndependentPointMetricsGet")
    RETURN 1
    
  END SUBROUTINE EquationsInterpolation_IndependentPointMetricsGet

  !
  !================================================================================================================================
  !

  !>Gets the materials interpolation parameters for an equations interpolation.
  SUBROUTINE EquationsInterpolation_MaterialsParametersGet(equationsInterpolation,variableType,materialsParameters,err,error,*)

    !Argument variables
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation !<A pointer to the equations interpolation to get the materials parameters for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to get the materials parameters for. \see FIELD_ROUTINES_VariableTypes
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: materialsParameters !<On exit, a pointer to the materials interpolation parameters for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsInterpolation_MaterialsParametersGet",err,error,*998)

    IF(ASSOCIATED(materialsParameters)) CALL FlagError("Materials parameters is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
    IF(variableType<1.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    materialsParameters=>equationsInterpolation%materialsInterpParameters(variableType)%ptr
    IF(.NOT.ASSOCIATED(materialsParameters)) THEN
      localError="Equations interpolation materials parameters is not associated for field variable type "// &
        & TRIM(NumberToVString(variableType,"*",err,error))//" of the equations interpolation."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("EquationsInterpolation_MaterialsParametersGet")
    RETURN
999 NULLIFY(materialsParameters)
998 ERRORS("EquationsInterpolation_MaterialsParametersGet",err,error)
    EXITS("EquationsInterpolation_MaterialsParametersGet")
    RETURN 1
    
  END SUBROUTINE EquationsInterpolation_MaterialsParametersGet

  !
  !================================================================================================================================
  !

  !>Gets the materials interpolated point for an equations interpolation.
  SUBROUTINE EquationsInterpolation_MaterialsPointGet(equationsInterpolation,variableType,materialsPoint,err,error,*)

    !Argument variables
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation !<A pointer to the equations interpolation to get the materials point for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to get the materials point for. \see FIELD_ROUTINES_VariableTypes
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: materialsPoint !<On exit, a pointer to the materials interpolated point for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsInterpolation_MaterialsPointGet",err,error,*998)

    IF(ASSOCIATED(materialsPoint)) CALL FlagError("Materials point is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
    IF(variableType<1.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    materialsPoint=>equationsInterpolation%materialsInterpPoint(variableType)%ptr
    IF(.NOT.ASSOCIATED(materialsPoint)) THEN
      localError="Equations interpolated materials point is not associated for field variable type "// &
        & TRIM(NumberToVString(variableType,"*",err,error))//" of the equations interpolation."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("EquationsInterpolation_MaterialsPointGet")
    RETURN
999 NULLIFY(materialsPoint)
998 ERRORS("EquationsInterpolation_MaterialsPointGet",err,error)
    EXITS("EquationsInterpolation_MaterialsPointGet")
    RETURN 1
    
  END SUBROUTINE EquationsInterpolation_MaterialsPointGet

  !
  !================================================================================================================================
  !

  !>Gets the previous dependent interpolation parameters for an equations interpolation.
  SUBROUTINE EquationsInterpolation_PreviousDependentParametersGet(equationsInterpolation,variableType,prevDependentParameters, &
    & err,error,*)

    !Argument variables
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation !<A pointer to the equations interpolation to get the previous dependent parameters for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to get the previous dependent parameters for. \see FIELD_ROUTINES_VariableTypes
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: prevDependentParameters !<On exit, a pointer to the previous dependent interpolation parameters for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsInterpolation_PreviousDependentParametersGet",err,error,*998)

    IF(ASSOCIATED(prevDependentParameters)) CALL FlagError("Previous dependent parameters is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
    IF(variableType<1.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    prevDependentParameters=>equationsInterpolation%prevDependentInterpParameters(variableType)%ptr
    IF(.NOT.ASSOCIATED(prevDependentParameters)) THEN
      localError="Equations interpolation previous dependent parameters is not associated for field variable type "// &
        & TRIM(NumberToVString(variableType,"*",err,error))//" of the equations interpolation."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("EquationsInterpolation_PreviousDependentParametersGet")
    RETURN
999 NULLIFY(prevDependentParameters)
998 ERRORS("EquationsInterpolation_PreviousDependentParametersGet",err,error)
    EXITS("EquationsInterpolation_PreviousDependentParametersGet")
    RETURN 1
    
  END SUBROUTINE EquationsInterpolation_PreviousDependentParametersGet

  !
  !================================================================================================================================
  !

  !>Gets the previous dependent interpolated point for an equations interpolation.
  SUBROUTINE EquationsInterpolation_PreviousDependentPointGet(equationsInterpolation,variableType,prevDependentPoint,err,error,*)

    !Argument variables
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation !<A pointer to the equations interpolation to get the previous dependent point for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to get the prevous dependent point for. \see FIELD_ROUTINES_VariableTypes
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: prevDependentPoint !<On exit, a pointer to the previous dependent interpolated point for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsInterpolation_PreviousDependentPointGet",err,error,*998)

    IF(ASSOCIATED(prevDependentPoint)) CALL FlagError("Previous dependent point is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
    IF(variableType<1.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    prevDependentPoint=>equationsInterpolation%prevDependentInterpPoint(variableType)%ptr
    IF(.NOT.ASSOCIATED(prevDependentPoint)) THEN
      localError="Equations interpolated previous dependent point is not associated for field variable type "// &
        & TRIM(NumberToVString(variableType,"*",err,error))//" of the equations interpolation."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("EquationsInterpolation_PreviousDependentPointGet")
    RETURN
999 NULLIFY(prevDependentPoint)
998 ERRORS("EquationsInterpolation_PreviousDependentPointGet",err,error)
    EXITS("EquationsInterpolation_PreviousDependentPointGet")
    RETURN 1
    
  END SUBROUTINE EquationsInterpolation_PreviousDependentPointGet

  !
  !================================================================================================================================
  !

  !>Gets the previous dependent interpolated point metrics for an equations interpolation.
  SUBROUTINE EquationsInterpolation_PreviousDependentPointMetricsGet(equationsInterpolation,variableType, &
    & prevDependentPointMetrics,err,error,*)

    !Argument variables
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation !<A pointer to the equations interpolation to get the previous dependent point metrics for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to get the prevous dependent point for. \see FIELD_ROUTINES_VariableTypes
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: prevDependentPointMetrics !<On exit, a pointer to the previous dependent interpolated point metrics for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsInterpolation_PreviousDependentPointMetricsGet",err,error,*998)

    IF(ASSOCIATED(prevDependentPointMetrics)) CALL FlagError("Previous dependent point metrics is already associated.", &
      & err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
    IF(variableType<1.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    prevDependentPointMetrics=>equationsInterpolation%prevDependentInterpPointMetrics(variableType)%ptr
    IF(.NOT.ASSOCIATED(prevDependentPointMetrics)) THEN
      localError="Equations interpolated previous dependent point metrics is not associated for field variable type "// &
        & TRIM(NumberToVString(variableType,"*",err,error))//" of the equations interpolation."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("EquationsInterpolation_PreviousDependentPointMetricsGet")
    RETURN
999 NULLIFY(prevDependentPointMetrics)
998 ERRORS("EquationsInterpolation_PreviousDependentPointMetricsGet",err,error)
    EXITS("EquationsInterpolation_PreviousDependentPointMetricsGet")
    RETURN 1
    
  END SUBROUTINE EquationsInterpolation_PreviousDependentPointMetricsGet

  !
  !================================================================================================================================
  !

  !>Gets the source interpolation parameters for an equations interpolation.
  SUBROUTINE EquationsInterpolation_SourceParametersGet(equationsInterpolation,variableType,sourceParameters,err,error,*)

    !Argument variables
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation !<A pointer to the equations interpolation to get the source parameters for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to get the source parameters for. \see FIELD_ROUTINES_VariableTypes
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: sourceParameters !<On exit, a pointer to the source interpolation parameters for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsInterpolation_SourceParametersGet",err,error,*998)

    IF(ASSOCIATED(sourceParameters)) CALL FlagError("Source parameters is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
    IF(variableType<1.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    sourceParameters=>equationsInterpolation%sourceInterpParameters(variableType)%ptr
    IF(.NOT.ASSOCIATED(sourceParameters)) THEN
      localError="Equations interpolation source parameters is not associated for field variable type "// &
        & TRIM(NumberToVString(variableType,"*",err,error))//" of the equations interpolation."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("EquationsInterpolation_SourceParametersGet")
    RETURN
999 NULLIFY(sourceParameters)
998 ERRORS("EquationsInterpolation_SourceParametersGet",err,error)
    EXITS("EquationsInterpolation_SourceParametersGet")
    RETURN 1
    
  END SUBROUTINE EquationsInterpolation_SourceParametersGet

  !
  !================================================================================================================================
  !

  !>Gets the source interpolated point for an equations interpolation.
  SUBROUTINE EquationsInterpolation_SourcePointGet(equationsInterpolation,variableType,sourcePoint,err,error,*)

    !Argument variables
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation !<A pointer to the equations interpolation to get the source point for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to get the source point for. \see FIELD_ROUTINES_VariableTypes
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: sourcePoint !<On exit, a pointer to the source interpolated point for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsInterpolation_SourcePointGet",err,error,*998)

    IF(ASSOCIATED(sourcePoint)) CALL FlagError("Source point is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
    IF(variableType<1.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    sourcePoint=>equationsInterpolation%sourceInterpPoint(variableType)%ptr
    IF(.NOT.ASSOCIATED(sourcePoint)) THEN
      localError="Equations interpolated source point is not associated for field variable type "// &
        & TRIM(NumberToVString(variableType,"*",err,error))//" of the equations interpolation."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("EquationsInterpolation_SourcePointGet")
    RETURN
999 NULLIFY(sourcePoint)
998 ERRORS("EquationsInterpolation_SourcePointGet",err,error)
    EXITS("EquationsInterpolation_SourcePointGet")
    RETURN 1
    
  END SUBROUTINE EquationsInterpolation_SourcePointGet

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

  !>Gets the scalar mapping for a scalar equations.
  SUBROUTINE EquationsScalar_ScalarMappingGet(scalarEquations,scalarMapping,err,error,*)

    !Argument variables
    TYPE(EquationsScalarType), POINTER :: scalarEquations !<A pointer to the scalar equations to get the scalar mapping for
    TYPE(EquationsMappingScalarType), POINTER :: scalarMapping !<On exit, a pointer to the scalar mapping in the specified scalar equations. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsScalar_ScalarMappingGet",err,error,*998)

    IF(ASSOCIATED(scalarMapping)) CALL FlagError("Scalar mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(scalarEquations)) CALL FlagError("Scalar equations is not associated.",err,error,*999)

    scalarMapping=>scalarEquations%scalarMapping
    IF(.NOT.ASSOCIATED(scalarMapping)) CALL FlagError("Scalar mapping is not associated for the scalar equations.",err,error,*999)
       
    EXITS("EquationsScalar_ScalarMappingGet")
    RETURN
999 NULLIFY(scalarMapping)
998 ERRORSEXITS("EquationsScalar_ScalarMappingGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsScalar_ScalarMappingGet

  !
  !================================================================================================================================
  !

  !>Gets the scalar matrices for a scalar equations.
  SUBROUTINE EquationsScalar_ScalarMatricesGet(scalarEquations,scalarMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsScalarType), POINTER :: scalarEquations !<A pointer to the scalar equations to get the scalar matrices for
    TYPE(EquationsMatricesScalarType), POINTER :: scalarMatrices !<On exit, a pointer to the scalar matrices in the specified scalar equations. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsScalar_ScalarMatricesGet",err,error,*998)

    IF(ASSOCIATED(scalarMatrices)) CALL FlagError("Scalar matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(scalarEquations)) CALL FlagError("Scalar equations is not associated.",err,error,*999)

    scalarMatrices=>scalarEquations%scalarMatrices
    IF(.NOT.ASSOCIATED(scalarMatrices)) CALL FlagError("Scalar matrices is not associated for the scalar equations.",err,error,*999)
       
    EXITS("EquationsScalar_ScalarMatricesGet")
    RETURN
999 NULLIFY(scalarMatrices)
998 ERRORSEXITS("EquationsScalar_ScalarMatricesGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsScalar_ScalarMatricesGet

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

  !>Gets the vector mapping for a vector equations.
  SUBROUTINE EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*)

    !Argument variables
    TYPE(EquationsVectorType), POINTER :: vectorEquations !<A pointer to the vector equations to get the vector mapping for
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping !<On exit, a pointer to the vector mapping in the specified vector equations. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsVector_VectorMappingGet",err,error,*998)

    IF(ASSOCIATED(vectorMapping)) CALL FlagError("Vector mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorEquations)) CALL FlagError("Vector equations is not associated.",err,error,*999)

    vectorMapping=>vectorEquations%vectorMapping
    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector mapping is not associated for the vector equations.",err,error,*999)
       
    EXITS("EquationsVector_VectorMappingGet")
    RETURN
999 NULLIFY(vectorMapping)
998 ERRORSEXITS("EquationsVector_VectorMappingGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsVector_VectorMappingGet

  !
  !================================================================================================================================
  !

  !>Gets the vector matrices for a vector equations.
  SUBROUTINE EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsVectorType), POINTER :: vectorEquations !<A pointer to the vector equations to get the vector matrices for
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices !<On exit, a pointer to the vector matrices in the specified vector equations. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsVector_VectorMatricesGet",err,error,*998)

    IF(ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorEquations)) CALL FlagError("Vector equations is not associated.",err,error,*999)

    vectorMatrices=>vectorEquations%vectorMatrices
    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is not associated for the vector equations.",err,error,*999)
       
    EXITS("EquationsVector_VectorMatricesGet")
    RETURN
999 NULLIFY(vectorMatrices)
998 ERRORSEXITS("EquationsVector_VectorMatricesGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsVector_VectorMatricesGet

  !
  !================================================================================================================================
  !

END MODULE EquationsAccessRoutines
