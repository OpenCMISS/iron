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
  USE ISO_VARYING_STRING
  USE Kinds
  USE Strings
  USE Types

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup EquationsRoutines_LinearityTypes Equations::Constants::LinearityTypes
  !> \brief The equations linearity types
  !> \see EquationsRoutines,OpenCMISS_EquationsLinearityTypes
  !>@{
  INTEGER(INTG), PARAMETER :: NUMBER_OF_EQUATIONS_LINEARITY_TYPES=3 !<The number of equations linearity types defined. \see EquationsRoutines_LinearityTypes,EquationsRoutines
  INTEGER(INTG), PARAMETER :: EQUATIONS_LINEAR=1 !<The equations are linear. \see EquationsRoutines_LinearityTypes,EquationsRoutines
  INTEGER(INTG), PARAMETER :: EQUATIONS_NONLINEAR=2 !<The equations are non-linear. \see EquationsRoutines_LinearityTypes,EquationsRoutines
  INTEGER(INTG), PARAMETER :: EQUATIONS_NONLINEAR_BCS=3 !<The equations have non-linear boundary conditions. \see EquationsRoutines_LinearityTypes,EquationsRoutines
  !>@}

 
  !> \addtogroup EquationsRoutines_TimeDependenceTypes Equations::Constants::TimeDependenceTypes
  !> \brief The equations time dependence type parameters
  !> \see EquationsRoutines,OpenCMISS_EquationsTimeDependenceTypes
  !>@{
  INTEGER(INTG), PARAMETER :: NUMBER_OF_EQUATIONS_TIME_TYPES=5 !<The number of equations time dependence types defined. \see EquationsRoutines_TimeDependenceTypes,EquationRoutines
  INTEGER(INTG), PARAMETER :: EQUATIONS_STATIC=1 !<The equations are static and have no time dependence. \see EquationsRoutines_TimeDependenceTypes,EquationsRoutines
  INTEGER(INTG), PARAMETER :: EQUATIONS_QUASISTATIC=2 !<The equations are quasi-static. \see EquationsRoutines_TimeDependenceTypes,EquationsRoutines
  INTEGER(INTG), PARAMETER :: EQUATIONS_FIRST_ORDER_DYNAMIC=3 !<The equations are first order dynamic. \see EquationsRoutines_TimeDependenceTypes,EquationsRoutines
  INTEGER(INTG), PARAMETER :: EQUATIONS_SECOND_ORDER_DYNAMIC=4 !<The equations are a second order dynamic. \see EquationsRoutines_TimeDependenceTypes,EquationsRoutines
  INTEGER(INTG), PARAMETER :: EQUATIONS_TIME_STEPPING=5 !<The equations are for time stepping. \see EquationsRoutines_TimeDependenceTypes,EquationsRoutines
  !>@}

  !Module types

  !Module variables

  !Interfaces

  PUBLIC NUMBER_OF_EQUATIONS_LINEARITY_TYPES,EQUATIONS_LINEAR,EQUATIONS_NONLINEAR,EQUATIONS_NONLINEAR_BCS

  PUBLIC NUMBER_OF_EQUATIONS_TIME_TYPES,EQUATIONS_STATIC,EQUATIONS_QUASISTATIC,EQUATIONS_FIRST_ORDER_DYNAMIC, &
    & EQUATIONS_SECOND_ORDER_DYNAMIC,EQUATIONS_TIME_STEPPING

  PUBLIC Equations_AssertIsFinished,Equations_AssertNotFinished

  PUBLIC Equations_AssertIsDynamic,Equations_AssertIsStatic

  PUBLIC Equations_AssertIsLinear,Equations_AssertIsNonlinear

  PUBLIC Equations_EquationsSetGet
  
  PUBLIC Equations_InterpolationGet

  PUBLIC Equations_LinearityTypeGet
  
  PUBLIC Equations_LumpingTypeGet
  
  PUBLIC Equations_OutputTypeGet
  
  PUBLIC Equations_ScalarEquationsGet

  PUBLIC Equations_SparsityTypeGet
  
  PUBLIC Equations_TimeDependenceTypeGet
  
  PUBLIC Equations_VectorEquationsGet

  PUBLIC EquationsInterpolation_DependentFieldExists
  
  PUBLIC EquationsInterpolation_DependentFieldGet

  PUBLIC EquationsInterpolation_DependentParametersGet
  
  PUBLIC EquationsInterpolation_DependentPointGet
  
  PUBLIC EquationsInterpolation_DependentPointMetricsGet
  
  PUBLIC EquationsInterpolation_FibreFieldExists
  
  PUBLIC EquationsInterpolation_FibreFieldGet

  PUBLIC EquationsInterpolation_FibreParametersGet
  
  PUBLIC EquationsInterpolation_FibrePointGet
  
  PUBLIC EquationsInterpolation_FibrePointMetricsGet
  
  PUBLIC EquationsInterpolation_GeometricFieldExists
  
  PUBLIC EquationsInterpolation_GeometricFieldGet

  PUBLIC EquationsInterpolation_GeometricParametersGet
  
  PUBLIC EquationsInterpolation_GeometricPointGet
  
  PUBLIC EquationsInterpolation_GeometricPointMetricsGet
  
  PUBLIC EquationsInterpolation_IndependentFieldExists
  
  PUBLIC EquationsInterpolation_IndependentFieldGet

  PUBLIC EquationsInterpolation_IndependentParametersGet
  
  PUBLIC EquationsInterpolation_IndependentPointGet
  
  PUBLIC EquationsInterpolation_IndependentPointMetricsGet
  
  PUBLIC EquationsInterpolation_MaterialsFieldExists
  
  PUBLIC EquationsInterpolation_MaterialsFieldGet

  PUBLIC EquationsInterpolation_MaterialsParametersGet
  
  PUBLIC EquationsInterpolation_MaterialsPointGet
  
  PUBLIC EquationsInterpolation_PreviousDependentParametersGet
  
  PUBLIC EquationsInterpolation_PreviousDependentPointGet
  
  PUBLIC EquationsInterpolation_PreviousDependentPointMetricsGet
  
  PUBLIC EquationsInterpolation_SourceFieldExists
  
  PUBLIC EquationsInterpolation_SourceFieldGet

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

  !>Assert that an equations has been finished
  SUBROUTINE Equations_AssertIsFinished(equations,err,error,*)

    !Argument Variables
    TYPE(EquationsType), POINTER, INTENT(IN) :: equations !<The equations to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("Equations_AssertIsFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
#endif    

    IF(.NOT.equations%equationsFinished) CALL FlagError("Equations has not been finished.",err,error,*999)
    
    EXITS("Equations_AssertIsFinished")
    RETURN
999 ERRORSEXITS("Equations_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE Equations_AssertIsFinished

  !
  !================================================================================================================================
  !

  !>Assert that an equations has not been finished
  SUBROUTINE Equations_AssertNotFinished(equations,err,error,*)

    !Argument Variables
    TYPE(EquationsType), POINTER, INTENT(IN) :: equations !<The equations to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Equations_AssertNotFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
#endif    

    IF(equations%equationsFinished) CALL FlagError("Equations has already been finished.",err,error,*999)
    
    EXITS("Equations_AssertNotFinished")
    RETURN
999 ERRORSEXITS("Equations_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE Equations_AssertNotFinished

  !
  !================================================================================================================================
  !

  !>Assert that an equations time dependence is dynamic
  SUBROUTINE Equations_AssertIsDynamic(equations,err,error,*)

    !Argument Variables
    TYPE(EquationsType), POINTER, INTENT(IN) :: equations !<The equations to assert the dynamic status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Equations_AssertIsDynamic",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
#endif    

    IF(equations%timeDependence/=EQUATIONS_FIRST_ORDER_DYNAMIC.AND. &
      & equations%timeDependence/=EQUATIONS_SECOND_ORDER_DYNAMIC) THEN
      localError="The equations time dependence type of "// &
        & TRIM(NumberToVString(equations%timeDependence,"*",err,error))// &
        & " does not correspond to first or second order dynamic equations."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Equations_AssertIsDynamic")
    RETURN
999 ERRORSEXITS("Equations_AssertIsDynamic",err,error)
    RETURN 1
    
  END SUBROUTINE Equations_AssertIsDynamic

  !
  !================================================================================================================================
  !

  !>Assert that an equations time dependence is static
  SUBROUTINE Equations_AssertIsStatic(equations,err,error,*)

    !Argument Variables
    TYPE(EquationsType), POINTER, INTENT(IN) :: equations !<The equations to assert the static status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Equations_AssertIsStatic",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
#endif    

    IF(equations%timeDependence/=EQUATIONS_STATIC.AND. &
      & equations%timeDependence/=EQUATIONS_QUASISTATIC) THEN
      localError="The equations time dependence type of "// &
        & TRIM(NumberToVString(equations%timeDependence,"*",err,error))// &
        & " does not correspond to static or quasi static equations."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Equations_AssertIsStatic")
    RETURN
999 ERRORSEXITS("Equations_AssertIsStatic",err,error)
    RETURN 1
    
  END SUBROUTINE Equations_AssertIsStatic

  !
  !================================================================================================================================
  !

  !>Assert that an equations linearity is linear
  SUBROUTINE Equations_AssertIsLinear(equations,err,error,*)

    !Argument Variables
    TYPE(EquationsType), POINTER, INTENT(IN) :: equations !<The equations to assert the linear status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Equations_AssertIsLinear",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
#endif    

    IF(equations%linearity/=EQUATIONS_LINEAR) THEN
      localError="The equations linearity type of "//TRIM(NumberToVString(equations%linearity,"*",err,error))// &
        & " does not correspond to linear equations."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Equations_AssertIsLinear")
    RETURN
999 ERRORSEXITS("Equations_AssertIsLinear",err,error)
    RETURN 1
    
  END SUBROUTINE Equations_AssertIsLinear

  !
  !================================================================================================================================
  !

  !>Assert that an equations linearity is nonlinear
  SUBROUTINE Equations_AssertIsNonlinear(equations,err,error,*)

    !Argument Variables
    TYPE(EquationsType), POINTER, INTENT(IN) :: equations !<The equations to assert the nonlinear status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Equations_AssertIsNonlinear",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
#endif    

    IF(equations%linearity/=EQUATIONS_NONLINEAR) THEN
      localError="The equations linearity type of "//TRIM(NumberToVString(equations%linearity,"*",err,error))// &
        & " does not correspond to nonlinear equations."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Equations_AssertIsNonlinear")
    RETURN
999 ERRORSEXITS("Equations_AssertIsNonlinear",err,error)
    RETURN 1
    
  END SUBROUTINE Equations_AssertIsNonlinear

  !
  !================================================================================================================================
  !

  !>Gets the equations set for an equations.
  SUBROUTINE Equations_EquationsSetGet(equations,equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER :: equations !<A pointer to the equations to get the equations set for
    TYPE(EquationsSetType), POINTER :: equationsSet !<On exit, a pointer to the equations set for the specified equations. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Equations_EquationsSetGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(equationsSet)) CALL FlagError("Equations set is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
#endif    

    equationsSet=>equations%equationsSet

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated for the equations.",err,error,*999)
#endif    
       
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

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(interpolation)) CALL FlagError("Interpolation is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
#endif

    interpolation=>equations%interpolation

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(interpolation)) CALL FlagError("Interpolation is not associated for the equations.",err,error,*999)
#endif    
       
    EXITS("Equations_InterpolationGet")
    RETURN
999 NULLIFY(interpolation)
998 ERRORSEXITS("Equations_InterpolationGet",err,error)
    RETURN 1
    
  END SUBROUTINE Equations_InterpolationGet

  !
  !================================================================================================================================
  !

  !>Gets the linearity type for an equations.
  SUBROUTINE Equations_LinearityTypeGet(equations,linearityType,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER :: equations !<A pointer to the equations to get the linearity type for
    INTEGER(INTG), INTENT(OUT) :: linearityType !<On exit, the linearity type for the specified equations. \see EquationsSetRoutines_LinearityTypes,EquationsSetRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Equations_LinearityTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
#endif

    linearityType=equations%linearity

    EXITS("Equations_LinearityTypeGet")
    RETURN
999 ERRORSEXITS("Equations_LinearityTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE Equations_LinearityTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the lumping type for an equations.
  SUBROUTINE Equations_LumpingTypeGet(equations,lumpingType,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER :: equations !<A pointer to the equations to get the lumping type for
    INTEGER(INTG), INTENT(OUT) :: lumpingType !<On exit, the lumping type for the specified equations. \see EquationsRoutines_EquationsLumpingTypes,EquationsRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Equations_LumpingTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
#endif

    lumpingType=equations%lumpingType

    EXITS("Equations_LumpingTypeGet")
    RETURN
999 ERRORSEXITS("Equations_LumpingTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE Equations_LumpingTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the output type for an equations.
  SUBROUTINE Equations_OutputTypeGet(equations,outputType,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER :: equations !<A pointer to the equations to get the output type for
    INTEGER(INTG), INTENT(OUT) :: outputType !<On exit, the output type for the specified equations.  \see EquationsRoutines_EquationsOutputTypes,EquationsRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Equations_OutputTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
#endif

    outputType=equations%outputType

    EXITS("Equations_OutputTypeGet")
    RETURN
999 ERRORSEXITS("Equations_OutputTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE Equations_OutputTypeGet

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

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(scalarEquations)) CALL FlagError("Scalar equations is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
#endif    

    scalarEquations=>equations%scalarEquations

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(scalarEquations)) CALL FlagError("Scalar equations is not associated for the equations.",err,error,*999)
#endif    
       
    EXITS("Equations_ScalarEquationsGet")
    RETURN
999 NULLIFY(scalarEquations)
998 ERRORSEXITS("Equations_ScalarEquationsGet",err,error)
    RETURN 1
    
  END SUBROUTINE Equations_ScalarEquationsGet

  !
  !================================================================================================================================
  !

  !>Gets the sparsity type for an equations.
  SUBROUTINE Equations_SparsityTypeGet(equations,sparsityType,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER :: equations !<A pointer to the equations to get the sparsity type for
    INTEGER(INTG), INTENT(OUT) :: sparsityType !<On exit, the sparsity type for the specified equations. \see EquationsRoutines_EquationsSparsityTypes,EquationsRoutines 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Equations_SparsityTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
#endif

    sparsityType=equations%sparsityType

    EXITS("Equations_SparsityTypeGet")
    RETURN
999 ERRORSEXITS("Equations_SparsityTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE Equations_SparsityTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the time dependence type for an equations.
  SUBROUTINE Equations_TimeDependenceTypeGet(equations,timeDependenceType,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER :: equations !<A pointer to the equations to get the time dependence type for
    INTEGER(INTG), INTENT(OUT) :: timeDependenceType !<On exit, the time dependence type for the specified equations. \see EquationsSetRoutines_TimeDependenceTypes,EquationsSetRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Equations_TimeDependenceTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
#endif

    timeDependenceType=equations%timeDependence

    EXITS("Equations_TimeDependenceTypeGet")
    RETURN
999 ERRORSEXITS("Equations_TimeDependenceTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE Equations_TimeDependenceTypeGet

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

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(vectorEquations)) CALL FlagError("Vector equations is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
#endif    

    vectorEquations=>equations%vectorEquations

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(vectorEquations)) CALL FlagError("Vector equations is not associated for the equations.",err,error,*999)
#endif    
       
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
    TYPE(FieldInterpolationParametersType), POINTER :: dependentParameters !<On exit, a pointer to the dependent interpolation parameters for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("EquationsInterpolation_DependentParametersGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dependentParameters)) CALL FlagError("Dependent parameters is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equationsInterpolation%dependentInterpParameters)) &
      & CALL FlagError("Equations interpolation dependent interpolation parameters is not associated.",err,error,*999)
    IF(variableType<1.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    dependentParameters=>equationsInterpolation%dependentInterpParameters(variableType)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(dependentParameters)) THEN
      localError="Equations interpolation dependent parameters is not associated for field variable type "// &
        & TRIM(NumberToVString(variableType,"*",err,error))//" of the equations interpolation."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
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

  !>Checks the dependent field for an equations interpolation.
  SUBROUTINE EquationsInterpolation_DependentFieldExists(equationsInterpolation,dependentField,err,error,*)

    !Argument variables
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation !<A pointer to the equations interpolation to get the dependent field for
    TYPE(FieldType), POINTER :: dependentField !<On exit, a pointer to the dependent field for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsInterpolation_DependentFieldExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dependentField)) CALL FlagError("Dependent field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
#endif    

    dependentField=>equationsInterpolation%dependentField
       
    EXITS("EquationsInterpolation_DependentFieldExists")
    RETURN
999 NULLIFY(dependentField)
998 ERRORS("EquationsInterpolation_DependentFieldExists",err,error)
    EXITS("EquationsInterpolation_DependentFieldExists")
    RETURN 1
    
  END SUBROUTINE EquationsInterpolation_DependentFieldExists

  !
  !================================================================================================================================
  !

  !>Gets the dependent field for an equations interpolation.
  SUBROUTINE EquationsInterpolation_DependentFieldGet(equationsInterpolation,dependentField,err,error,*)

    !Argument variables
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation !<A pointer to the equations interpolation to get the dependent field for
    TYPE(FieldType), POINTER :: dependentField !<On exit, a pointer to the dependent field for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsInterpolation_DependentFieldGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dependentField)) CALL FlagError("Dependent field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
#endif    

    dependentField=>equationsInterpolation%dependentField

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(dependentField)) &
      & CALL FlagError("Dependent field is not assocaited for the equations interpolation.",err,error,*999)
#endif    
       
    EXITS("EquationsInterpolation_DependentFieldGet")
    RETURN
999 NULLIFY(dependentField)
998 ERRORS("EquationsInterpolation_DependentFieldGet",err,error)
    EXITS("EquationsInterpolation_DependentFieldGet")
    RETURN 1
    
  END SUBROUTINE EquationsInterpolation_DependentFieldGet

  !
  !================================================================================================================================
  !

  !>Gets the dependent interpolated point for an equations interpolation.
  SUBROUTINE EquationsInterpolation_DependentPointGet(equationsInterpolation,variableType,dependentPoint,err,error,*)

    !Argument variables
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation !<A pointer to the equations interpolation to get the dependent point for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to get the dependent point for. \see FIELD_ROUTINES_VariableTypes
    TYPE(FieldInterpolatedPointType), POINTER :: dependentPoint !<On exit, a pointer to the dependent interpolated point for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("EquationsInterpolation_DependentPointGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dependentPoint)) CALL FlagError("Dependent point is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equationsInterpolation%dependentInterpPoint)) &
      & CALL FlagError("Equations interpolation dependent interpolated point is not associated.",err,error,*999)
    IF(variableType<1.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    dependentPoint=>equationsInterpolation%dependentInterpPoint(variableType)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(dependentPoint)) THEN
      localError="Equations interpolated dependent point is not associated for field variable type "// &
        & TRIM(NumberToVString(variableType,"*",err,error))//" of the equations interpolation."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
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
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: dependentPointMetrics !<On exit, a pointer to the dependent interpolated point metrics for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("EquationsInterpolation_DependentPointMetricsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dependentPointMetrics)) CALL FlagError("Dependent point metrics is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equationsInterpolation%dependentInterpPointMetrics)) &
      & CALL FlagError("Equations interpolation dependent interpolated point metrics is not associated.",err,error,*999)
    IF(variableType<1.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    dependentPointMetrics=>equationsInterpolation%dependentInterpPointMetrics(variableType)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(dependentPointMetrics)) THEN
      localError="Equations interpolated dependent point metrics is not associated for field variable type "// &
        & TRIM(NumberToVString(variableType,"*",err,error))//" of the equations interpolation."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
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

  !>Checks the fibre field for an equations interpolation.
  SUBROUTINE EquationsInterpolation_FibreFieldExists(equationsInterpolation,fibreField,err,error,*)

    !Argument variables
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation !<A pointer to the equations interpolation to get the fibre field for
    TYPE(FieldType), POINTER :: fibreField !<On exit, a pointer to the fibre field for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsInterpolation_FibreFieldExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(fibreField)) CALL FlagError("Fibre field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
#endif    

    fibreField=>equationsInterpolation%fibreField
       
    EXITS("EquationsInterpolation_FibreFieldExists")
    RETURN
999 NULLIFY(fibreField)
998 ERRORS("EquationsInterpolation_FibreFieldExists",err,error)
    EXITS("EquationsInterpolation_FibreFieldExists")
    RETURN 1
    
  END SUBROUTINE EquationsInterpolation_FibreFieldExists

  !
  !================================================================================================================================
  !

  !>Gets the fibre field for an equations interpolation.
  SUBROUTINE EquationsInterpolation_FibreFieldGet(equationsInterpolation,fibreField,err,error,*)

    !Argument variables
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation !<A pointer to the equations interpolation to get the fibre field for
    TYPE(FieldType), POINTER :: fibreField !<On exit, a pointer to the fibre field for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsInterpolation_FibreFieldGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(fibreField)) CALL FlagError("Fibre field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
#endif    

    fibreField=>equationsInterpolation%fibreField

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(fibreField)) &
      & CALL FlagError("Fibre field is not assocaited for the equations interpolation.",err,error,*999)
#endif    
       
    EXITS("EquationsInterpolation_FibreFieldGet")
    RETURN
999 NULLIFY(fibreField)
998 ERRORS("EquationsInterpolation_FibreFieldGet",err,error)
    EXITS("EquationsInterpolation_FibreFieldGet")
    RETURN 1
    
  END SUBROUTINE EquationsInterpolation_FibreFieldGet

  !
  !================================================================================================================================
  !

  !>Gets the fibre interpolation parameters for an equations interpolation.
  SUBROUTINE EquationsInterpolation_FibreParametersGet(equationsInterpolation,variableType,fibreParameters,err,error,*)

    !Argument variables
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation !<A pointer to the equations interpolation to get the fibre parameters for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to get the fibre parameters for. \see FIELD_ROUTINES_VariableTypes
    TYPE(FieldInterpolationParametersType), POINTER :: fibreParameters !<On exit, a pointer to the fibre interpolation parameters for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("EquationsInterpolation_FibreParametersGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(fibreParameters)) CALL FlagError("Fibre parameters is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equationsInterpolation%fibreInterpParameters)) &
      & CALL FlagError("Equations interpolation fibre interpolation parameters is not associated.",err,error,*999)
    IF(variableType<1.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    fibreParameters=>equationsInterpolation%fibreInterpParameters(variableType)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(fibreParameters)) THEN
      localError="Equations interpolation fibre parameters is not associated for field variable type "// &
        & TRIM(NumberToVString(variableType,"*",err,error))//" of the equations interpolation."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
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
    TYPE(FieldInterpolatedPointType), POINTER :: fibrePoint !<On exit, a pointer to the fibre interpolated point for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("EquationsInterpolation_FibrePointGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(fibrePoint)) CALL FlagError("Fibre point is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equationsInterpolation%fibreInterpPoint)) &
      & CALL FlagError("Equations interpolation fibre interpolated point is not associated.",err,error,*999)
    IF(variableType<1.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    fibrePoint=>equationsInterpolation%fibreInterpPoint(variableType)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(fibrePoint)) THEN
      localError="Equations interpolated fibre point is not associated for field variable type "// &
        & TRIM(NumberToVString(variableType,"*",err,error))//" of the equations interpolation."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
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
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: fibrePointMetrics !<On exit, a pointer to the fibre interpolated point metrics for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("EquationsInterpolation_FibrePointMetricsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(fibrePointMetrics)) CALL FlagError("Fibre point metrics is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equationsInterpolation%fibreInterpPointMetrics)) &
      & CALL FlagError("Equations interpolation fibre interpolated point metrics is not associated.",err,error,*999)
    IF(variableType<1.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    fibrePointMetrics=>equationsInterpolation%fibreInterpPointMetrics(variableType)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(fibrePointMetrics)) THEN
      localError="Equations interpolated fibre point metrics is not associated for field variable type "// &
        & TRIM(NumberToVString(variableType,"*",err,error))//" of the equations interpolation."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
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

  !>Checks the geometric field for an equations interpolation.
  SUBROUTINE EquationsInterpolation_GeometricFieldExists(equationsInterpolation,geometricField,err,error,*)

    !Argument variables
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation !<A pointer to the equations interpolation to get the geometric field for
    TYPE(FieldType), POINTER :: geometricField !<On exit, a pointer to the geometric field for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsInterpolation_GeometricFieldExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(geometricField)) CALL FlagError("Geometric field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
#endif    

    geometricField=>equationsInterpolation%geometricField
       
    EXITS("EquationsInterpolation_GeometricFieldExists")
    RETURN
999 NULLIFY(geometricField)
998 ERRORS("EquationsInterpolation_GeometricFieldExists",err,error)
    EXITS("EquationsInterpolation_GeometricFieldExists")
    RETURN 1
    
  END SUBROUTINE EquationsInterpolation_GeometricFieldExists

  !
  !================================================================================================================================
  !

  !>Gets the geometric field for an equations interpolation.
  SUBROUTINE EquationsInterpolation_GeometricFieldGet(equationsInterpolation,geometricField,err,error,*)

    !Argument variables
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation !<A pointer to the equations interpolation to get the geometric field for
    TYPE(FieldType), POINTER :: geometricField !<On exit, a pointer to the geometric field for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsInterpolation_GeometricFieldGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(geometricField)) CALL FlagError("Geometric field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
#endif    

    geometricField=>equationsInterpolation%geometricField

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(geometricField)) &
      & CALL FlagError("Geometric field is not assocaited for the equations interpolation.",err,error,*999)
#endif    
       
    EXITS("EquationsInterpolation_GeometricFieldGet")
    RETURN
999 NULLIFY(geometricField)
998 ERRORS("EquationsInterpolation_GeometricFieldGet",err,error)
    EXITS("EquationsInterpolation_GeometricFieldGet")
    RETURN 1
    
  END SUBROUTINE EquationsInterpolation_GeometricFieldGet

  !
  !================================================================================================================================
  !

  !>Gets the geometric interpolation parameters for an equations interpolation.
  SUBROUTINE EquationsInterpolation_GeometricParametersGet(equationsInterpolation,variableType,geometricParameters,err,error,*)

    !Argument variables
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation !<A pointer to the equations interpolation to get the geometric parameters for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to get the geometric parameters for. \see FIELD_ROUTINES_VariableTypes
    TYPE(FieldInterpolationParametersType), POINTER :: geometricParameters !<On exit, a pointer to the geometric interpolation parameters for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("EquationsInterpolation_GeometricParametersGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(geometricParameters)) CALL FlagError("Geometric parameters is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equationsInterpolation%geometricInterpParameters)) &
      & CALL FlagError("Equations interpolation geometric interpolation parameters is not associated.",err,error,*999)
    IF(variableType<1.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    geometricParameters=>equationsInterpolation%geometricInterpParameters(variableType)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(geometricParameters)) THEN
      localError="Equations interpolation geometric parameters is not associated for field variable type "// &
        & TRIM(NumberToVString(variableType,"*",err,error))//" of the equations interpolation."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
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
    TYPE(FieldInterpolatedPointType), POINTER :: geometricPoint !<On exit, a pointer to the geometric interpolated point for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("EquationsInterpolation_GeometricPointGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(geometricPoint)) CALL FlagError("Geometric point is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equationsInterpolation%geometricInterpPoint)) &
      & CALL FlagError("Equations interpolation geometric interpolated point is not associated.",err,error,*999)
    IF(variableType<1.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    geometricPoint=>equationsInterpolation%geometricInterpPoint(variableType)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(geometricPoint)) THEN
      localError="Equations interpolated geometric point is not associated for field variable type "// &
        & TRIM(NumberToVString(variableType,"*",err,error))//" of the equations interpolation."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
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
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricPointMetrics !<On exit, a pointer to the geometric interpolated point metrics for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("EquationsInterpolation_GeometricPointMetricsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(geometricPointMetrics)) CALL FlagError("Geometric point metrics is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equationsInterpolation%geometricInterpPointMetrics)) &
      & CALL FlagError("Equations interpolation geometric interpolated point metrics is not associated.",err,error,*999)
    IF(variableType<1.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    geometricPointMetrics=>equationsInterpolation%geometricInterpPointMetrics(variableType)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(geometricPointMetrics)) THEN
      localError="Equations interpolated geometric point metrics is not associated for field variable type "// &
        & TRIM(NumberToVString(variableType,"*",err,error))//" of the equations interpolation."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
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

  !>Checks the independent field for an equations interpolation.
  SUBROUTINE EquationsInterpolation_IndependentFieldExists(equationsInterpolation,independentField,err,error,*)

    !Argument variables
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation !<A pointer to the equations interpolation to get the independent field for
    TYPE(FieldType), POINTER :: independentField !<On exit, a pointer to the independent field for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsInterpolation_IndependentFieldExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(independentField)) CALL FlagError("Independent field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
#endif    

    independentField=>equationsInterpolation%independentField
       
    EXITS("EquationsInterpolation_IndependentFieldExists")
    RETURN
999 NULLIFY(independentField)
998 ERRORS("EquationsInterpolation_IndependentFieldExists",err,error)
    EXITS("EquationsInterpolation_IndependentFieldExists")
    RETURN 1
    
  END SUBROUTINE EquationsInterpolation_IndependentFieldExists

  !
  !================================================================================================================================
  !

  !>Gets the independent field for an equations interpolation.
  SUBROUTINE EquationsInterpolation_IndependentFieldGet(equationsInterpolation,independentField,err,error,*)

    !Argument variables
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation !<A pointer to the equations interpolation to get the independent field for
    TYPE(FieldType), POINTER :: independentField !<On exit, a pointer to the independent field for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsInterpolation_IndependentFieldGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(independentField)) CALL FlagError("Independent field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
#endif    

    independentField=>equationsInterpolation%independentField

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(independentField)) &
      & CALL FlagError("Independent field is not assocaited for the equations interpolation.",err,error,*999)
#endif    
       
    EXITS("EquationsInterpolation_IndependentFieldGet")
    RETURN
999 NULLIFY(independentField)
998 ERRORS("EquationsInterpolation_IndependentFieldGet",err,error)
    EXITS("EquationsInterpolation_IndependentFieldGet")
    RETURN 1
    
  END SUBROUTINE EquationsInterpolation_IndependentFieldGet

  !
  !================================================================================================================================
  !

  !>Gets the independent interpolation parameters for an equations interpolation.
  SUBROUTINE EquationsInterpolation_IndependentParametersGet(equationsInterpolation,variableType,independentParameters,err,error,*)

    !Argument variables
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation !<A pointer to the equations interpolation to get the independent parameters for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to get the independent parameters for. \see FIELD_ROUTINES_VariableTypes
    TYPE(FieldInterpolationParametersType), POINTER :: independentParameters !<On exit, a pointer to the independent interpolation parameters for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("EquationsInterpolation_IndependentParametersGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(independentParameters)) CALL FlagError("Independent parameters is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equationsInterpolation%independentInterpParameters)) &
      & CALL FlagError("Equations interpolation independent interpolation parameters is not associated.",err,error,*999)
    IF(variableType<1.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    independentParameters=>equationsInterpolation%independentInterpParameters(variableType)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(independentParameters)) THEN
      localError="Equations interpolation independent parameters is not associated for field variable type "// &
        & TRIM(NumberToVString(variableType,"*",err,error))//" of the equations interpolation."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
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
    TYPE(FieldInterpolatedPointType), POINTER :: independentPoint !<On exit, a pointer to the independent interpolated point for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("EquationsInterpolation_IndependentPointGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(independentPoint)) CALL FlagError("Independent point is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equationsInterpolation%independentInterpPoint)) &
      & CALL FlagError("Equations interpolation independent interpolated point is not associated.",err,error,*999)
    IF(variableType<1.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    independentPoint=>equationsInterpolation%independentInterpPoint(variableType)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(independentPoint)) THEN
      localError="Equations interpolated independent point is not associated for field variable type "// &
        & TRIM(NumberToVString(variableType,"*",err,error))//" of the equations interpolation."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
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
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: independentPointMetrics !<On exit, a pointer to the independent interpolated point metrics for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("EquationsInterpolation_IndependentPointMetricsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(independentPointMetrics)) CALL FlagError("Independent point metrics is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equationsInterpolation%independentInterpPointMetrics)) &
      & CALL FlagError("Equations interpolation independent interpolated point metrics is not associated.",err,error,*999)
    IF(variableType<1.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    independentPointMetrics=>equationsInterpolation%independentInterpPointMetrics(variableType)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(independentPointMetrics)) THEN
      localError="Equations interpolated independent point metrics is not associated for field variable type "// &
        & TRIM(NumberToVString(variableType,"*",err,error))//" of the equations interpolation."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
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

  !>Checks the materials field for an equations interpolation.
  SUBROUTINE EquationsInterpolation_MaterialsFieldExists(equationsInterpolation,materialsField,err,error,*)

    !Argument variables
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation !<A pointer to the equations interpolation to get the materials field for
    TYPE(FieldType), POINTER :: materialsField !<On exit, a pointer to the materials field for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsInterpolation_MaterialsFieldExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(materialsField)) CALL FlagError("Materials field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
#endif    

    materialsField=>equationsInterpolation%materialsField
       
    EXITS("EquationsInterpolation_MaterialsFieldExists")
    RETURN
999 NULLIFY(materialsField)
998 ERRORS("EquationsInterpolation_MaterialsFieldExists",err,error)
    EXITS("EquationsInterpolation_MaterialsFieldExists")
    RETURN 1
    
  END SUBROUTINE EquationsInterpolation_MaterialsFieldExists

  !
  !================================================================================================================================
  !

  !>Gets the materials field for an equations interpolation.
  SUBROUTINE EquationsInterpolation_MaterialsFieldGet(equationsInterpolation,materialsField,err,error,*)

    !Argument variables
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation !<A pointer to the equations interpolation to get the materials field for
    TYPE(FieldType), POINTER :: materialsField !<On exit, a pointer to the materials field for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsInterpolation_MaterialsFieldGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(materialsField)) CALL FlagError("Materials field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
#endif    

    materialsField=>equationsInterpolation%materialsField

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(materialsField)) &
      & CALL FlagError("Materials field is not assocaited for the equations interpolation.",err,error,*999)
#endif    
       
    EXITS("EquationsInterpolation_MaterialsFieldGet")
    RETURN
999 NULLIFY(materialsField)
998 ERRORS("EquationsInterpolation_MaterialsFieldGet",err,error)
    EXITS("EquationsInterpolation_MaterialsFieldGet")
    RETURN 1
    
  END SUBROUTINE EquationsInterpolation_MaterialsFieldGet

  !
  !================================================================================================================================
  !

  !>Gets the materials interpolation parameters for an equations interpolation.
  SUBROUTINE EquationsInterpolation_MaterialsParametersGet(equationsInterpolation,variableType,materialsParameters,err,error,*)

    !Argument variables
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation !<A pointer to the equations interpolation to get the materials parameters for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to get the materials parameters for. \see FIELD_ROUTINES_VariableTypes
    TYPE(FieldInterpolationParametersType), POINTER :: materialsParameters !<On exit, a pointer to the materials interpolation parameters for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("EquationsInterpolation_MaterialsParametersGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(materialsParameters)) CALL FlagError("Materials parameters is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equationsInterpolation%materialsInterpParameters)) &
      & CALL FlagError("Equations interpolation materials interpolation parameters is not associated.",err,error,*999)
    IF(variableType<1.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    materialsParameters=>equationsInterpolation%materialsInterpParameters(variableType)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(materialsParameters)) THEN
      localError="Equations interpolation materials parameters is not associated for field variable type "// &
        & TRIM(NumberToVString(variableType,"*",err,error))//" of the equations interpolation."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
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
    TYPE(FieldInterpolatedPointType), POINTER :: materialsPoint !<On exit, a pointer to the materials interpolated point for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("EquationsInterpolation_MaterialsPointGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(materialsPoint)) CALL FlagError("Materials point is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equationsInterpolation%materialsInterpPoint)) &
      & CALL FlagError("Equations interpolation materials interpolated point is not associated.",err,error,*999)
    IF(variableType<1.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    materialsPoint=>equationsInterpolation%materialsInterpPoint(variableType)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(materialsPoint)) THEN
      localError="Equations interpolated materials point is not associated for field variable type "// &
        & TRIM(NumberToVString(variableType,"*",err,error))//" of the equations interpolation."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
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
    TYPE(FieldInterpolationParametersType), POINTER :: prevDependentParameters !<On exit, a pointer to the previous dependent interpolation parameters for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("EquationsInterpolation_PreviousDependentParametersGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(prevDependentParameters)) CALL FlagError("Previous dependent parameters is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equationsInterpolation%prevDependentInterpParameters)) &
      & CALL FlagError("Equations interpolation previous dependent interpolation parameters is not associated.",err,error,*999)
    IF(variableType<1.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    prevDependentParameters=>equationsInterpolation%prevDependentInterpParameters(variableType)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(prevDependentParameters)) THEN
      localError="Equations interpolation previous dependent parameters is not associated for field variable type "// &
        & TRIM(NumberToVString(variableType,"*",err,error))//" of the equations interpolation."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
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
    TYPE(FieldInterpolatedPointType), POINTER :: prevDependentPoint !<On exit, a pointer to the previous dependent interpolated point for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("EquationsInterpolation_PreviousDependentPointGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(prevDependentPoint)) CALL FlagError("Previous dependent point is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equationsInterpolation%prevDependentInterpPoint)) &
      & CALL FlagError("Equations interpolation previous dependent interpolated point is not associated.",err,error,*999)
    IF(variableType<1.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    prevDependentPoint=>equationsInterpolation%prevDependentInterpPoint(variableType)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(prevDependentPoint)) THEN
      localError="Equations interpolated previous dependent point is not associated for field variable type "// &
        & TRIM(NumberToVString(variableType,"*",err,error))//" of the equations interpolation."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
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
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: prevDependentPointMetrics !<On exit, a pointer to the previous dependent interpolated point metrics for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("EquationsInterpolation_PreviousDependentPointMetricsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(prevDependentPointMetrics)) CALL FlagError("Previous dependent point metrics is already associated.", &
      & err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equationsInterpolation%prevDependentInterpPointMetrics)) &
      & CALL FlagError("Equations interpolation previous dependent interpolated point metrics is not associated.",err,error,*999)
    IF(variableType<1.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    prevDependentPointMetrics=>equationsInterpolation%prevDependentInterpPointMetrics(variableType)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(prevDependentPointMetrics)) THEN
      localError="Equations interpolated previous dependent point metrics is not associated for field variable type "// &
        & TRIM(NumberToVString(variableType,"*",err,error))//" of the equations interpolation."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
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

  !>Checks the source field for an equations interpolation.
  SUBROUTINE EquationsInterpolation_SourceFieldExists(equationsInterpolation,sourceField,err,error,*)

    !Argument variables
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation !<A pointer to the equations interpolation to get the source field for
    TYPE(FieldType), POINTER :: sourceField !<On exit, a pointer to the source field for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsInterpolation_SourceFieldExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(sourceField)) CALL FlagError("Source field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
#endif    

    sourceField=>equationsInterpolation%sourceField
       
    EXITS("EquationsInterpolation_SourceFieldExists")
    RETURN
999 NULLIFY(sourceField)
998 ERRORS("EquationsInterpolation_SourceFieldExists",err,error)
    EXITS("EquationsInterpolation_SourceFieldExists")
    RETURN 1
    
  END SUBROUTINE EquationsInterpolation_SourceFieldExists

  !
  !================================================================================================================================
  !

  !>Gets the source field for an equations interpolation.
  SUBROUTINE EquationsInterpolation_SourceFieldGet(equationsInterpolation,sourceField,err,error,*)

    !Argument variables
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation !<A pointer to the equations interpolation to get the source field for
    TYPE(FieldType), POINTER :: sourceField !<On exit, a pointer to the source field for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsInterpolation_SourceFieldGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(sourceField)) CALL FlagError("Source field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
#endif    

    sourceField=>equationsInterpolation%sourceField

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(sourceField)) &
      & CALL FlagError("Source field is not assocaited for the equations interpolation.",err,error,*999)
#endif    
       
    EXITS("EquationsInterpolation_SourceFieldGet")
    RETURN
999 NULLIFY(sourceField)
998 ERRORS("EquationsInterpolation_SourceFieldGet",err,error)
    EXITS("EquationsInterpolation_SourceFieldGet")
    RETURN 1
    
  END SUBROUTINE EquationsInterpolation_SourceFieldGet

  !
  !================================================================================================================================
  !

  !>Gets the source interpolation parameters for an equations interpolation.
  SUBROUTINE EquationsInterpolation_SourceParametersGet(equationsInterpolation,variableType,sourceParameters,err,error,*)

    !Argument variables
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation !<A pointer to the equations interpolation to get the source parameters for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to get the source parameters for. \see FIELD_ROUTINES_VariableTypes
    TYPE(FieldInterpolationParametersType), POINTER :: sourceParameters !<On exit, a pointer to the source interpolation parameters for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("EquationsInterpolation_SourceParametersGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(sourceParameters)) CALL FlagError("Source parameters is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equationsInterpolation%sourceInterpParameters)) &
      & CALL FlagError("Equations interpolation source interpolation parameters is not associated.",err,error,*999)
    IF(variableType<1.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    sourceParameters=>equationsInterpolation%sourceInterpParameters(variableType)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(sourceParameters)) THEN
      localError="Equations interpolation source parameters is not associated for field variable type "// &
        & TRIM(NumberToVString(variableType,"*",err,error))//" of the equations interpolation."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
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
    TYPE(FieldInterpolatedPointType), POINTER :: sourcePoint !<On exit, a pointer to the source interpolated point for the specified equations interpolation. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("EquationsInterpolation_SourcePointGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(sourcePoint)) CALL FlagError("Source point is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equationsInterpolation%sourceInterpPoint)) &
      & CALL FlagError("Equations interpolation source interpolated point is not associated.",err,error,*999)
    IF(variableType<1.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type needs to be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    sourcePoint=>equationsInterpolation%sourceInterpPoint(variableType)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(sourcePoint)) THEN
      localError="Equations interpolated source point is not associated for field variable type "// &
        & TRIM(NumberToVString(variableType,"*",err,error))//" of the equations interpolation."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
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

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(equations)) CALL FlagError("Equations is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(scalarEquations)) CALL FlagError("Scalar equations is not associated.",err,error,*999)
#endif    

    equations=>scalarEquations%equations

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated for the scalar equations.",err,error,*999)
#endif    
       
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

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(scalarMapping)) CALL FlagError("Scalar mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(scalarEquations)) CALL FlagError("Scalar equations is not associated.",err,error,*999)
#endif    

    scalarMapping=>scalarEquations%scalarMapping

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(scalarMapping)) CALL FlagError("Scalar mapping is not associated for the scalar equations.",err,error,*999)
#endif    
       
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

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(scalarMatrices)) CALL FlagError("Scalar matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(scalarEquations)) CALL FlagError("Scalar equations is not associated.",err,error,*999)
#endif    

    scalarMatrices=>scalarEquations%scalarMatrices

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(scalarMatrices)) CALL FlagError("Scalar matrices is not associated for the scalar equations.",err,error,*999)
#endif    
       
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

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(equations)) CALL FlagError("Equations is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorEquations)) CALL FlagError("Vector equations is not associated.",err,error,*999)
#endif    

    equations=>vectorEquations%equations

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated for the vector equations.",err,error,*999)
#endif    
       
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

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(vectorMapping)) CALL FlagError("Vector mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorEquations)) CALL FlagError("Vector equations is not associated.",err,error,*999)
#endif    

    vectorMapping=>vectorEquations%vectorMapping

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(vectorMapping)) CALL FlagError("Vector mapping is not associated for the vector equations.",err,error,*999)
#endif    
       
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

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(vectorEquations)) CALL FlagError("Vector equations is not associated.",err,error,*999)
#endif    

    vectorMatrices=>vectorEquations%vectorMatrices

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(vectorMatrices)) CALL FlagError("Vector matrices is not associated for the vector equations.",err,error,*999)
#endif    
       
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
