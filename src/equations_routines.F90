!> \file
!> \author Chris Bradley
!> \brief This module handles all equations routines.
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

!> This module handles all equations routines.
MODULE EquationsRoutines

  USE BaseRoutines
  USE EquationsAccessRoutines
  USE EquationsMappingRoutines
  USE EquationsMappingAccessRoutines
  USE EquationsMatricesRoutines
  USE EquationsMatricesAccessRoutines
  USE EquationsSetAccessRoutines
  USE FieldRoutines
  USE FieldAccessRoutines
  USE ISO_VARYING_STRING
  USE Kinds
  USE Strings
  USE Types

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

 
  !Module types

  !Module variables

  !Interfaces
 
  PUBLIC Equations_Initialise,Equations_Finalise

  PUBLIC Equations_CreateStart,Equations_CreateFinish

  PUBLIC Equations_Destroy

  PUBLIC Equations_EqualityTypeSet
  
  PUBLIC Equations_EquationTypeSet

  PUBLIC Equations_LinearityTypeSet

  PUBLIC Equations_LumpingTypeSet

  PUBLIC Equations_OutputTypeSet

  PUBLIC Equations_SparsityTypeSet

  PUBLIC Equations_TimeDependenceTypeSet

CONTAINS

  !
  !================================================================================================================================
  !

  !>Finish the creation of equations.
  SUBROUTINE Equations_CreateFinish(equations,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER :: equations !<A pointer to the equations to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Equations_CreateFinish",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
    IF(equations%equationsFinished) CALL FlagError("Equations have already been finished.",err,error,*999)
    
    !Initialise the eqations interpolation
    CALL Equations_InterpolationInitialise(equations,err,error,*999)
    !Set the finished flag
    equations%equationsFinished=.TRUE.
       
    EXITS("Equations_CreateFinish")
    RETURN
999 ERRORSEXITS("Equations_CreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE Equations_CreateFinish

  !
  !================================================================================================================================
  !

  !>Start the creation of equations for the equation set.
  SUBROUTINE Equations_CreateStart(equationsSet,equations,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to create equations for
    TYPE(EquationsType), POINTER :: equations !<On exit, a pointer to the created equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Equations_CreateStart",err,error,*998)

    IF(ASSOCIATED(equations)) CALL FlagError("Equations is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*998)
    IF(ASSOCIATED(equationsSet%equations)) CALL FlagError("Equations are already associated for the equations set.",err,error,*998)

    !Initialise the equations
    CALL Equations_Initialise(equationsSet,err,error,*999)
    !Return the pointer
    equations=>equationsSet%equations
       
    EXITS("Equations_CreateStart")
    RETURN
999 NULLIFY(equations)
998 ERRORSEXITS("Equations_CreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE Equations_CreateStart

  !
  !================================================================================================================================
  !

  !>Destroys equations
  SUBROUTINE Equations_Destroy(equations,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER :: equations !<A pointer to the equations to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Equations_Destroy",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
    
    CALL Equations_Finalise(equations,err,error,*999)
       
    EXITS("Equations_Destroy")
    RETURN
999 ERRORSEXITS("Equations_Destroy",err,error)
    RETURN 1
    
  END SUBROUTINE Equations_Destroy

  !
  !================================================================================================================================
  !

  !>Sets/changes the equality type for equations.
  SUBROUTINE Equations_EqualityTypeSet(equations,equalityType,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER :: equations !<A pointer to the equations to set the equality for
    INTEGER(INTG), INTENT(IN) :: equalityType !<The equality type to set \see EquationsRoutines_EquationEqualityTypes,EquationsRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Equations_EqualityTypeSet",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
    IF(equations%equationsFinished) CALL FlagError("Equations has already been finished.",err,error,*999)

    SELECT CASE(equalityType)
    CASE(EQUATIONS_EQUALS_TYPE)
      equations%equalityType=EQUATIONS_EQUALS_TYPE
    CASE(EQUATIONS_LESS_THAN_TYPE)
      equations%equalityType=EQUATIONS_LESS_THAN_TYPE
    CASE(EQUATIONS_LESS_THAN_EQUALS_TYPE)
      equations%equalityType=EQUATIONS_LESS_THAN_EQUALS_TYPE
    CASE(EQUATIONS_GREATER_THAN_TYPE)
      equations%equalityType=EQUATIONS_GREATER_THAN_TYPE
    CASE(EQUATIONS_GREATER_THAN_EQUALS_TYPE)
      equations%equalityType=EQUATIONS_GREATER_THAN_EQUALS_TYPE
    CASE DEFAULT
      localError="The specified equation equality type of "//TRIM(NumberToVString(equalityType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("Equations_EqualityTypeSet")
    RETURN
999 ERRORSEXITS("Equations_EqualityTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE Equations_EqualityTypeSet
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the equation type for equations.
  SUBROUTINE Equations_EquationTypeSet(equations,equationType,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER :: equations !<A pointer to the equations to set the equation for
    INTEGER(INTG), INTENT(IN) :: equationType !<The equation type to set \see EquationsRoutines_EquationsTypes,EquationsRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError,localError
 
    ENTERS("Equations_EquationTypeSet",err,error,*998)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*998)
    IF(equations%equationsFinished) CALL FlagError("Equations has already been finished.",err,error,*998)

    IF(equationType/=equations%equationType) THEN
      !Initialise the new equation type
      SELECT CASE(equationType)
      CASE(EQUATIONS_SCALAR_TYPE)
        CALL Equations_ScalarInitialise(equations,err,error,*999)
      CASE(EQUATIONS_VECTOR_TYPE)
        CALL Equations_VectorInitialise(equations,err,error,*999)
      CASE(EQUATIONS_FUNCTIONAL_TYPE)
       CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The specified equation type of "//TRIM(NumberToVString(equationType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Finalise the old equation type
      SELECT CASE(equations%equationType)
      CASE(EQUATIONS_SCALAR_TYPE)
        CALL Equations_ScalarFinalise(equations%scalarEquations,err,error,*999)
      CASE(EQUATIONS_VECTOR_TYPE)
        CALL Equations_VectorFinalise(equations%vectorEquations,err,error,*999)
      CASE(EQUATIONS_FUNCTIONAL_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The specified equation type of "//TRIM(NumberToVString(equationType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      equations%equationType=equationType
    ENDIF
      
    EXITS("Equations_EquationTypeSet")
    RETURN
999 SELECT CASE(equationType)
    CASE(EQUATIONS_SCALAR_TYPE)
      CALL Equations_ScalarFinalise(equations%scalarEquations,dummyErr,dummyError,*998)
    CASE(EQUATIONS_VECTOR_TYPE)
      CALL Equations_VectorFinalise(equations%vectorEquations,dummyErr,dummyError,*998)
    END SELECT
998 ERRORSEXITS("Equations_EquationTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE Equations_EquationTypeSet
  
  !
  !================================================================================================================================
  !

  !>Finalise the equations and deallocate all memory.
  SUBROUTINE Equations_Finalise(equations,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER :: equations !<A pointer to the equations to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Equations_Finalise",err,error,*999)

    IF(ASSOCIATED(equations)) THEN
      CALL Equations_InterpolationFinalise(equations%interpolation,err,error,*999)
      CALL Equations_ScalarFinalise(equations%scalarEquations,err,error,*999)
      CALL Equations_VectorFinalise(equations%vectorEquations,err,error,*999)
      DEALLOCATE(equations)
    ENDIF
       
    EXITS("Equations_Finalise")
    RETURN
999 ERRORSEXITS("Equations_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE Equations_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the equations for an equations set.
  SUBROUTINE Equations_Initialise(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to initialise the equations for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
 
    ENTERS("Equations_Initialise",err,error,*998)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated",err,error,*998)
    IF(ASSOCIATED(equationsSet%equations)) CALL FlagError("Equations is already associated for this equations set.",err,error,*998)

    ALLOCATE(equationsSet%equations,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate equations.",err,error,*999)
    equationsSet%equations%equationsSet=>equationsSet
    equationsSet%equations%equationsFinished=.FALSE.
    equationsSet%equations%equationType=EQUATIONS_VECTOR_TYPE
    equationsSet%equations%equalityType=EQUATIONS_EQUALS_TYPE
    equationsSet%equations%linearity=EQUATIONS_LINEAR
    equationsSet%equations%timeDependence=EQUATIONS_STATIC
    equationsSet%equations%outputType=EQUATIONS_NO_OUTPUT
    equationsSet%equations%sparsityType=EQUATIONS_SPARSE_MATRICES
    equationsSet%equations%lumpingType=EQUATIONS_UNLUMPED_MATRICES
    NULLIFY(equationsSet%equations%interpolation)
    NULLIFY(equationsSet%equations%scalarEquations)
    NULLIFY(equationsSet%equations%vectorEquations)
    CALL Equations_VectorInitialise(equationsSet%equations,err,error,*999)
       
    EXITS("Equations_Initialise")
    RETURN
999 CALL Equations_Finalise(equationsSet%equations,dummyErr,dummyError,*998)
998 ERRORSEXITS("Equations_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE Equations_Initialise
  
  !
  !================================================================================================================================
  !
  
  !>Finalises the interpolation information for equations and deallocates all memory
  SUBROUTINE Equations_InterpolationFinalise(equationsInterpolation,err,error,*)

    !Argument variables
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation !<A pointer to the equations interpolation to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Equations_InterpolationFinalise",err,error,*999)

    IF(ASSOCIATED(equationsInterpolation)) THEN
      CALL Field_InterpolationParametersFinalise(equationsInterpolation%geometricInterpParameters,err,error,*999)
      CALL Field_InterpolationParametersFinalise(equationsInterpolation%fibreInterpParameters,err,error,*999)
      CALL Field_InterpolationParametersFinalise(equationsInterpolation%dependentInterpParameters,err,error,*999)
      CALL Field_InterpolationParametersFinalise(equationsInterpolation%prevDependentInterpParameters,err,error,*999)
      CALL Field_InterpolationParametersFinalise(equationsInterpolation%independentInterpParameters,err,error,*999)
      CALL Field_InterpolationParametersFinalise(equationsInterpolation%materialsInterpParameters,err,error,*999)
      CALL Field_InterpolationParametersFinalise(equationsInterpolation%sourceInterpParameters,err,error,*999)
      CALL Field_InterpolatedPointsFinalise(equationsInterpolation%geometricInterpPoint,err,error,*999)
      CALL Field_InterpolatedPointsFinalise(equationsInterpolation%dependentInterpPoint,err,error,*999)
      CALL Field_InterpolatedPointsFinalise(equationsInterpolation%prevDependentInterpPoint,err,error,*999)
      CALL Field_InterpolatedPointsFinalise(equationsInterpolation%independentInterpPoint,err,error,*999)
      CALL Field_InterpolatedPointsFinalise(equationsInterpolation%fibreInterpPoint,err,error,*999)
      CALL Field_InterpolatedPointsFinalise(equationsInterpolation%materialsInterpPoint,err,error,*999)
      CALL Field_InterpolatedPointsFinalise(equationsInterpolation%sourceInterpPoint,err,error,*999)
      CALL Field_PhysicalPointsFinalise(equationsInterpolation%dependentPhysicalPoint,err,error,*999)
      CALL Field_InterpolatedPointsMetricsFinalise(equationsInterpolation%dependentInterpPointMetrics,err,error,*999)
      CALL Field_InterpolatedPointsMetricsFinalise(equationsInterpolation%prevDependentInterpPointMetrics,err,error,*999)
      CALL Field_InterpolatedPointsMetricsFinalise(equationsInterpolation%independentInterpPointMetrics,err,error,*999)
      CALL Field_InterpolatedPointsMetricsFinalise(equationsInterpolation%geometricInterpPointMetrics,err,error,*999)
      CALL Field_InterpolatedPointsMetricsFinalise(equationsInterpolation%fibreInterpPointMetrics,err,error,*999)
      DEALLOCATE(equationsInterpolation)
    ENDIF
       
    EXITS("Equations_InterpolationFinalise")
    RETURN
999 ERRORSEXITS("Equations_InterpolationFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE Equations_InterpolationFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the interpolation information for equations
  SUBROUTINE Equations_InterpolationInitialise(equations,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER :: equations !<The pointer to the equations to initialise the interpolation for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("Equations_InterpolationInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated",err,error,*998)
    IF(ASSOCIATED(equations%interpolation)) &
      & CALL FlagError("Interpolation is already associated for these equations.",err,error,*998)
    NULLIFY(equationsSet)
    CALL Equations_EquationsSetGet(equations,equationsSet,err,error,*999)    

    ALLOCATE(equations%interpolation,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate equations interpolation",err,error,*999)
    equations%interpolation%equations=>equations
    NULLIFY(equations%interpolation%geometricInterpParameters)
    NULLIFY(equations%interpolation%fibreInterpParameters)
    NULLIFY(equations%interpolation%dependentInterpParameters)
    NULLIFY(equations%interpolation%prevDependentInterpParameters)
    NULLIFY(equations%interpolation%independentInterpParameters)
    NULLIFY(equations%interpolation%materialsInterpParameters)
    NULLIFY(equations%interpolation%sourceInterpParameters)
    NULLIFY(equations%interpolation%geometricInterpPoint)
    NULLIFY(equations%interpolation%fibreInterpPoint)
    NULLIFY(equations%interpolation%dependentInterpPoint)
    NULLIFY(equations%interpolation%prevDependentInterpPoint)
    NULLIFY(equations%interpolation%independentInterpPoint)
    NULLIFY(equations%interpolation%materialsInterpPoint)
    NULLIFY(equations%interpolation%sourceInterpPoint)
    NULLIFY(equations%interpolation%dependentPhysicalPoint)
    NULLIFY(equations%interpolation%dependentInterpPointMetrics)
    NULLIFY(equations%interpolation%prevDependentInterpPointMetrics)
    NULLIFY(equations%interpolation%independentInterpPointMetrics)
    NULLIFY(equations%interpolation%geometricInterpPointMetrics)
    NULLIFY(equations%interpolation%fibreInterpPointMetrics)
    
    equations%interpolation%geometricField=>equationsSet%geometry%geometricField
    equations%interpolation%fibreField=>equationsSet%geometry%fibreField
    equations%interpolation%dependentField=>equationsSet%dependent%dependentField
    IF(ASSOCIATED(equationsSet%independent)) THEN
      equations%interpolation%independentField=>equationsSet%independent%independentField
    ELSE
      NULLIFY(equations%interpolation%independentField)
    ENDIF
    IF(ASSOCIATED(equationsSet%materials)) THEN
      equations%interpolation%materialsField=>equationsSet%materials%materialsField
    ELSE
      NULLIFY(equations%interpolation%materialsField)
    ENDIF
    IF(ASSOCIATED(equationsSet%source)) THEN
      equations%interpolation%sourceField=>equationsSet%source%sourceField
    ELSE
      NULLIFY(equations%interpolation%sourceField)
    ENDIF
    
    CALL Field_InterpolationParametersInitialise(equations%interpolation%geometricField, &
      & equations%interpolation%geometricInterpParameters,err,error,*999)
    CALL Field_InterpolatedPointsInitialise(equations%interpolation%geometricInterpParameters, &
      & equations%interpolation%geometricInterpPoint,err,error,*999)
    CALL Field_InterpolatedPointsMetricsInitialise(equations%interpolation%geometricInterpPoint, &
      & equations%interpolation%geometricInterpPointMetrics,err,error,*999)
    CALL Field_InterpolationParametersInitialise(equations%interpolation%dependentField, &
      & equations%interpolation%dependentInterpParameters,err,error,*999)
    CALL Field_InterpolatedPointsInitialise(equations%interpolation%dependentInterpParameters, &
      & equations%interpolation%dependentInterpPoint,err,error,*999)
    !CALL Field_PhysicalPointsInitialise(equations%interpolation%dependentInterpPoint, &
    !  & equations%interpolation%geometricInterpPoint,equations%interpolation%dependentPhysicalPoint, &
    !  & err,error,*999)
    IF(equations%interpolation%dependentField%type==FIELD_GEOMETRIC_TYPE.OR. &
      & equations%interpolation%dependentField%type==FIELD_FIBRE_TYPE.OR. &
      & equations%interpolation%dependentField%type==FIELD_GEOMETRIC_GENERAL_TYPE) THEN
      CALL Field_InterpolatedPointsMetricsInitialise(equations%interpolation%dependentInterpPoint, &
        & equations%interpolation%dependentInterpPointMetrics,err,error,*999)
    ENDIF
    IF(equations%timeDependence/=EQUATIONS_STATIC) THEN
      CALL Field_InterpolationParametersInitialise(equations%interpolation%dependentField, &
        & equations%interpolation%prevDependentInterpParameters,err,error,*999)
      CALL Field_InterpolatedPointsInitialise(equations%interpolation%prevDependentInterpParameters, &
        & equations%interpolation%prevDependentInterpPoint,err,error,*999)
      IF(equations%interpolation%dependentField%type==FIELD_GEOMETRIC_TYPE.OR. &
        & equations%interpolation%dependentField%type==FIELD_FIBRE_TYPE.OR. &
        & equations%interpolation%dependentField%type==FIELD_GEOMETRIC_GENERAL_TYPE) THEN
        CALL Field_InterpolatedPointsMetricsInitialise(equations%interpolation%prevDependentInterpPoint, &
          & equations%interpolation%prevDependentInterpPointMetrics,err,error,*999)
      ENDIF
    ENDIF
    IF(ASSOCIATED(equations%interpolation%fibreField)) THEN
      CALL Field_InterpolationParametersInitialise(equations%interpolation%fibreField, &
        & equations%interpolation%fibreInterpParameters,err,error,*999)
      CALL Field_InterpolatedPointsInitialise(equations%interpolation%fibreInterpParameters,  &
        & equations%interpolation%fibreInterpPoint,err,error,*999)
      CALL Field_InterpolatedPointsMetricsInitialise(equations%interpolation%fibreInterpPoint,  &
        & equations%interpolation%fibreInterpPointMetrics,err,error,*999)
    ENDIF
    IF(ASSOCIATED(equations%interpolation%independentField)) THEN
      CALL Field_InterpolationParametersInitialise(equations%interpolation%independentField, &
        & equations%interpolation%independentInterpParameters,err,error,*999)
      CALL Field_InterpolatedPointsInitialise(equations%interpolation%independentInterpParameters,  &
        & equations%interpolation%independentInterpPoint,err,error,*999)
      IF(equations%interpolation%independentField%TYPE==FIELD_GEOMETRIC_TYPE.OR. &
        & equations%interpolation%independentField%type==FIELD_GEOMETRIC_GENERAL_TYPE.OR. &
        & equations%interpolation%independentField%type==FIELD_FIBRE_TYPE) THEN
        CALL Field_InterpolatedPointsMetricsInitialise(equations%interpolation%independentInterpPoint,  &
          &  equations%interpolation%independentInterpPointMetrics,err,error,*999)
      END IF
    ENDIF
    IF(ASSOCIATED(equations%interpolation%materialsField)) THEN
      CALL Field_InterpolationParametersInitialise(equations%interpolation%materialsField, &
        & equations%interpolation%materialsInterpParameters,err,error,*999)
      CALL Field_InterpolatedPointsInitialise(equations%interpolation%materialsInterpParameters,  &
        & equations%interpolation%materialsInterpPoint,err,error,*999)
    ENDIF
    IF(ASSOCIATED(equations%interpolation%sourceField)) THEN
      CALL Field_InterpolationParametersInitialise(equations%interpolation%sourceField, &
        & equations%interpolation%sourceInterpParameters,err,error,*999)
      CALL Field_InterpolatedPointsInitialise(equations%interpolation%sourceInterpParameters, &
        & equations%interpolation%sourceInterpPoint,err,error,*999)
    ENDIF
       
    EXITS("Equations_InterpolationInitialise")
    RETURN
999 CALL Equations_InterpolationFinalise(equations%interpolation,dummyErr,dummyError,*998)
998 ERRORSEXITS("Equations_InterpolationInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE Equations_InterpolationInitialise

  !
  !================================================================================================================================
  !

  !>Sets/changes the linearity type for equations.
  SUBROUTINE Equations_LinearityTypeSet(equations,linearityType,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER :: equations !<A pointer to the equations to set the linearity for
    INTEGER(INTG), INTENT(IN) :: linearityType !<The linearity type to set \see EquationsSetRoutines_LinearityTypes,EquationsSetRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Equations_LinearityTypeSet",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
    IF(equations%equationsFinished) CALL FlagError("Equations has already been finished.",err,error,*999)

    SELECT CASE(linearityType)
    CASE(EQUATIONS_LINEAR)
      equations%linearity=EQUATIONS_LINEAR
    CASE(EQUATIONS_NONLINEAR)
      equations%linearity=EQUATIONS_NONLINEAR
    CASE(EQUATIONS_NONLINEAR_BCS)
      equations%linearity=EQUATIONS_NONLINEAR_BCS
    CASE DEFAULT
      localError="The specified linearity type of "//TRIM(NumberToVString(linearityType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("Equations_LinearityTypeSet")
    RETURN
999 ERRORSEXITS("Equations_LinearityTypeSet",err,error)
    RETURN 1
  END SUBROUTINE Equations_LinearityTypeSet
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the matrix lumping for the equations.
  SUBROUTINE Equations_LumpingTypeSet(equations,lumpingType,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER :: equations !<A pointer to the equations to set the lumping for
    INTEGER(INTG), INTENT(IN) :: lumpingType !<The lumping type to set \see EquationsRoutines_LumpingTypes,EquationsRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Equations_LumpingTypeSet",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
    IF(equations%equationsFinished) CALL FlagError("Equations has already been finished.",err,error,*999)

    IF(equations%timeDependence==EQUATIONS_FIRST_ORDER_DYNAMIC.OR. &
      & equations%timeDependence==EQUATIONS_SECOND_ORDER_DYNAMIC) THEN
      SELECT CASE(lumpingType)
      CASE(EQUATIONS_UNLUMPED_MATRICES)
        equations%lumpingType=EQUATIONS_UNLUMPED_MATRICES
      CASE(EQUATIONS_LUMPED_MATRICES)
        equations%lumpingType=EQUATIONS_LUMPED_MATRICES
      CASE DEFAULT
        localError="The specified lumping type of "//TRIM(NumberToVString(lumpingType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      localError="Invalid equations time dependence. The equations time dependence of "// &
        & TRIM(NumberToVString(equations%timeDependence,"*",err,error))// &
        & " does not correspond to dynamic equations. You can only set lumping for dynamic equations."
      CALL FlagError(localError,err,error,*999)
    ENDIF
        
    EXITS("Equations_LumpingTypeSet")
    RETURN
999 ERRORSEXITS("Equations_LumpingTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE Equations_LumpingTypeSet
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the output type for the equations.
  SUBROUTINE Equations_OutputTypeSet(equations,outputType,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER :: equations !<A pointer to the equations to set the output type for
    INTEGER(INTG), INTENT(IN) :: outputType !<The output type to set \see EquationsRoutines_OutputTypes,EquationsRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Equations_OutputTypeSet",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
    IF(equations%equationsFinished) CALL FlagError("Equations has already been finished.",err,error,*999)

    SELECT CASE(outputType)
    CASE(EQUATIONS_NO_OUTPUT)
      equations%outputType=EQUATIONS_NO_OUTPUT
    CASE(EQUATIONS_TIMING_OUTPUT)
      equations%outputType=EQUATIONS_TIMING_OUTPUT
    CASE(EQUATIONS_MATRIX_OUTPUT)
      equations%outputType=EQUATIONS_MATRIX_OUTPUT
    CASE(EQUATIONS_ELEMENT_MATRIX_OUTPUT)
      equations%outputType=EQUATIONS_ELEMENT_MATRIX_OUTPUT
    CASE(EQUATIONS_NODAL_MATRIX_OUTPUT)
      equations%outputType=EQUATIONS_NODAL_MATRIX_OUTPUT
    CASE DEFAULT
      localError="The specified output type of "//TRIM(NumberToVString(outputType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("Equations_OutputTypeSet")
    RETURN
999 ERRORSEXITS("Equations_OutputTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE Equations_OutputTypeSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the sparsity type for the equations.
  SUBROUTINE Equations_SparsityTypeSet(equations,sparsityType,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER :: equations !<A pointer to the equations to set the sparsity type for
    INTEGER(INTG), INTENT(IN) :: sparsityType !<The sparsity type to set \see EquationsRoutines_SparsityTypes,EquationsRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Equations_SparsityTypeSet",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
    IF(equations%equationsFinished) CALL FlagError("Equations has already been finished.",err,error,*999)
    
    SELECT CASE(sparsityType)
    CASE(EQUATIONS_SPARSE_MATRICES)
      equations%sparsityType=EQUATIONS_SPARSE_MATRICES
    CASE(EQUATIONS_FULL_MATRICES)
      equations%sparsityType=EQUATIONS_FULL_MATRICES
    CASE DEFAULT
      localError="The specified sparsity type of "//TRIM(NumberToVString(sparsityType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("Equations_SparsityTypeSet")
    RETURN
999 ERRORSEXITS("Equations_SparsityTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE Equations_SparsityTypeSet
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the time dependence type for equations.
  SUBROUTINE Equations_TimeDependenceTypeSet(equations,timeDependenceType,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER :: equations !<A pointer to the equations to set the linearity for
    INTEGER(INTG), INTENT(IN) :: timeDependenceType !<The time dependence type to set \see EquationsRoutines_TimeDependenceTypes,EquationsRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Equations_TimeDependenceTypeSet",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
    IF(equations%equationsFinished) CALL FlagError("Equations has already been finished.",err,error,*999)

    SELECT CASE(timeDependenceType)
    CASE(EQUATIONS_STATIC)
      equations%timeDependence=EQUATIONS_STATIC
    CASE(EQUATIONS_QUASISTATIC)
      equations%timeDependence=EQUATIONS_QUASISTATIC
    CASE(EQUATIONS_FIRST_ORDER_DYNAMIC)
      equations%timeDependence=EQUATIONS_FIRST_ORDER_DYNAMIC
    CASE(EQUATIONS_SECOND_ORDER_DYNAMIC)
      equations%timeDependence=EQUATIONS_SECOND_ORDER_DYNAMIC
    CASE DEFAULT
      localError="The specified time dependence type of "//TRIM(NumberToVString(timeDependenceType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("Equations_TimeDependenceTypeSet")
    RETURN
999 ERRORSEXITS("Equations_TimeDependenceTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE Equations_TimeDependenceTypeSet

   !
   !================================================================================================================================
   !

   !>Finalise the scalar equations and deallocate all memory.
   SUBROUTINE Equations_ScalarFinalise(scalarEquations,err,error,*)

     !Argument variables
     TYPE(EquationsScalarType), POINTER :: scalarEquations !<A pointer to the scalar equations to finalise
     INTEGER(INTG), INTENT(OUT) :: err !<The error code
     TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
     !Local Variables

     ENTERS("Equations_ScalarFinalise",err,error,*999)

     IF(ASSOCIATED(scalarEquations)) THEN
       IF(ASSOCIATED(scalarEquations%scalarMapping)) &
         & CALL EquationsMapping_ScalarDestroy(scalarEquations%scalarMapping,err,error,*999)
       IF(ASSOCIATED(scalarEquations%scalarMatrices)) &
         & CALL EquationsMatrices_ScalarDestroy(scalarEquations%scalarMatrices,err,error,*999)
       DEALLOCATE(scalarEquations)
     ENDIF
       
     EXITS("Equations_ScalarFinalise")
     RETURN
 999 ERRORSEXITS("Equations_ScalarFinalise",err,error)
     RETURN 1
    
   END SUBROUTINE Equations_ScalarFinalise

   !
   !================================================================================================================================
   !

   !>Initialise the scalar equations for an equations
   SUBROUTINE Equations_ScalarInitialise(equations,err,error,*)

     !Argument variables
     TYPE(EquationsType), POINTER :: equations !<A pointer to the equations to initialise the scalar equations for
     INTEGER(INTG), INTENT(OUT) :: err !<The error code
     TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
     !Local Variables
     INTEGER(INTG) :: dummyErr
     TYPE(VARYING_STRING) :: dummyError

     ENTERS("Equations_ScalarInitialise",err,error,*998)

     IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*998)
     IF(ASSOCIATED(equations%scalarEquations)) CALL FlagError("Equations scalar equations are already associated.",err,error,*998)

     ALLOCATE(equations%scalarEquations,STAT=err)
     IF(err/=0) CALL FlagError("Could not allocate equations scalar equations.",err,error,*999)
     equations%scalarEquations%equations=>equations
     NULLIFY(equations%scalarEquations%scalarMapping)
     NULLIFY(equations%scalarEquations%scalarMatrices)
        
     EXITS("Equations_ScalarInitialise")
     RETURN
 999 CALL Equations_ScalarFinalise(equations%scalarEquations,dummyErr,dummyError,*998)
 998 ERRORSEXITS("Equations_ScalarInitialise",err,error)
     RETURN 1
    
   END SUBROUTINE Equations_ScalarInitialise

  !
  !================================================================================================================================
  !

  !>Finalise the vector equations and deallocate all memory.
  SUBROUTINE Equations_VectorFinalise(vectorEquations,err,error,*)

    !Argument variables
    TYPE(EquationsVectorType), POINTER :: vectorEquations !<A pointer to the vector equations to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Equations_VectorFinalise",err,error,*999)

    IF(ASSOCIATED(vectorEquations)) THEN
      IF(ASSOCIATED(vectorEquations%vectorMapping)) &
        & CALL EquationsMapping_VectorDestroy(vectorEquations%vectorMapping,err,error,*999)
      IF(ASSOCIATED(vectorEquations%vectorMatrices)) &
        & CALL EquationsMatrices_VectorDestroy(vectorEquations%vectorMatrices,err,error,*999)
      DEALLOCATE(vectorEquations)
    ENDIF
       
    EXITS("Equations_VectorFinalise")
    RETURN
999 ERRORSEXITS("Equations_VectorFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE Equations_VectorFinalise

  !
  !================================================================================================================================
  !

  !>Initialise the vector equations for an equations
  SUBROUTINE Equations_VectorInitialise(equations,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER :: equations !<A pointer to the equations to initialise the vector equations for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("Equations_VectorInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*998)
    IF(ASSOCIATED(equations%vectorEquations)) CALL FlagError("Equations vector equations are already associated.",err,error,*998)

    ALLOCATE(equations%vectorEquations,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate equations vector equations.",err,error,*999)
    equations%vectorEquations%equations=>equations
    NULLIFY(equations%vectorEquations%vectorMapping)
    NULLIFY(equations%vectorEquations%vectorMatrices)
        
    EXITS("Equations_VectorInitialise")
    RETURN
999 CALL Equations_VectorFinalise(equations%vectorEquations,dummyErr,dummyError,*998)
998 ERRORSEXITS("Equations_VectorInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE Equations_VectorInitialise

  !
  !================================================================================================================================
  !

END MODULE EquationsRoutines
!
