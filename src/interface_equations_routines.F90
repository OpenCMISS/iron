!> \file
!> \author Chris Bradley
!> \brief This module handles all interface equations routines.
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

!>This module handles all interface equations routines.
MODULE InterfaceEquationsRoutines

  USE BaseRoutines
  USE EquationsRoutines
  USE FieldRoutines
  USE FieldAccessRoutines
  USE INTERFACE_CONDITIONS_CONSTANTS
  USE InterfaceConditionAccessRoutines
  USE InterfaceEquationsAccessRoutines
  USE InterfaceMappingRoutines
  USE InterfaceMatricesRoutines
  USE InterfaceMatricesAccessRoutines
  USE INTERFACE_MATRICES_CONSTANTS
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

  !>Gets the time dependence type for an interface matrix
  INTERFACE InterfaceEquations_MatrixTimeDependenceTypeGet
    MODULE PROCEDURE InterfaceEquations_MatrixTimeDependenceTypeGet0
    MODULE PROCEDURE InterfaceEquations_MatrixTimeDependenceTypeGet1
  END INTERFACE InterfaceEquations_MatrixTimeDependenceTypeGet

  !>Sets the time dependence type for an interface matrix
  INTERFACE InterfaceEquations_MatrixTimeDependenceTypeSet
    MODULE PROCEDURE InterfaceEquations_MatrixTimeDependenceTypeSet0
    MODULE PROCEDURE InterfaceEquations_MatrixTimeDependenceTypeSet1
  END INTERFACE InterfaceEquations_MatrixTimeDependenceTypeSet

  PUBLIC InterfaceEquations_CreateFinish,InterfaceEquations_CreateStart

  PUBLIC InterfaceEquations_Destroy

  PUBLIC InterfaceEquations_InterpolationSetsNumberSet

  PUBLIC InterfaceEquations_MatrixTimeDependenceTypeSet,InterfaceEquations_MatrixTimeDependenceTypeGet
  
  PUBLIC InterfaceEquations_OutputTypeGet,InterfaceEquations_OutputTypeSet
  
  PUBLIC InterfaceEquations_SparsityTypeGet,InterfaceEquations_SparsityTypeSet

  PUBLIC InterfaceEquations_LinearityTypeGet,InterfaceEquations_LinearityTypeSet

  PUBLIC InterfaceEquations_TimeDependenceTypeGet,InterfaceEquations_TimeDependenceTypeSet

CONTAINS

  !
  !================================================================================================================================
  !

  !>Finish the creation of interface equations.
  SUBROUTINE InterfaceEquations_CreateFinish(interfaceEquations,err,error,*)

    !Argument variables
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations !<A pointer to the interface equations to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: variableIdx
    TYPE(FieldType), POINTER :: dependentField,geometricField,lagrangeField,nullField,penaltyField
    TYPE(FieldVariableType), POINTER :: dependentVariable
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(InterfaceDependentType), POINTER :: interfaceDependent
    TYPE(InterfaceEquationsInterpolationType), POINTER :: interfaceEquationsInterpolation
    TYPE(VARYING_STRING) :: localError

    ENTERS("InterfaceEquations_CreateFinish",err,error,*999)

    CALL InterfaceEquations_AssertNotFinished(interfaceEquations,err,error,*999)

    !Create the interpolation sets
    NULLIFY(nullField)
    NULLIFY(interfaceCondition)
    CALL InterfaceEquations_InterfaceConditionGet(interfaceEquations,interfaceCondition,err,error,*999)
    SELECT CASE(interfaceCondition%method)
    CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
      NULLIFY(interfaceDependent)
      CALL InterfaceCondition_InterfaceDependentGet(interfaceCondition,interfaceDependent,err,error,*999)
      NULLIFY(interfaceEquationsInterpolation)
      CALL InterfaceEquations_EquationsInterpolationGet(interfaceEquations,interfaceEquationsInterpolation,err,error,*999)
      NULLIFY(geometricField)
      CALL InterfaceCondition_GeometricFieldGet(interfaceCondition,geometricField,err,error,*999)
      NULLIFY(lagrangeField)
      CALL InterfaceCondition_LagrangeFieldGet(interfaceCondition,lagrangeField,err,error,*999)
      NULLIFY(penaltyField)
      IF(ASSOCIATED(interfaceCondition%penalty)) &
        & CALL InterfaceCondition_PenaltyFieldGet(interfaceCondition,penaltyField,err,error,*999)
      CALL InterfaceEquations_DomainInterpolationSetup(interfaceEquationsInterpolation%interfaceInterpolation, &
        & geometricField,lagrangeField,penaltyField,err,error,*999)
      DO variableIdx=1,interfaceDependent%numberOfDependentVariables
        NULLIFY(dependentVariable)
        CALL InterfaceDependent_DependentVariableGet(interfaceDependent,variableIdx,dependentVariable,err,error,*999)
        NULLIFY(dependentField)
        CALL FieldVariable_FieldGet(dependentVariable,dependentField,err,error,*999)
        NULLIFY(geometricField)
        CALL Field_GeometricFieldGet(dependentField,geometricField,err,error,*999)
        CALL InterfaceEquations_DomainInterpolationSetup(interfaceEquationsInterpolation% &
          & variableInterpolation(variableIdx),geometricField,dependentField,nullField,err,error,*999)
      ENDDO !variableIdx
    CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The interface condition method of "//TRIM(NumberToVString(interfaceCondition%method,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    !Set the finished flag
    interfaceEquations%interfaceEquationsFinished=.TRUE.
       
    EXITS("InterfaceEquations_CreateFinish")
    RETURN
999 ERRORSEXITS("InterfaceEquations_CreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_CreateFinish

  !
  !================================================================================================================================
  !

  !>Start the creation of interface equations for an interface condition.
  SUBROUTINE InterfaceEquations_CreateStart(interfaceCondition,interfaceEquations,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition to create interface equations for
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations !<On exit, a pointer to the created interface equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("InterfaceEquations_CreateStart",err,error,*998)

    IF(ASSOCIATED(interfaceEquations)) CALL FlagError("Interface equations is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)
    IF(ASSOCIATED(interfaceCondition%interfaceEquations)) &
      & CALL FlagError("Interface equations are already associated for the interface condition.",err,error,*999)
    
    !Initialise the equations
    CALL InterfaceEquations_Initialise(interfaceCondition,err,error,*999)
    !Return the pointer
    interfaceEquations=>interfaceCondition%interfaceEquations
       
    EXITS("InterfaceEquations_CreateStart")
    RETURN
999 NULLIFY(interfaceEquations)
998 ERRORSEXITS("InterfaceEquations_CreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_CreateStart

  !
  !================================================================================================================================
  !

  !>Destroys interface equations.
  SUBROUTINE InterfaceEquations_Destroy(interfaceEquations,err,error,*)

    !Argument variables
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations !<A pointer to the interface equations to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("InterfaceEquations_Destroy",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceEquations)) CALL FlagError("Interface equations is not associated.",err,error,*999)
    
    CALL InterfaceEquations_Finalise(interfaceEquations,err,error,*999)
       
    EXITS("InterfaceEquations_Destroy")
    RETURN
999 ERRORSEXITS("InterfaceEquations_Destroy",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_Destroy

  !
  !================================================================================================================================
  !

  !>Finalises the interface equations domain interpolation and deallocates all memory.
  SUBROUTINE InterfaceEquations_DomainInterpolationFinalise(domainInterpolation,err,error,*)

    !Argument variables
    TYPE(InterfaceEquationsDomainInterpolationType) :: domainInterpolation !<The domain interpolation to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: interpolationSetIdx
 
    ENTERS("InterfaceEquations_DomainInterpolationFinalise",err,error,*999)

    NULLIFY(domainInterpolation%geometricField)
    IF(ALLOCATED(domainInterpolation%geometricInterpolation)) THEN
      DO interpolationSetIdx=1,SIZE(domainInterpolation%geometricInterpolation,1)
        CALL InterfaceEquations_InterpolationSetFinalise(domainInterpolation%geometricInterpolation(interpolationSetIdx), &
          & err,error,*999)
      ENDDO !interpolationSetIdx
      DEALLOCATE(domainInterpolation%geometricInterpolation)
    ENDIF
    domainInterpolation%numberOfGeometricInterpolationSets=0
    NULLIFY(domainInterpolation%dependentField)
    IF(ALLOCATED(domainInterpolation%dependentInterpolation)) THEN
     DO interpolationSetIdx=1,SIZE(domainInterpolation%dependentInterpolation,1)
        CALL InterfaceEquations_InterpolationSetFinalise(domainInterpolation%dependentInterpolation(interpolationSetIdx), &
          & err,error,*999)
      ENDDO !interpolationSetIdx
      DEALLOCATE(domainInterpolation%dependentInterpolation)
    ENDIF
    domainInterpolation%numberOfDependentInterpolationSets=0
       
    EXITS("InterfaceEquations_DomainInterpolationFinalise")
    RETURN
999 ERRORS("InterfaceEquations_DomainInterpolationFinalise",err,error)
    EXITS("InterfaceEquations_DomainInterpolationFinalise")
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_DomainInterpolationFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the interface equations domain interpolation.
  SUBROUTINE InterfaceEquations_DomainInterpolationInitialise(domainInterpolation,err,error,*)

    !Argument variables
    TYPE(InterfaceEquationsDomainInterpolationType) :: domainInterpolation !<The domain interpolation to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceEquations_DomainInterpolationInitialise",err,error,*999)

    NULLIFY(domainInterpolation%penaltyField)
    domainInterpolation%numberOfPenaltyInterpolationSets=1
    NULLIFY(domainInterpolation%geometricField)
    domainInterpolation%numberOfGeometricInterpolationSets=1
    NULLIFY(domainInterpolation%dependentField)
    domainInterpolation%numberOfDependentInterpolationSets=1
       
    EXITS("InterfaceEquations_DomainInterpolationInitialise")
    RETURN
999 ERRORS("InterfaceEquations_DomainInterpolationInitialise",err,error)
    EXITS("InterfaceEquations_DomainInterpolationInitialise")
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_DomainInterpolationInitialise

  !
  !================================================================================================================================
  !

  !>Sets up the interface equations domain interpolation. 
  SUBROUTINE InterfaceEquations_DomainInterpolationSetup(domainInterpolation,geometricField,lagrangeField, &
    & penaltyField,err,error,*)

    !Argument variables
    TYPE(InterfaceEquationsDomainInterpolationType) :: domainInterpolation !<The domain interpolation to initialise
    TYPE(FieldType), POINTER :: geometricField !<A pointer to the geometric field to set up the domain interpolation for
    TYPE(FieldType), POINTER :: lagrangeField !<A pointer to the Lagrange field to set up the domain interpoaltion for
    TYPE(FieldType), POINTER :: penaltyField !<A pointer to the penalty field to set up the domain interpoaltion for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,interpolationSetIdx
    TYPE(VARYING_STRING) :: dummyError
 
    ENTERS("InterfaceEquations_DomainInterpolationSetup",err,error,*998)

    IF(.NOT.ASSOCIATED(geometricField)) CALL FlagError("Geometric field is not associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(lagrangeField)) CALL FlagError("Lagrange field is not associated.",err,error,*998)
    
    domainInterpolation%geometricField=>geometricField
    ALLOCATE(domainInterpolation%geometricInterpolation(domainInterpolation%numberOfGeometricInterpolationSets),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate domain interpolation geometric interpolation.",err,error,*999)
    DO interpolationSetIdx=1,domainInterpolation%numberOfGeometricInterpolationSets
      CALL InterfaceEquations_InterpolationSetInitialise(domainInterpolation%geometricInterpolation(interpolationSetIdx), &
        & err,error,*999)
      CALL Field_InterpolationParametersInitialise(domainInterpolation%geometricField,domainInterpolation% &
        & geometricInterpolation(interpolationSetIdx)%interpolationParameters,err,error,*999)
      CALL Field_InterpolatedPointsInitialise(domainInterpolation%geometricInterpolation(interpolationSetIdx)% &
        & interpolationParameters,domainInterpolation%geometricInterpolation(interpolationSetIdx)%interpolatedPoint, &
        & err,error,*999)
      IF(geometricField%TYPE==FIELD_GEOMETRIC_TYPE.OR.geometricField%TYPE==FIELD_FIBRE_TYPE) &
        & CALL Field_InterpolatedPointsMetricsInitialise(domainInterpolation%geometricInterpolation(interpolationSetIdx)% &
        & interpolatedPoint,domainInterpolation%geometricInterpolation(interpolationSetIdx)%interpolatedPointMetrics, &
        & err,error,*999)
    ENDDO !interpolationSetIdx
    domainInterpolation%dependentField=>lagrangeField
    ALLOCATE(domainInterpolation%dependentInterpolation(domainInterpolation%numberOfDependentInterpolationSets),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate domain interpolation dependent interpolation.",err,error,*999)
    DO interpolationSetIdx=1,domainInterpolation%numberOfDependentInterpolationSets
      CALL InterfaceEquations_InterpolationSetInitialise(domainInterpolation%dependentInterpolation(interpolationSetIdx), &
        & err,error,*999)
      CALL Field_InterpolationParametersInitialise(domainInterpolation%dependentField,domainInterpolation% &
        & dependentInterpolation(interpolationSetIdx)%interpolationParameters,err,error,*999)
      CALL Field_InterpolatedPointsInitialise(domainInterpolation%dependentInterpolation(interpolationSetIdx)% &
        & interpolationParameters,domainInterpolation%dependentInterpolation(interpolationSetIdx)%interpolatedPoint, &
        & err,error,*999)
    ENDDO !interpolationSetIdx
    IF(ASSOCIATED(penaltyField)) THEN
      domainInterpolation%penaltyField=>penaltyField
      ALLOCATE(domainInterpolation%penaltyInterpolation(domainInterpolation%numberOfPenaltyInterpolationSets),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate domain interpolation dependent interpolation.",err,error,*999)
      DO interpolationSetIdx=1,domainInterpolation%numberOfPenaltyInterpolationSets
        CALL InterfaceEquations_InterpolationSetInitialise(domainInterpolation%penaltyInterpolation(interpolationSetIdx), &
          & err,error,*999)
        CALL Field_InterpolationParametersInitialise(domainInterpolation%penaltyField,domainInterpolation% &
          & penaltyInterpolation(interpolationSetIdx)%interpolationParameters,err,error,*999)
        CALL Field_InterpolatedPointsInitialise(domainInterpolation%penaltyInterpolation(interpolationSetIdx)% &
          & interpolationParameters,domainInterpolation%penaltyInterpolation(interpolationSetIdx)%interpolatedPoint, &
          & err,error,*999)
      ENDDO !interpolationSetIdx
    ENDIF
    
    EXITS("InterfaceEquations_DomainInterpolationSetup")
    RETURN
999 CALL InterfaceEquations_DomainInterpolationFinalise(domainInterpolation,dummyErr,dummyError,*998)
998 ERRORS("InterfaceEquations_DomainInterpolationSetup",err,error)
    EXITS("InterfaceEquations_DomainInterpolationSetup")
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_DomainInterpolationSetup

  !
  !================================================================================================================================
  !

  SUBROUTINE InterfaceEquations_InterpolationSetsNumberSet(domainInterpolation,numberOfGeometricSets, &
     & numberOfDependentSets,numberOfPenaltySets,err,error,*)

    !Argument variables
    TYPE(InterfaceEquationsDomainInterpolationType) :: domainInterpolation !<The domain interpolation to set the number of interpolation sets for
    INTEGER(INTG), INTENT(IN) :: numberOfGeometricSets !<The number of geometric interface interpolation sets to set
    INTEGER(INTG), INTENT(IN) :: numberOfDependentSets !<The number of dependent interface interpolation sets to set
    INTEGER(INTG), INTENT(IN) :: numberOfPenaltySets !<The number of penalty interface interpolation sets to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("InterfaceEquations_InterpolationSetsNumberSet",err,error,*999)
    
    IF(numberOfGeometricSets<=0) THEN
      localError="The specified number of geometric sets of "//TRIM(NumberToVString(numberOfGeometricSets,"*",err,error))// &
        & " is invalid. The number of geometric sets must be >= 1."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(numberOfDependentSets<=0) THEN
      localError="The specified number of dependent sets of "//TRIM(NumberToVString(numberOfDependentSets,"*",err,error))// &
        & " is invalid. The number of dependent sets must be >= 1."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(numberOfPenaltySets<0) THEN
      localError="The specified number of penalty sets of "//TRIM(NumberToVString(numberOfPenaltySets,"*",err,error))// &
        & " is invalid. The number of penalty sets must be >= 0."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    domainInterpolation%numberOfGeometricInterpolationSets=numberOfGeometricSets
    domainInterpolation%numberOfDependentInterpolationSets=numberOfDependentSets
    domainInterpolation%numberOfPenaltyInterpolationSets=numberOfPenaltySets
       
    EXITS("InterfaceEquations_InterpolationSetsNumberSet")
    RETURN
999 ERRORS("InterfaceEquations_InterpolationSetsNumberSet",err,error)
    EXITS("InterfaceEquations_InterpolationSetsNumberSet")
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_InterpolationSetsNumberSet

  !
  !================================================================================================================================
  !

  !>Finalise the interface equations and deallocate all memory.
  SUBROUTINE InterfaceEquations_Finalise(interfaceEquations,err,error,*)

    !Argument variables
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations !<A pointer to the interface equations to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("InterfaceEquations_Finalise",err,error,*999)

    IF(ASSOCIATED(interfaceEquations)) THEN
      CALL InterfaceEquations_InterpolationFinalise(interfaceEquations%interpolation,err,error,*999)
      IF(ASSOCIATED(interfaceEquations%interfaceMapping)) &
        & CALL InterfaceMapping_Destroy(interfaceEquations%interfaceMapping,err,error,*999)
      IF(ASSOCIATED(interfaceEquations%interfaceMatrices)) &
        & CALL InterfaceMatrices_Destroy(interfaceEquations%interfaceMatrices,err,error,*999)
      DEALLOCATE(interfaceEquations)
    ENDIF
       
    EXITS("InterfaceEquations_Finalise")
    RETURN
999 ERRORSEXITS("InterfaceEquations_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the interface equations for an interface condition
  SUBROUTINE InterfaceEquations_Initialise(interfaceCondition,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition to initialise the interface equations for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
 
    ENTERS("InterfaceEquations_Initialise",err,error,*998)

    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*998)
    IF(ASSOCIATED(interfaceCondition%interfaceEquations)) &
      & CALL FlagError("Interface equations is already associated for this interface condition.",err,error,*998)
    
    ALLOCATE(interfaceCondition%interfaceEquations,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate interface equations.",err,error,*999)
    interfaceCondition%interfaceEquations%interfaceCondition=>interfaceCondition
    interfaceCondition%interfaceEquations%linearity=INTERFACE_EQUATIONS_LINEAR
    interfaceCondition%interfaceEquations%timeDependence=INTERFACE_EQUATIONS_STATIC
    interfaceCondition%interfaceEquations%outputType=INTERFACE_EQUATIONS_NO_OUTPUT
    interfaceCondition%interfaceEquations%sparsityType=INTERFACE_EQUATIONS_SPARSE_MATRICES
    NULLIFY(interfaceCondition%interfaceEquations%interpolation)
    NULLIFY(interfaceCondition%interfaceEquations%interfaceMapping)
    NULLIFY(interfaceCondition%interfaceEquations%interfaceMatrices)
    interfaceCondition%interfaceEquations%interfaceEquationsFinished=.FALSE.
    CALL InterfaceEquations_InterpolationInitialise(interfaceCondition%interfaceEquations,err,error,*999)
        
    EXITS("InterfaceEquations_Initialise")
    RETURN
999 CALL InterfaceEquations_Finalise(interfaceCondition%interfaceEquations,dummyErr,dummyError,*998)
998 ERRORSEXITS("InterfaceEquations_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises the interface equations interpolation and deallocates all memory.
  SUBROUTINE InterfaceEquations_InterpolationFinalise(interfaceEquationsInterpolation,err,error,*)

    !Argument variables
    TYPE(InterfaceEquationsInterpolationType), POINTER :: interfaceEquationsInterpolation!<A pointer to the interface equations interpolation to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: variableIdx
 
    ENTERS("InterfaceEquations_InterpolationFinalise",err,error,*999)

    IF(ASSOCIATED(interfaceEquationsInterpolation)) THEN
      CALL InterfaceEquations_DomainInterpolationFinalise(interfaceEquationsInterpolation%interfaceInterpolation,err,error,*999)
      IF(ALLOCATED(interfaceEquationsInterpolation%variableInterpolation)) THEN
        DO variableIdx=1,SIZE(interfaceEquationsInterpolation%variableInterpolation,1)
          CALL InterfaceEquations_DomainInterpolationFinalise(interfaceEquationsInterpolation% &
            & variableInterpolation(variableIdx),err,error,*999)
        ENDDO !variableIdx
        DEALLOCATE(interfaceEquationsInterpolation%variableInterpolation)
      ENDIF
    ENDIF
       
    EXITS("InterfaceEquations_InterpolationFinalise")
    RETURN
999 ERRORSEXITS("InterfaceEquations_InterpolationFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_InterpolationFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the interface equations interpolation.
  SUBROUTINE InterfaceEquations_InterpolationInitialise(interfaceEquations,err,error,*)

    !Argument variables
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations !<A pointer to the interface equations to initialise the interpolation for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,variableIdx
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(InterfaceDependentType), POINTER :: interfaceDependent
    TYPE(VARYING_STRING) :: dummyError
 
    ENTERS("InterfaceEquations_InterpolationInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(interfaceEquations)) CALL FlagError("Interface equations is not associated.",err,error,*998)
    IF(ASSOCIATED(interfaceEquations%interpolation)) &
      & CALL FlagError("Interface equations interpolation is already associated.",err,error,*998)

    NULLIFY(interfaceCondition)
    CALL InterfaceEquations_InterfaceConditionGet(interfaceEquations,interfaceCondition,err,error,*998)
    NULLIFY(interfaceDependent)
    CALL InterfaceCondition_InterfaceDependentGet(interfaceCondition,interfaceDependent,err,error,*998)
    
    ALLOCATE(interfaceEquations%interpolation,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate interface equations interpolation.",err,error,*999)
    interfaceEquations%interpolation%interfaceEquations=>interfaceEquations
    CALL InterfaceEquations_DomainInterpolationInitialise(interfaceEquations%interpolation%interfaceInterpolation,err,error,*999)
    interfaceEquations%interpolation%interfaceInterpolation%interpolation=>interfaceEquations%interpolation
    ALLOCATE(interfaceEquations%interpolation%variableInterpolation(interfaceDependent%numberOfDependentVariables),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate interface equations interpolation variable interpolation.",err,error,*999)
    DO variableIdx=1,interfaceDependent%numberOfDependentVariables
      CALL InterfaceEquations_DomainInterpolationInitialise(interfaceEquations%interpolation% &
        & variableInterpolation(variableIdx),err,error,*999)
      interfaceEquations%interpolation%variableInterpolation(variableIdx)%interpolation=>interfaceEquations%interpolation
    ENDDO !variableIdx
       
    EXITS("InterfaceEquations_InterpolationInitialise")
    RETURN
999 CALL InterfaceEquations_InterpolationFinalise(interfaceEquations%interpolation,dummyErr,dummyError,*998)
998 ERRORSEXITS("InterfaceEquations_InterpolationInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_InterpolationInitialise

  !
  !================================================================================================================================
  !

  !>Finalises the interface equations interpolation set and deallocates all memory.
  SUBROUTINE InterfaceEquations_InterpolationSetFinalise(interpolationSet,err,error,*)

    !Argument variables
    TYPE(InterfaceEquationsInterpolationSetType) :: interpolationSet !<The interpolation set to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceEquations_InterpolationSetFinalise",err,error,*999)

    CALL Field_InterpolationParametersFinalise(interpolationSet%interpolationParameters,err,error,*999)
    CALL Field_InterpolatedPointsFinalise(interpolationSet%interpolatedPoint,err,error,*999)
    CALL Field_InterpolatedPointsMetricsFinalise(interpolationSet%interpolatedPointMetrics,err,error,*999)
       
    EXITS("InterfaceEquations_InterpolationSetFinalise")
    RETURN
999 ERRORSEXITS("InterfaceEquations_InterpolationSetFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_InterpolationSetFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the interface equations interpolation set.
  SUBROUTINE InterfaceEquations_InterpolationSetInitialise(interpolationSet,err,error,*)

    !Argument variables
    TYPE(InterfaceEquationsInterpolationSetType) :: interpolationSet !<The interpolation set to intialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceEquations_InterpolationSetInitialise",err,error,*999)

    NULLIFY(interpolationSet%interpolationParameters)
    NULLIFY(interpolationSet%interpolatedPoint)
    NULLIFY(interpolationSet%interpolatedPointMetrics)
       
    EXITS("InterfaceEquations_InterpolationSetInitialise")
    RETURN
999 ERRORS("InterfaceEquations_InterpolationSetInitialise",err,error)
    EXITS("InterfaceEquations_InterpolationSetInitialise")
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_InterpolationSetInitialise

  !
  !================================================================================================================================
  !

  !>Gets the output type for interface equations.
  SUBROUTINE InterfaceEquations_OutputTypeGet(interfaceEquations,outputType,err,error,*)

    !Argument variables
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations !<A pointer to the interface equations to get the output type for
    INTEGER(INTG), INTENT(OUT) :: outputType !<On exit, the output type of the interface equations \see InterfaceEquationsRoutines_OutputTypes,InterfaceEquationsRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceEquations_OutputTypeGet",err,error,*999)

    CALL InterfaceEquations_AssertIsFinished(interfaceEquations,err,error,*999)
    
    outputType=interfaceEquations%outputType
       
    EXITS("InterfaceEquations_OutputTypeGet")
    RETURN
999 ERRORSEXITS("InterfaceEquations_OutputTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_OutputTypeGet
  
  !
  !================================================================================================================================
  !

  !>Gets the time dependence type for an interface matrix
  SUBROUTINE InterfaceEquations_MatrixTimeDependenceTypeGet0(interfaceEquations,interfaceMatrixIdx,hasTranspose, &
    & timeDependenceType,err,error,*)
    
    !Argument variables
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations !<A pointer to the interface equations containing the interface matrix to get the time dependence type for.
    INTEGER(INTG), INTENT(IN) :: interfaceMatrixIdx !<The index of the interface matrix in the interface equations to get the time dependence for.
    LOGICAL, INTENT(IN) :: hasTranspose !<Is .TRUE. if the interface matrix has a transpose. .FALSE. if not. 
    INTEGER(INTG), INTENT(OUT) :: timeDependenceType !<The interface matrix time dependence type to get. If hasTranspose is .TRUE. then two timeDependenceTypes are required. The first one for the the interface matrix and the second one for the transposed matrix. \see InterfaceMatricesRoutines_InterfaceMatricesTimeDependenceTypes
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string    
    !Local variables
    INTEGER(INTG) :: timeDependenceTypes(1)
    
    ENTERS("InterfaceEquations_MatrixTimeDependenceTypeGet0",err,error,*999)

    CALL InterfaceEquations_MatrixTimeDependenceTypeGet1(interfaceEquations,interfaceMatrixIdx,hasTranspose,timeDependenceTypes, &
      & err,error,*999)
    timeDependenceType=timeDependenceTypes(1)
    
    EXITS("InterfaceEquations_MatrixTimeDependenceTypeGet0")
    RETURN
999 ERRORS("InterfaceEquations_MatrixTimeDependenceTypeGet0",Err,Error)
    EXITS("InterfaceEquations_MatrixTimeDependenceTypeGet0")
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_MatrixTimeDependenceTypeGet0

  !
  !================================================================================================================================
  !

  !>Gets the time dependence type for an interface matrix
  SUBROUTINE InterfaceEquations_MatrixTimeDependenceTypeGet1(interfaceEquations,interfaceMatrixIdx,hasTranspose, &
    & timeDependenceTypes,err,error,*)
    
    !Argument variables
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations !<A pointer to the interface equations containing the interface matrix to set the time dependence type for.
    INTEGER(INTG), INTENT(IN) :: interfaceMatrixIdx !<The index of the interface matrix in the interface equations to set the time dependence for.
    LOGICAL, INTENT(IN) :: hasTranspose !<Is .TRUE. if the interface matrix has a transpose. .FALSE. if not. 
    INTEGER(INTG), INTENT(OUT) :: timeDependenceTypes(:) !<timeDependenceTypes(transposeIdx). The interface matrix time dependence type to set. If hasTranspose is .TRUE. then two timeDependenceTypes are required. The first one for the the interface matrix and the second one for the transposed matrix. \see InterfaceMatricesRoutines_InterfaceMatricesTimeDependenceTypes
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string    
    !Local variables
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("InterfaceEquations_MatrixTimeDependenceTypeGet",err,error,*999)
    
    IF(.NOT.ASSOCIATED(interfaceEquations)) CALL FlagError("Interface equations not associated.",err,error,*999)
    NULLIFY(interfaceMatrices)
    CALL InterfaceEquations_InterfaceMatricesGet(interfaceEquations,interfaceMatrices,err,error,*999)
    NULLIFY(interfaceMatrix)
    CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,interfaceMatrixIdx,interfaceMatrix,err,error,*999)
    IF(hasTranspose) THEN
      IF(.NOT.SIZE(timeDependenceTypes,1)==2) THEN
        localError="The size of the time dependence types array of "// &
          & TRIM(NumberToVString(SIZE(timeDependenceTypes,1),"*",err,error))// &
          & " is invalid. The size needs to be 2 for an interface matrix that has a transpose."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      IF(.NOT.SIZE(timeDependenceTypes,1)==1) THEN
        localError="The size of the time dependence types array of "// &
          & TRIM(NumberToVString(SIZE(timeDependenceTypes,1),"*",err,error))// &
          & " is invalid. The size needs to be 1 for an interface matrix that has a transpose."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDIF

    timeDependenceTypes(1)=interfaceMatrix%interfaceMatrixTimeDependenceType
    IF(hasTranspose) THEN
      IF(interfaceMatrix%hasTranspose) THEN
        timeDependenceTypes(2)=interfaceMatrix%interfaceMatrixTransposeTimeDependenceType
      ELSE
        localError="Can not set the tranpose time dependence time for interface matrix index "// &
          & TRIM(NumberToVString(interfaceMatrixIdx,"*",err,error))//" as that matrix does not have a transpose."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDIF
    
    EXITS("InterfaceEquations_MatrixTimeDependenceTypeGet1")
    RETURN
999 ERRORS("InterfaceEquations_MatrixTimeDependenceTypeGet1",Err,Error)
    EXITS("InterfaceEquations_MatrixTimeDependenceTypeGet1")
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_MatrixTimeDependenceTypeGet1

  !
  !================================================================================================================================
  !

  !>Sets the time dependence type for an interface matrix
  SUBROUTINE InterfaceEquations_MatrixTimeDependenceTypeSet0(interfaceEquations,interfaceMatrixIdx,hasTranspose, &
    & timeDependenceType,err,error,*)
    
    !Argument variables
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations !<A pointer to the interface equations containing the interface matrix to set the time dependence type for.
    INTEGER(INTG), INTENT(IN) :: interfaceMatrixIdx !<The index of the interface matrix in the interface equations to set the time dependence for.
    LOGICAL, INTENT(IN) :: hasTranspose !<Is .TRUE. if the interface matrix has a transpose. .FALSE. if not. 
    INTEGER(INTG), INTENT(IN) :: timeDependenceType !<The interface matrix time dependence type to set. If hasTranspose is .TRUE. then two timeDependenceTypes are required. The first one for the the interface matrix and the second one for the transposed matrix. \see InterfaceMatricesRoutines_InterfaceMatricesTimeDependenceTypes
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string    
    !Local variables
    
    ENTERS("InterfaceEquations_MatrixTimeDependenceTypeSet0",err,error,*999)

    CALL InterfaceEquations_MatrixTimeDependenceTypeSet1(interfaceEquations,interfaceMatrixIdx,hasTranspose,[timeDependenceType], &
      & err,error,*999)
    
    EXITS("InterfaceEquations_MatrixTimeDependenceTypeSet0")
    RETURN
999 ERRORS("InterfaceEquations_MatrixTimeDependenceTypeSet0",Err,Error)
    EXITS("InterfaceEquations_MatrixTimeDependenceTypeSet0")
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_MatrixTimeDependenceTypeSet0

  !
  !================================================================================================================================
  !

  !>Sets the time dependence type for an interface matrix
  SUBROUTINE InterfaceEquations_MatrixTimeDependenceTypeSet1(interfaceEquations,interfaceMatrixIdx,hasTranspose, &
    & timeDependenceTypes,err,error,*)
    
    !Argument variables
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations !<A pointer to the interface equations containing the interface matrix to set the time dependence type for.
    INTEGER(INTG), INTENT(IN) :: interfaceMatrixIdx !<The index of the interface matrix in the interface equations to set the time dependence for.
    LOGICAL, INTENT(IN) :: hasTranspose !<Is .TRUE. if the interface matrix has a transpose. .FALSE. if not. 
    INTEGER(INTG), INTENT(IN) :: timeDependenceTypes(:) !<timeDependenceTypes(transposeIdx). The interface matrix time dependence type to set. If hasTranspose is .TRUE. then two timeDependenceTypes are required. The first one for the the interface matrix and the second one for the transposed matrix. \see InterfaceMatricesRoutines_InterfaceMatricesTimeDependenceTypes
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string    
    !Local variables
    INTEGER(INTG) :: transposeIdx
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("InterfaceEquations_MatrixTimeDependenceTypeSet",err,error,*999)
    
    IF(.NOT.ASSOCIATED(interfaceEquations)) CALL FlagError("Interface equations not associated.",err,error,*999)
    
    NULLIFY(interfaceMatrices)
    CALL InterfaceEquations_InterfaceMatricesGet(interfaceEquations,interfaceMatrices,err,error,*999)
    NULLIFY(interfaceMatrix)
    CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,interfaceMatrixIdx,interfaceMatrix,err,error,*999)
    IF(hasTranspose) THEN
      IF(.NOT.SIZE(timeDependenceTypes,1)==2) THEN
        localError="The size of the time dependence types array of "// &
          & TRIM(NumberToVString(SIZE(timeDependenceTypes,1),"*",err,error))// &
          & " is invalid. The size needs to be 2 for an interface matrix that has a transpose."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      IF(.NOT.SIZE(timeDependenceTypes,1)==1) THEN
        localError="The size of the time dependence types array of "// &
          & TRIM(NumberToVString(SIZE(timeDependenceTypes,1),"*",err,error))// &
          & " is invalid. The size needs to be 1 for an interface matrix that has a transpose."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDIF
    DO transposeIdx=1,SIZE(timeDependenceTypes,1)
      SELECT CASE(timeDependenceTypes(transposeIdx))
      CASE(INTERFACE_MATRIX_STATIC,INTERFACE_MATRIX_QUASI_STATIC,INTERFACE_MATRIX_FIRST_ORDER_DYNAMIC, &
        & INTERFACE_MATRIX_SECOND_ORDER_DYNAMIC)
        !Good input
      CASE DEFAULT
        localError="The specified interface matrix time dependence type of "// &
          & TRIM(NumberToVString(timeDependenceTypes(transposeIdx),"*",err,error))//" at position "// &
          & TRIM(NumberToVString(transposeIdx,"*",err,error))//" in the array is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDDO !transposeIdx

    interfaceMatrix%interfaceMatrixTimeDependenceType=timeDependenceTypes(1)
    IF(hasTranspose) THEN
      IF(interfaceMatrix%hasTranspose) THEN
        interfaceMatrix%interfaceMatrixTransposeTimeDependenceType=timeDependenceTypes(2)
      ELSE
        localError="Can not set the tranpose time dependence time for interface matrix index "// &
          & TRIM(NumberToVString(interfaceMatrixIdx,"*",err,error))//" as that matrix does not have a transpose."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDIF
    
    EXITS("InterfaceEquations_MatrixTimeDependenceTypeSet1")
    RETURN
999 ERRORS("InterfaceEquations_MatrixTimeDependenceTypeSet1",Err,Error)
    EXITS("InterfaceEquations_MatrixTimeDependenceTypeSet1")
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_MatrixTimeDependenceTypeSet1

  !
  !================================================================================================================================
  !

  !>Sets/changes the output type for the interface equations.
  SUBROUTINE InterfaceEquations_OutputTypeSet(interfaceEquations,outputType,err,error,*)

    !Argument variables
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations !<A pointer to the interface equations to set the output type for
    INTEGER(INTG), INTENT(IN) :: outputType !<The output type to set \see InterfaceEquationsRoutines_OutputTypes,InterfaceEquationsRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("InterfaceEquations_OutputTypeSet",err,error,*999)

    CALL InterfaceEquations_AssertNotFinished(interfaceEquations,err,error,*999)
    
    SELECT CASE(outputType)
    CASE(INTERFACE_EQUATIONS_NO_OUTPUT)
      interfaceEquations%outputType=INTERFACE_EQUATIONS_NO_OUTPUT
    CASE(INTERFACE_EQUATIONS_TIMING_OUTPUT)
      interfaceEquations%outputType=INTERFACE_EQUATIONS_TIMING_OUTPUT
    CASE(INTERFACE_EQUATIONS_MATRIX_OUTPUT)
      interfaceEquations%outputType=INTERFACE_EQUATIONS_MATRIX_OUTPUT
    CASE(INTERFACE_EQUATIONS_ELEMENT_MATRIX_OUTPUT)
      interfaceEquations%outputType=INTERFACE_EQUATIONS_ELEMENT_MATRIX_OUTPUT
    CASE DEFAULT
      localError="The specified output type of "//TRIM(NumberToVString(outputType,"*",err,error))//" is invalid"
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("InterfaceEquations_OutputTypeSet")
    RETURN
999 ERRORSEXITS("InterfaceEquations_OutputTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_OutputTypeSet

  !
  !================================================================================================================================
  !

  !>Gets the sparsity type for interface equations.
  SUBROUTINE InterfaceEquations_SparsityTypeGet(interfaceEquations,sparsityType,err,error,*)

    !Argument variables
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations !<A pointer to the interface equations to get the sparsity type for
    INTEGER(INTG), INTENT(OUT) :: sparsityType !<On exit, the sparsity type of the interface equations. \see InterfaceEquationsRoutines_SparsityTypes,InterfaceEquationsRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceEquations_SparsityTypeGet",err,error,*999)

    CALL InterfaceEquations_AssertIsFinished(interfaceEquations,err,error,*999)
    
    sparsityType=interfaceEquations%sparsityType
       
    EXITS("InterfaceEquations_SparsityTypeGet")
    RETURN
999 ERRORSEXITS("InterfaceEquations_SparsityTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_SparsityTypeGet
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the sparsity type for the interface equations.
  SUBROUTINE InterfaceEquations_SparsityTypeSet(interfaceEquations,sparsityType,err,error,*)

    !Argument variables
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations !<A pointer to the interface equations to set the sparsity type for
    INTEGER(INTG), INTENT(IN) :: sparsityType !<The sparsity type to set \see InterfaceEquationsRoutines_SparsityTypes,InterfaceEquationsRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("InterfaceEquations_SparsityTypeSet",err,error,*999)

    CALL InterfaceEquations_AssertNotFinished(interfaceEquations,err,error,*999)

    SELECT CASE(sparsityType)
    CASE(INTERFACE_EQUATIONS_SPARSE_MATRICES)
      interfaceEquations%sparsityType=INTERFACE_EQUATIONS_SPARSE_MATRICES
    CASE(INTERFACE_EQUATIONS_FULL_MATRICES)
      interfaceEquations%sparsityType=INTERFACE_EQUATIONS_FULL_MATRICES
    CASE DEFAULT
      localError="The specified sparsity type of "//TRIM(NumberToVString(sparsityType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
        
    EXITS("InterfaceEquations_SparsityTypeSet")
    RETURN
999 ERRORSEXITS("InterfaceEquations_SparsityTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_SparsityTypeSet
  
  !
  !================================================================================================================================
  !

  !>Gets the linearity type for interface equations.
  SUBROUTINE InterfaceEquations_LinearityTypeGet(interfaceEquations,linearityType,err,error,*)

    !Argument variables
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations !<A pointer to the interface equations to get the linearity for
    INTEGER(INTG), INTENT(OUT) :: linearityType !<On exit, the linearity type of the interface equations. \see InterfaceEquations_LinearityTypes
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceEquations_LinearityTypeGet",err,error,*999)

    CALL InterfaceEquations_AssertIsFinished(interfaceEquations,err,error,*999)
    
    linearityType=interfaceEquations%linearity
       
    EXITS("InterfaceEquations_LinearityTypeGet")
    RETURN
999 ERRORSEXITS("InterfaceEquations_LinearityTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_LinearityTypeGet
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the linearity type for interface equations.
  SUBROUTINE InterfaceEquations_LinearityTypeSet(interfaceEquations,linearityType,err,error,*)

    !Argument variables
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations !<A pointer to the interface equations to set the linearity for
    INTEGER(INTG), INTENT(IN) :: linearityType !<The linearity type to set \see InterfaceEquations_LinearityTypes
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("InterfaceEquations_LinearityTypeSet",err,error,*999)

    CALL InterfaceEquations_AssertNotFinished(interfaceEquations,err,error,*999)
    
    SELECT CASE(linearityType)
    CASE(INTERFACE_EQUATIONS_LINEAR)
      interfaceEquations%linearity=INTERFACE_EQUATIONS_LINEAR
    CASE(INTERFACE_EQUATIONS_NONLINEAR)
      interfaceEquations%linearity=INTERFACE_EQUATIONS_NONLINEAR
    CASE(INTERFACE_EQUATIONS_NONLINEAR_BCS)
      interfaceEquations%linearity=INTERFACE_EQUATIONS_NONLINEAR_BCS
    CASE DEFAULT
      localError="The specified linearity type of "//TRIM(NumberToVString(linearityType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("InterfaceEquations_LinearityTypeSet")
    RETURN
999 ERRORSEXITS("InterfaceEquations_LinearityTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_LinearityTypeSet
  
  !
  !================================================================================================================================
  !

  !>Gets the time dependence type for interface equations.
  SUBROUTINE InterfaceEquations_TimeDependenceTypeGet(interfaceEquations,timeDependenceType,err,error,*)

    !Argument variables
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations !<A pointer to the interface equations to get the output type for
    INTEGER(INTG), INTENT(OUT) :: timeDependenceType !<On exit, the time dependence type of the interface equations \see InterfaceEquations_TimeDependenceTypes
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceEquations_TimeDependenceTypeGet",err,error,*999)

    CALL InterfaceEquations_AssertIsFinished(interfaceEquations,err,error,*999)

    timeDependenceType=interfaceEquations%timeDependence
       
    EXITS("InterfaceEquations_TimeDependenceTypeGet")
    RETURN
999 ERRORSEXITS("InterfaceEquations_TimeDependenceTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_TimeDependenceTypeGet
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the time dependence type for interface equations.
  SUBROUTINE InterfaceEquations_TimeDependenceTypeSet(interfaceEquations,timeDependenceType,err,error,*)

    !Argument variables
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations !<A pointer to the interface equations to set the linearity for
    INTEGER(INTG), INTENT(IN) :: timeDependenceType !<The time dependence type to set \see InterfaceEquations_TimeDependenceTypes
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("InterfaceEquations_TimeDependenceTypeSet",err,error,*999)

    CALL InterfaceEquations_AssertNotFinished(interfaceEquations,err,error,*999)
    
    SELECT CASE(timeDependenceType)
    CASE(INTERFACE_EQUATIONS_STATIC)
      interfaceEquations%timeDependence=INTERFACE_EQUATIONS_STATIC
    CASE(INTERFACE_EQUATIONS_QUASISTATIC)
      interfaceEquations%timeDependence=INTERFACE_EQUATIONS_QUASISTATIC
    CASE(INTERFACE_EQUATIONS_FIRST_ORDER_DYNAMIC)
      interfaceEquations%timeDependence=INTERFACE_EQUATIONS_FIRST_ORDER_DYNAMIC
    CASE(INTERFACE_EQUATIONS_SECOND_ORDER_DYNAMIC)
      interfaceEquations%timeDependence=INTERFACE_EQUATIONS_SECOND_ORDER_DYNAMIC
    CASE DEFAULT
      localError="The specified time dependence type of "//TRIM(NumberToVString(timeDependenceType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("InterfaceEquations_TimeDependenceTypeSet")
    RETURN
999 ERRORSEXITS("InterfaceEquations_TimeDependenceTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_TimeDependenceTypeSet

  !
  !================================================================================================================================
  !
  
END MODULE InterfaceEquationsRoutines
