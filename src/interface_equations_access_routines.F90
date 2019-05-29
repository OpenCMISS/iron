!> \file
!> \author Chris Bradley
!> \brief This module contains all interface equations access method routines.
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

!> This module contains all interface equations access method routines.
MODULE InterfaceEquationsAccessRoutines
  
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

  !> \addtogroup OpenCMISS_InterfaceEquationsConstants OpenCMISS::Iron::InterfaceEquations::Constants
  !> \brief Matrix vector constants.
  !>@{
  !> \addtogroup InterfaceEquations_LinearityTypes InterfaceEquations::Constants::LinearityTypes
  !> \brief The interface equations linearity types
  !> \see InterfaceEquations
  !>@{
  INTEGER(INTG), PARAMETER :: NUMBER_OF_INTERFACE_EQUATIONS_LINEARITY_TYPES=3 !<The number of interface equations linearity types defined. \see InterfaceEquations_LinearityTypes,OpenCMISS_InterfaceEquationsConstants
  INTEGER(INTG), PARAMETER :: INTERFACE_EQUATIONS_LINEAR=1 !<The interface equations are linear. \see InterfaceEquations_LinearityTypes,OpenCMISS_InterfaceEquationsConstants
  INTEGER(INTG), PARAMETER :: INTERFACE_EQUATIONS_NONLINEAR=2 !<The interface equations are non-linear. \see InterfaceEquations_LinearityTypes,OpenCMISS_InterfaceEquationsConstants
  INTEGER(INTG), PARAMETER :: INTERFACE_EQUATIONS_NONLINEAR_BCS=3 !<The interface equations have non-linear boundary conditions. \see InterfaceEquations_LinearityTypes,OpenCMISS_InterfaceEquationsConstants
  !>@}
  !> \addtogroup InterfaceEquations_OutputTypes InterfaceEquations::Constants::OutputTypes
  !> \brief The interface equations output types
  !> \see InterfaceEquations
  !>@{
  INTEGER(INTG), PARAMETER :: INTERFACE_EQUATIONS_NO_OUTPUT=0 !<No output. \see InterfaceEquations_OutputTypes,OpenCMISS_InterfaceEquationsConstants
  INTEGER(INTG), PARAMETER :: INTERFACE_EQUATIONS_TIMING_OUTPUT=1 !<Timing information output. \see InterfaceEquations_OutputTypes,OpenCMISS_InterfaceEquationsConstants
  INTEGER(INTG), PARAMETER :: INTERFACE_EQUATIONS_MATRIX_OUTPUT=2 !<All below and equation matrices output. \see InterfaceEquations_OutputTypes,OpenCMISS_InterfaceEquationsConstants
  INTEGER(INTG), PARAMETER :: INTERFACE_EQUATIONS_ELEMENT_MATRIX_OUTPUT=3 !<All below and element matrices output. \see InterfaceEquations_OutputTypes,OpenCMISS_InterfaceEquationsConstants
  !>@}
  !> \addtogroup InterfaceEquations_SparsityTypes InterfaceEquations::Constants::SparsityTypes
  !> \brief Interface equations matrices sparsity types
  !> \see InterfaceEquations
  !>@{
  INTEGER(INTG), PARAMETER :: INTERFACE_EQUATIONS_SPARSE_MATRICES=1 !<Use sparse matrices for the interface equations. \see InterfaceEquations_SparsityTypes,OpenCMISS_InterfaceEquationsConstants
  INTEGER(INTG), PARAMETER :: INTERFACE_EQUATIONS_FULL_MATRICES=2 !<Use fully populated matrices for the interface equations. \see InterfaceEquations_SparsityTypes,OpenCMISS_InterfaceEquationsConstants
  !>@}
  !> \addtogroup InterfaceEquations_TimeDependenceTypes InterfaceEquations::Constants::TimeDependenceTypes
  !> \brief The interface equations time dependence type parameters
  !> \see InterfaceEquations
  !>@{
  INTEGER(INTG), PARAMETER :: NUMBER_OF_INTERFACE_EQUATIONS_TIME_TYPES=5 !<The number of interface equations time dependence types defined. \see InterfaceEquations_TimeDependenceTypes,OpenCMISS_InterfaceEquationsConstants
  INTEGER(INTG), PARAMETER :: INTERFACE_EQUATIONS_STATIC=1 !<The interface conditions are static and have no time dependence. \see InterfaceEquations_TimeDependenceTypes,OpenCMISS_InterfaceEquationsConstants 
  INTEGER(INTG), PARAMETER :: INTERFACE_EQUATIONS_QUASISTATIC=2 !<The interface conditions are quasi-static. \see InterfaceEquations_TimeDependenceTypes,OpenCMISS_InterfaceEquationsConstants
  INTEGER(INTG), PARAMETER :: INTERFACE_EQUATIONS_FIRST_ORDER_DYNAMIC=3 !<The interface conditions are first order dynamic. \see InterfaceEquations_TimeDependenceTypes,OpenCMISS_InterfaceEquationsConstants
  INTEGER(INTG), PARAMETER :: INTERFACE_EQUATIONS_SECOND_ORDER_DYNAMIC=4 !<The interface conditions are a second order dynamic. \see InterfaceEquations_TimeDependenceTypes,OpenCMISS_InterfaceEquationsConstants
  INTEGER(INTG), PARAMETER :: INTERFACE_EQUATIONS_TIME_STEPPING=5 !<The interface conditions are for time stepping. \see InterfaceEquations_TimeDependenceTypes,OpenCMISS_InterfaceEquationsConstants
  !>@}
  !>@}  
  
  !Module types

  !Module variables

  !Interfaces

  PUBLIC NUMBER_OF_INTERFACE_EQUATIONS_LINEARITY_TYPES,INTERFACE_EQUATIONS_LINEAR,INTERFACE_EQUATIONS_NONLINEAR, &
    & INTERFACE_EQUATIONS_NONLINEAR_BCS
  
  PUBLIC INTERFACE_EQUATIONS_NO_OUTPUT,INTERFACE_EQUATIONS_TIMING_OUTPUT,INTERFACE_EQUATIONS_MATRIX_OUTPUT, &
    & INTERFACE_EQUATIONS_ELEMENT_MATRIX_OUTPUT

  PUBLIC INTERFACE_EQUATIONS_SPARSE_MATRICES,INTERFACE_EQUATIONS_FULL_MATRICES

  PUBLIC NUMBER_OF_INTERFACE_EQUATIONS_TIME_TYPES,INTERFACE_EQUATIONS_STATIC,INTERFACE_EQUATIONS_QUASISTATIC, &
    & INTERFACE_EQUATIONS_FIRST_ORDER_DYNAMIC,INTERFACE_EQUATIONS_SECOND_ORDER_DYNAMIC,INTERFACE_EQUATIONS_TIME_STEPPING

  PUBLIC InterfaceDomainInterpolation_DependentFieldGet
  
  PUBLIC InterfaceDomainInterpolation_DependentInterpSetGet
  
  PUBLIC InterfaceDomainInterpolation_GeometricFieldGet
  
  PUBLIC InterfaceDomainInterpolation_GeometricInterpSetGet
  
  PUBLIC InterfaceDomainInterpolation_PenaltyFieldGet
  
  PUBLIC InterfaceDomainInterpolation_PenaltyInterpSetGet

  PUBLIC InterfaceEquations_AssertIsFinished,InterfaceEquations_AssertNotFinished

  PUBLIC InterfaceEquations_EquationsInterpolationGet

  PUBLIC InterfaceEquations_InterfaceConditionGet

  PUBLIC InterfaceEquations_InterfaceMappingGet

  PUBLIC InterfaceEquations_InterfaceMatricesGet

  PUBLIC InterfaceEquationsInterpolation_InterfaceInterpGet

  PUBLIC InterfaceEquationsInterpolation_VariableInterpGet

  PUBLIC InterfaceInterpolationSet_InterpolationParametersGet

  PUBLIC InterfaceInterpolationSet_InterpolatedPointGet

  PUBLIC InterfaceInterpolationSet_InterpolatedPointMetricsGet


CONTAINS

  !
  !================================================================================================================================
  !

  !>Gets the dependent field for an interface domain interpolation.
  SUBROUTINE InterfaceDomainInterpolation_DependentFieldGet(domainInterpolation,dependentField,err,error,*)

    !Argument variables
    TYPE(InterfaceEquationsDomainInterpolationType), POINTER :: domainInterpolation !<A pointer to the interface domain interpolation to get the dependent field for
    TYPE(FieldType), POINTER :: dependentField !<On exit, a pointer to the dependent field for the specified interface domain interpolation. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceDomainInterpolation_DependentFieldGet",err,error,*998)

    IF(ASSOCIATED(dependentField)) CALL FlagError("Dependent field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(domainInterpolation)) CALL FlagError("Interface domain interpolation is not associated.",err,error,*999)

    dependentField=>domainInterpolation%dependentField
    IF(.NOT.ASSOCIATED(dependentField)) &
      & CALL FlagError("Interface domain interpolation dependent field is not associated.",err,error,*999)
       
    EXITS("InterfaceDomainInterface_DependentFieldGet")
    RETURN
999 NULLIFY(dependentField)
998 ERRORS("InterfaceDomainInterpolation_DependentFieldGet",err,error)
    EXITS("InterfaceDomainInterpolation_DependentFieldGet")
    RETURN 1
    
  END SUBROUTINE InterfaceDomainInterpolation_DependentFieldGet

  !
  !================================================================================================================================
  !

  !>Gets the dependent interpolation set for an interface domain interpolation.
  SUBROUTINE InterfaceDomainInterpolation_DependentInterpSetGet(domainInterpolation,interpolationSetIdx, &
    & dependentInterpolationSet,err,error,*)

    !Argument variables
    TYPE(InterfaceEquationsDomainInterpolationType), POINTER :: domainInterpolation !<A pointer to the interface domain interpolation to get the dependent field for
    INTEGER(INTG), INTENT(IN) :: interpolationSetIdx !<The interpoltion set index to get.
    TYPE(InterfaceEquationsInterpolationSetType), POINTER :: dependentInterpolationSet !<On exit, a pointer to the dependent interpolation set for the specified interface domain interpolation. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("InterfaceDomainInterpolation_DependentInterpSetGet",err,error,*998)

    IF(ASSOCIATED(dependentInterpolationSet)) CALL FlagError("Dependent interpolation set is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(domainInterpolation)) CALL FlagError("Interface domain interpolation is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(domainInterpolation%dependentInterpolation)) &
      & CALL FlagError("Interface domain interpolation dependent interpolation is not allocated.",err,error,*999)
    IF(interpolationSetIdx<=0.OR.interpolationSetIdx>domainInterpolation%numberOfDependentInterpolationSets) THEN
      localError="The specified interpolation set index of "//TRIM(NumberToVString(interpolationSetIdx,"*",err,error))// &
        & " is invalid. The interpolation set index should be >= 1 and <= "// & 
        & TRIM(NumberToVString(domainInterpolation%numberOfDependentInterpolationSets,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    dependentInterpolationSet=>domainInterpolation%dependentInterpolation(interpolationSetIdx)
    IF(.NOT.ASSOCIATED(dependentInterpolationSet)) THEN
      localError="The dependent interpolation set is not associated for interpolation set index "// &
        & TRIM(NumberToVString(interpolationSetIdx,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("InterfaceDomainInterface_DependentInterpSetGet")
    RETURN
999 NULLIFY(dependentInterpolationSet)
998 ERRORS("InterfaceDomainInterpolation_DependentInterpSetGet",err,error)
    EXITS("InterfaceDomainInterpolation_DependentInterpSetGet")
    RETURN 1
    
  END SUBROUTINE InterfaceDomainInterpolation_DependentInterpSetGet

  !
  !================================================================================================================================
  !

  !>Gets the geometric field for an interface domain interpolation.
  SUBROUTINE InterfaceDomainInterpolation_GeometricFieldGet(domainInterpolation,geometricField,err,error,*)

    !Argument variables
    TYPE(InterfaceEquationsDomainInterpolationType), POINTER :: domainInterpolation !<A pointer to the interface domain interpolation to get the geometric field for
    TYPE(FieldType), POINTER :: geometricField !<On exit, a pointer to the geometric field for the specified interface domain interpolation. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceDomainInterpolation_GeometricFieldGet",err,error,*998)

    IF(ASSOCIATED(geometricField)) CALL FlagError("Geometric field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(domainInterpolation)) CALL FlagError("Interface domain interpolation is not associated.",err,error,*999)

    geometricField=>domainInterpolation%geometricField
    IF(.NOT.ASSOCIATED(geometricField)) &
      & CALL FlagError("Interface domain interpolation geometric field is not associated.",err,error,*999)
       
    EXITS("InterfaceDomainInterface_GeometricFieldGet")
    RETURN
999 NULLIFY(geometricField)
998 ERRORS("InterfaceDomainInterpolation_GeometricFieldGet",err,error)
    EXITS("InterfaceDomainInterpolation_GeometricFieldGet")
    RETURN 1
    
  END SUBROUTINE InterfaceDomainInterpolation_GeometricFieldGet

  !
  !================================================================================================================================
  !

  !>Gets the geometric interpolation set for an interface domain interpolation.
  SUBROUTINE InterfaceDomainInterpolation_GeometricInterpSetGet(domainInterpolation,interpolationSetIdx, &
    & geometricInterpolationSet,err,error,*)

    !Argument variables
    TYPE(InterfaceEquationsDomainInterpolationType), POINTER :: domainInterpolation !<A pointer to the interface domain interpolation to get the geometric field for
    INTEGER(INTG), INTENT(IN) :: interpolationSetIdx !<The interpoltion set index to get.
    TYPE(InterfaceEquationsInterpolationSetType), POINTER :: geometricInterpolationSet !<On exit, a pointer to the geometric interpolation set for the specified interface domain interpolation. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("InterfaceDomainInterpolation_GeometricInterpSetGet",err,error,*998)

    IF(ASSOCIATED(geometricInterpolationSet)) CALL FlagError("Geometric interpolation set is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(domainInterpolation)) CALL FlagError("Interface domain interpolation is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(domainInterpolation%geometricInterpolation)) &
      & CALL FlagError("Interface domain interpolation geometric interpolation is not allocated.",err,error,*999)
    IF(interpolationSetIdx<=0.OR.interpolationSetIdx>domainInterpolation%numberOfGeometricInterpolationSets) THEN
      localError="The specified interpolation set index of "//TRIM(NumberToVString(interpolationSetIdx,"*",err,error))// &
        & " is invalid. The interpolation set index should be >= 1 and <= "// & 
        & TRIM(NumberToVString(domainInterpolation%numberOfGeometricInterpolationSets,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    geometricInterpolationSet=>domainInterpolation%geometricInterpolation(interpolationSetIdx)
    IF(.NOT.ASSOCIATED(geometricInterpolationSet)) THEN
      localError="The geometric interpolation set is not associated for interpolation set index "// &
        & TRIM(NumberToVString(interpolationSetIdx,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("InterfaceDomainInterface_GeometricInterpSetGet")
    RETURN
999 NULLIFY(geometricInterpolationSet)
998 ERRORS("InterfaceDomainInterpolation_GeometricInterpSetGet",err,error)
    EXITS("InterfaceDomainInterpolation_GeometricInterpSetGet")
    RETURN 1
    
  END SUBROUTINE InterfaceDomainInterpolation_GeometricInterpSetGet

  !
  !================================================================================================================================
  !

  !>Gets the penalty field for an interface domain interpolation.
  SUBROUTINE InterfaceDomainInterpolation_PenaltyFieldGet(domainInterpolation,penaltyField,err,error,*)

    !Argument variables
    TYPE(InterfaceEquationsDomainInterpolationType), POINTER :: domainInterpolation !<A pointer to the interface domain interpolation to get the penalty field for
    TYPE(FieldType), POINTER :: penaltyField !<On exit, a pointer to the penalty field for the specified interface domain interpolation. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceDomainInterpolation_PenaltyFieldGet",err,error,*998)

    IF(ASSOCIATED(penaltyField)) CALL FlagError("Penalty field is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(domainInterpolation)) CALL FlagError("Domain interpolation is not associated.",err,error,*999)

    penaltyField=>domainInterpolation%penaltyField
    IF(.NOT.ASSOCIATED(penaltyField)) &
      & CALL FlagError("Interface domain interpolation penalty field is not associated.",err,error,*999)
       
    EXITS("InterfaceDomainInterface_PenaltyFieldGet")
    RETURN
999 NULLIFY(penaltyField)
998 ERRORS("InterfaceDomainInterpolation_PenaltyFieldGet",err,error)
    EXITS("InterfaceDomainInterpolation_PenaltyFieldGet")
    RETURN 1
    
  END SUBROUTINE InterfaceDomainInterpolation_PenaltyFieldGet

  !
  !================================================================================================================================
  !

  !>Gets the penalty interpolation set for an interface domain interpolation.
  SUBROUTINE InterfaceDomainInterpolation_PenaltyInterpSetGet(domainInterpolation,interpolationSetIdx, &
    & penaltyInterpolationSet,err,error,*)

    !Argument variables
    TYPE(InterfaceEquationsDomainInterpolationType), POINTER :: domainInterpolation !<A pointer to the interface domain interpolation to get the penalty field for
    INTEGER(INTG), INTENT(IN) :: interpolationSetIdx !<The interpoltion set index to get.
    TYPE(InterfaceEquationsInterpolationSetType), POINTER :: penaltyInterpolationSet !<On exit, a pointer to the penalty interpolation set for the specified interface domain interpolation. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("InterfaceDomainInterpolation_PenaltyInterpSetGet",err,error,*998)

    IF(ASSOCIATED(penaltyInterpolationSet)) CALL FlagError("Penalty interpolation set is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(domainInterpolation)) CALL FlagError("Interface domain interpolation is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(domainInterpolation%penaltyInterpolation)) &
      & CALL FlagError("Interface domain interpolation penalty interpolation is not allocated.",err,error,*999)
    IF(interpolationSetIdx<=0.OR.interpolationSetIdx>domainInterpolation%numberOfPenaltyInterpolationSets) THEN
      localError="The specified interpolation set index of "//TRIM(NumberToVString(interpolationSetIdx,"*",err,error))// &
        & " is invalid. The interpolation set index should be >= 1 and <= "// & 
        & TRIM(NumberToVString(domainInterpolation%numberOfPenaltyInterpolationSets,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    penaltyInterpolationSet=>domainInterpolation%penaltyInterpolation(interpolationSetIdx)
    IF(.NOT.ASSOCIATED(penaltyInterpolationSet)) THEN
      localError="The penalty interpolation set is not associated for interpolation set index "// &
        & TRIM(NumberToVString(interpolationSetIdx,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("InterfaceDomainInterface_PenaltyInterpSetGet")
    RETURN
999 NULLIFY(penaltyInterpolationSet)
998 ERRORS("InterfaceDomainInterpolation_PenaltyInterpSetGet",err,error)
    EXITS("InterfaceDomainInterpolation_PenaltyInterpSetGet")
    RETURN 1
    
  END SUBROUTINE InterfaceDomainInterpolation_PenaltyInterpSetGet

  !
  !================================================================================================================================
  !

  !>Assert that an interface equations has been finished
  SUBROUTINE InterfaceEquations_AssertIsFinished(interfaceEquations,err,error,*)

    !Argument Variables
    TYPE(InterfaceEquationsType), POINTER, INTENT(IN) :: interfaceEquations !<The interface equations to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceEquations_AssertIsFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceEquations)) CALL FlagError("Interface equations is not associated.",err,error,*999)

    IF(.NOT.interfaceEquations%interfaceEquationsFinished) &
      & CALL FlagError("Interface equations has not been finished.",err,error,*999)
    
    EXITS("InterfaceEquations_AssertIsFinished")
    RETURN
999 ERRORSEXITS("InterfaceEquations_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_AssertIsFinished

  !
  !================================================================================================================================
  !

  !>Assert that an interface equations has not been finished
  SUBROUTINE InterfaceEquations_AssertNotFinished(interfaceEquations,err,error,*)

    !Argument Variables
    TYPE(InterfaceEquationsType), POINTER, INTENT(IN) :: interfaceEquations !<The interface equations to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceEquations_AssertNotFinished",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceEquations)) CALL FlagError("Interface equations is not associated.",err,error,*999)

    IF(interfaceEquations%interfaceEquationsFinished) &
      & CALL FlagError("Interface equations has already been finished.",err,error,*999)
    
    EXITS("InterfaceEquations_AssertNotFinished")
    RETURN
999 ERRORSEXITS("InterfaceEquations_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_AssertNotFinished

  !
  !================================================================================================================================
  !

  !>Gets the equations interpolation for an interface equations.
  SUBROUTINE InterfaceEquations_EquationsInterpolationGet(interfaceEquations,equationsInterpolation,err,error,*)

    !Argument variables
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations !<A pointer to the interface equations to get the equations interpolation for
    TYPE(InterfaceEquationsInterpolationType), POINTER :: equationsInterpolation !<On exit, a pointer to the equations interpolation for the specified interface equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceEquations_EquationsInterpolationGet",err,error,*998)

    IF(ASSOCIATED(equationsInterpolation)) CALL FlagError("Equations interpolation is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceEquations)) CALL FlagError("Interface equations is not associated.",err,error,*999)

    equationsInterpolation=>interfaceEquations%interpolation
    IF(.NOT.ASSOCIATED(equationsInterpolation)) &
      & CALL FlagError("Interface equations interpolation is not associated.",err,error,*999)
       
    EXITS("InterfaceEquations_EquationsInterpolationGet")
    RETURN
999 NULLIFY(equationsInterpolation)
998 ERRORS("InterfaceEquations_EquationsInterpolationGet",err,error)
    EXITS("InterfaceEquations_EquationsInterpolationGet")
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_EquationsInterpolationGet
  
  !
  !================================================================================================================================
  !

  !>Gets the interface condition for an interface equations.
  SUBROUTINE InterfaceEquations_InterfaceConditionGet(interfaceEquations,interfaceCondition,err,error,*)

    !Argument variables
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations !<A pointer to the interface equations to get the interface matrices for
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<On exit, a pointer to the interface condition for the specified interface equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceEquations_InterfaceConditionGet",err,error,*998)

    IF(ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceEquations)) CALL FlagError("Interface equations is not associated.",err,error,*999)

    interfaceCondition=>interfaceEquations%interfaceCondition
    IF(.NOT.ASSOCIATED(interfaceCondition)) &
      & CALL FlagError("Interface equations interface condition is not associated.",err,error,*999)
       
    EXITS("InterfaceEquations_InterfaceConditionGet")
    RETURN
999 NULLIFY(interfaceCondition)
998 ERRORSEXITS("InterfaceEquations_InterfaceConditionGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_InterfaceConditionGet
  
  !
  !================================================================================================================================
  !

  !>Gets the interface mapping for an interface equations.
  SUBROUTINE InterfaceEquations_InterfaceMappingGet(interfaceEquations,interfaceMapping,err,error,*)

    !Argument variables
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations !<A pointer to the interface equations to get the interface matrices for
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping !<On exit, a pointer to the interface mapping for the specified interface equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceEquations_InterfaceMappingGet",err,error,*998)

    IF(ASSOCIATED(interfaceMapping)) CALL FlagError("Interface mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceEquations)) CALL FlagError("Interface equations is not associated.",err,error,*999)

    interfaceMapping=>interfaceEquations%interfaceMapping
    IF(.NOT.ASSOCIATED(interfaceMapping)) &
      & CALL FlagError("Interface equations interface mapping is not associated.",err,error,*999)
       
    EXITS("InterfaceEquations_InterfaceMappingGet")
    RETURN
999 NULLIFY(interfaceMapping)
998 ERRORSEXITS("InterfaceEquations_InterfaceMappingGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_InterfaceMappingGet
  
  !
  !================================================================================================================================
  !

  !>Gets the interface matrices for an interface equations.
  SUBROUTINE InterfaceEquations_InterfaceMatricesGet(interfaceEquations,interfaceMatrices,err,error,*)

    !Argument variables
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations !<A pointer to the interface equations to get the interface matrices for
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices !<On exit, a pointer to the interface matrices for the specified interface equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceEquations_InterfaceMatricesGet",err,error,*998)

    IF(ASSOCIATED(interfaceMatrices)) CALL FlagError("Interface matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceEquations)) CALL FlagError("Interface equations is not associated.",err,error,*999)

    interfaceMatrices=>interfaceEquations%interfaceMatrices
    IF(.NOT.ASSOCIATED(interfaceMatrices)) &
      & CALL FlagError("Interface equations interface matrices is not associated.",err,error,*999)
       
    EXITS("InterfaceEquations_InterfaceMatricesGet")
    RETURN
999 NULLIFY(interfaceMatrices)
998 ERRORSEXITS("InterfaceEquations_InterfaceMatricesGet",err,error)
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_InterfaceMatricesGet

  !
  !================================================================================================================================
  !

  !>Gets the interface interpolation for an interface equations interpolation.
  SUBROUTINE InterfaceEquationsInterpolation_InterfaceInterpGet(interfaceEquationsInterpolation, &
    & interfaceInterpolation,err,error,*)

    !Argument variables
    TYPE(InterfaceEquationsInterpolationType), POINTER :: interfaceEquationsInterpolation !<A pointer to the interface equations interpolation to get the interface interpolation for
    TYPE(InterfaceEquationsDomainInterpolationType), POINTER :: interfaceInterpolation !<On exit, a pointer to the interface interpolation for the specified interface equations interpolation. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("InterfaceEquationsInterpolation_InterfaceInterpGet",err,error,*998)

    IF(ASSOCIATED(interfaceInterpolation)) CALL FlagError("Interface interpolation is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceEquationsInterpolation)) &
      & CALL FlagError("Interface equations interpolation is not associated.",err,error,*999)

    interfaceInterpolation=>interfaceEquationsInterpolation%interfaceInterpolation
    IF(.NOT.ASSOCIATED(interfaceInterpolation)) &
      & CALL FlagError("Interface equations interpolation interface interpolation is not associated.",err,error,*999)
       
    EXITS("InterfaceEquationsInterface_InterfaceInterpGet")
    RETURN
999 NULLIFY(interfaceInterpolation)
998 ERRORS("InterfaceEquationsInterpolation_InterfaceInterpGet",err,error)
    EXITS("InterfaceEquationsInterpolation_InterfaceInterpGet")
    RETURN 1
    
  END SUBROUTINE InterfaceEquationsInterpolation_InterfaceInterpGet

  !
  !================================================================================================================================
  !

  !>Gets the variable interpolation for an interface equations interpolation.
  SUBROUTINE InterfaceEquationsInterpolation_VariableInterpGet(interfaceEquationsInterpolation,variableIdx, &
    & variableInterpolation,err,error,*)

    !Argument variables
    TYPE(InterfaceEquationsInterpolationType), POINTER :: interfaceEquationsInterpolation !<A pointer to the interface equations interpolation to get the variable interpolation for
    INTEGER(INTG), INTENT(IN) :: variableIdx !<The index of the variable to get the interpolation for. 
    TYPE(InterfaceEquationsDomainInterpolationType), POINTER :: variableInterpolation !<On exit, a pointer to the variable interpolation for the specified interface equations interpolation. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("InterfaceEquationsInterpolation_VariableInterpGet",err,error,*998)

    IF(ASSOCIATED(variableInterpolation)) CALL FlagError("Variable interpolation is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceEquationsInterpolation)) &
      & CALL FlagError("Interface equations interpolation is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(interfaceEquationsInterpolation%variableInterpolation)) &
      & CALL FlagError("Interface equations variable interpolation is not allocated.",err,error,*999)
    IF(variableIdx<=0.OR.variableIdx>SIZE(interfaceEquationsInterpolation%variableInterpolation,1)) THEN
      localError="The specified variable index of "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
        & " is invalid. The index must be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(interfaceEquationsInterpolation%variableInterpolation,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    variableInterpolation=>interfaceEquationsInterpolation%variableInterpolation(variableIdx)
    IF(.NOT.ASSOCIATED(variableInterpolation)) THEN
      localError="The interface equations interpolation variable interpolation is not associated for variable index "// &
        & TRIM(NumberToVString(variableIdx,"*",err,error))//"."
      CALL FlagError("Interface equations interpolation interface interpolation is not associated.",err,error,*999)
    ENDIF
       
    EXITS("InterfaceEquationsInterface_VariableInterpGet")
    RETURN
999 NULLIFY(variableInterpolation)
998 ERRORS("InterfaceEquationsInterpolation_VariableInterpGet",err,error)
    EXITS("InterfaceEquationsInterpolation_VariableInterpGet")
    RETURN 1
    
  END SUBROUTINE InterfaceEquationsInterpolation_VariableInterpGet

  !
  !================================================================================================================================
  !

  !>Gets the interpolation parameters for an interface interpolation set.
  SUBROUTINE InterfaceInterpolationSet_InterpolationParametersGet(interfaceInterpolationSet,variableType,interpolationParameters, &
    & err,error,*)

    !Argument variables
    TYPE(InterfaceEquationsInterpolationSetType), POINTER :: interfaceInterpolationSet !<A pointer to the interface interpolation set to get the interpolation parameters for
    INTEGER(INTG), INTENT(IN) :: variableType !<The variable type to get the interpolation parameters for.
    TYPE(FieldInterpolationParametersType), POINTER :: interpolationParameters !<On exit, a pointer to the interpolation parameters for the specified interface interpolation set. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("InterfaceInterpolationSet_InterpolationParametersGet",err,error,*998)

    IF(ASSOCIATED(interpolationParameters)) CALL FlagError("Interpolation parameters is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceInterpolationSet)) CALL FlagError("Interface interpolation set is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(interfaceInterpolationSet%interpolationParameters)) &
      & CALL FlagError("Interface interpolation set interpolation parameters is not associated.",err,error,*999)
    IF(variableType<=0.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type should be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    interpolationParameters=>interfaceInterpolationSet%interpolationParameters(variableType)%ptr
    IF(.NOT.ASSOCIATED(interpolationParameters)) &
      & CALL FlagError("Interface interpolation set interpolation parameters is not associated.",err,error,*999)
       
    EXITS("InterfaceInterpolationSet_InterpolationParametersGet")
    RETURN
999 NULLIFY(interpolationParameters)
998 ERRORS("InterfaceInterpolationSet_InterpolationParametersGet",err,error)
    EXITS("InterfaceInterpolationSet_InterpolationParametersGet")
    RETURN 1
    
  END SUBROUTINE InterfaceInterpolationSet_InterpolationParametersGet

  !
  !================================================================================================================================
  !

  !>Gets the interpolated point for an interface interpolation set.
  SUBROUTINE InterfaceInterpolationSet_InterpolatedPointGet(interfaceInterpolationSet,variableType,interpolatedPoint, &
    & err,error,*)

    !Argument variables
    TYPE(InterfaceEquationsInterpolationSetType), POINTER :: interfaceInterpolationSet !<A pointer to the interface interpolation set to get the interpolated point for
    INTEGER(INTG), INTENT(IN) :: variableType !<The variable type to get the interpolated point for.
    TYPE(FieldInterpolatedPointType), POINTER :: interpolatedPoint !<On exit, a pointer to the interpolated point for the specified interface interpolation set. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("InterfaceInterpolationSet_InterpolatedPointGet",err,error,*998)

    IF(ASSOCIATED(interpolatedPoint)) CALL FlagError("Interpolated point is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceInterpolationSet)) CALL FlagError("Interface interpolation set is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(interfaceInterpolationSet%interpolatedPoint)) &
      & CALL FlagError("Interface interpolation set interpolated point is not associated.",err,error,*999)
    IF(variableType<=0.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type should be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    interpolatedPoint=>interfaceInterpolationSet%interpolatedPoint(variableType)%ptr
    IF(.NOT.ASSOCIATED(interpolatedPoint)) &
      & CALL FlagError("Interface interpolation set interpolated point is not associated.",err,error,*999)
       
    EXITS("InterfaceInterpolationSet_InterpolatedPointGet")
    RETURN
999 NULLIFY(interpolatedPoint)
998 ERRORS("InterfaceInterpolationSet_InterpolatedPointGet",err,error)
    EXITS("InterfaceInterpolationSet_InterpolatedPointGet")
    RETURN 1
    
  END SUBROUTINE InterfaceInterpolationSet_InterpolatedPointGet

  !
  !================================================================================================================================
  !

  !>Gets the interpolated point metrics for an interface interpolation set.
  SUBROUTINE InterfaceInterpolationSet_InterpolatedPointMetricsGet(interfaceInterpolationSet,variableType, &
    & interpolatedPointMetrics,err,error,*)

    !Argument variables
    TYPE(InterfaceEquationsInterpolationSetType), POINTER :: interfaceInterpolationSet !<A pointer to the interface interpolation set to get the interpolated point metrics for
    INTEGER(INTG), INTENT(IN) :: variableType !<The variable type to get the interpolated point metrics for.
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: interpolatedPointMetrics !<On exit, a pointer to the interpolated point metrics for the specified interface interpolation set. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("InterfaceInterpolationSet_InterpolatedPointMetricsGet",err,error,*998)

    IF(ASSOCIATED(interpolatedPointMetrics)) CALL FlagError("Interpolated point metrics is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(interfaceInterpolationSet)) CALL FlagError("Interface interpolation set is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(interfaceInterpolationSet%interpolatedPointMetrics)) &
      & CALL FlagError("Interface interpolation set interpolated point metrics is not associated.",err,error,*999)
    IF(variableType<=0.OR.variableType>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
      localError="The specified variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
        & " is invalid. The variable type should be >= 1 and <= "// &
        & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    interpolatedPointMetrics=>interfaceInterpolationSet%interpolatedPointMetrics(variableType)%ptr
    IF(.NOT.ASSOCIATED(interpolatedPointMetrics)) &
      & CALL FlagError("Interface interpolation set interpolated point metrics is not associated.",err,error,*999)
       
    EXITS("InterfaceInterpolationSet_InterpolatedPointMetricsGet")
    RETURN
999 NULLIFY(interpolatedPointMetrics)
998 ERRORS("InterfaceInterpolationSet_InterpolatedPointMetricsGet",err,error)
    EXITS("InterfaceInterpolationSet_InterpolatedPointMetricsGet")
    RETURN 1
    
  END SUBROUTINE InterfaceInterpolationSet_InterpolatedPointMetricsGet

  !
  !================================================================================================================================
  !

END MODULE InterfaceEquationsAccessRoutines
