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
  USE EquationsSetConstants
  USE FIELD_ROUTINES
  USE FieldAccessRoutines
  USE ISO_VARYING_STRING
  USE Kinds
  USE Strings
  USE Types

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !> \addtogroup EquationsRoutines_EquationTypes EquationsRoutines::EquationTypes
  !> \brief The types of equations
  !> \see EquationsRoutines,OPENCMISS_EquationsTypes
  !>@{
  INTEGER(INTG), PARAMETER :: EQUATIONS_SCALAR_TYPE=1 !<Single scalar equation. \see EquationsRoutines_EquationTypes,EquationsRoutines
  INTEGER(INTG), PARAMETER :: EQUATIONS_VECTOR_TYPE=2 !<Vector of multiple equations. \see EquationsRoutines_EquationsTypes,EquationsRoutines
  INTEGER(INTG), PARAMETER :: EQUATIONS_FUNCTIONAL_TYPE=3 !<Vector of functional equations. \see EquationsRoutines_EquationsTypes,EquationsRoutines
  !>@}

  !> \addtogroup EquationsRoutines_EquationEqualityTypes EquationsRoutines::EquationEqualityTypes
  !> \brief The types of equality for the equations
  !> \see EquationsRoutines,OPENCMISS_EquationsEqualityTypes
  !>@{
  INTEGER(INTG), PARAMETER :: EQUATIONS_EQUALS_TYPE=1 !<The equations equal zero \see EquationsRoutines_EquationEqualityTypes,EquationsRoutines
  INTEGER(INTG), PARAMETER :: EQUATIONS_LESS_THAN_TYPE=2 !<The equations are less than zero. \see EquationsRoutines_EquationsEqualityTypes,EquationsRoutines
  INTEGER(INTG), PARAMETER :: EQUATIONS_LESS_THAN_EQUALS_TYPE=3 !<The equations are less than or equal to zero. \see EquationsRoutines_EquationsEqualityTypes,EquationsRoutines
  INTEGER(INTG), PARAMETER :: EQUATIONS_GREATER_THAN_TYPE=4 !<The equations are greater than zero. \see EquationsRoutines_EquationsEqualityTypes,EquationsRoutines
  INTEGER(INTG), PARAMETER :: EQUATIONS_GREATER_THAN_EQUALS_TYPE=5 !<The equations are greater than or equal to zero. \see EquationsRoutines_EquationsEqualityTypes,EquationsRoutines
  !>@}

  !> \addtogroup EquationsRoutines_OutputTypes EquationsRoutines::OutputTypes
  !> \brief The equations output types
  !> \see EquationsRoutines,OPENCMISS_EquationsConstants
  !>@{
  INTEGER(INTG), PARAMETER :: EQUATIONS_NO_OUTPUT=0 !<No output. \see EquationsRoutines_OutputTypes,EquationsRoutines
  INTEGER(INTG), PARAMETER :: EQUATIONS_TIMING_OUTPUT=1 !<Timing information output. \see EquationsRoutines_OutputTypes,EquationsRoutines
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRIX_OUTPUT=2 !<All below and equation matrices output. \see EquationsRoutines_OutputTypes,EquationsRoutines
  INTEGER(INTG), PARAMETER :: EQUATIONS_ELEMENT_MATRIX_OUTPUT=3 !<All below and element matrices output. \see EquationsRoutines_OutputTypes,EquationsRoutines
  INTEGER(INTG), PARAMETER :: EQUATIONS_NODAL_MATRIX_OUTPUT=4 !<All below and nodal matrices output. \see EquationsRoutines_OutputTypes,EquationsRoutines
  !>@}

  !> \addtogroup EquationsRoutines_SparsityTypes EquationsRoutines::SparsityTypes
  !> \brief Equations matrices sparsity types
  !> \see EquationsRoutines,OPENCMISS_EquationsSparsityTypes
  !>@{
  INTEGER(INTG), PARAMETER :: EQUATIONS_SPARSE_MATRICES=1 !<Use sparse matrices for the equations. \see EquationsRoutines_SparsityTypes,EquationsRoutines
  INTEGER(INTG), PARAMETER :: EQUATIONS_FULL_MATRICES=2 !<Use fully populated matrices for the equations. \see EquationsRoutines_SparsityTypes,EquationsRoutines
 !>@}
 
  !> \addtogroup EquationsRoutines_LumpingTypes EquationsRoutines::LumpingTypes
  !> \brief Equations matrices lumping types
  !> \see EquationsRoutines,OPENCMISS_EquationsLumpingTypes
  !>@{
  INTEGER(INTG), PARAMETER :: EQUATIONS_UNLUMPED_MATRICES=1 !<The equations matrices are not lumped. \see EquationsRoutines_LumpingTypes,EquationsRoutines
  INTEGER(INTG), PARAMETER :: EQUATIONS_LUMPED_MATRICES=2 !<The equations matrices are "mass" lumped. \see EquationsRoutines_LumpingTypes,EquationsRoutines
 !>@}
 
  !Module types

  !Module variables

  !Interfaces
 
  PUBLIC EQUATIONS_SCALAR_TYPE,EQUATIONS_VECTOR_TYPE,EQUATIONS_FUNCTIONAL_TYPE

  PUBLIC EQUATIONS_EQUALS_TYPE,EQUATIONS_LESS_THAN_TYPE,EQUATIONS_LESS_THAN_EQUALS_TYPE,EQUATIONS_GREATER_THAN_TYPE, &
    & EQUATIONS_GREATER_THAN_EQUALS_TYPE
  
  PUBLIC EQUATIONS_NO_OUTPUT,EQUATIONS_TIMING_OUTPUT,EQUATIONS_MATRIX_OUTPUT,EQUATIONS_ELEMENT_MATRIX_OUTPUT, &
    & EQUATIONS_NODAL_MATRIX_OUTPUT

  PUBLIC EQUATIONS_SPARSE_MATRICES,EQUATIONS_FULL_MATRICES

  PUBLIC EQUATIONS_UNLUMPED_MATRICES,EQUATIONS_LUMPED_MATRICES
  
  PUBLIC Equations_Initialise,Equations_Finalise

  PUBLIC Equations_CreateStart,Equations_CreateFinish

  PUBLIC Equations_Destroy

  PUBLIC Equations_EqualityTypeGet,Equations_EqualityTypeSet
  
  PUBLIC Equations_EquationTypeGet,Equations_EquationTypeSet

  PUBLIC Equations_LinearityTypeGet,Equations_LinearityTypeSet

  PUBLIC Equations_LumpingTypeGet,Equations_LumpingTypeSet

  PUBLIC Equations_OutputTypeGet,Equations_OutputTypeSet

  PUBLIC Equations_SparsityTypeGet,Equations_SparsityTypeSet

  PUBLIC Equations_TimeDependenceTypeGet,Equations_TimeDependenceTypeSet

  PUBLIC Equations_DerivedVariableGet

  PUBLIC Equations_NumberOfLinearMatricesGet

  PUBLIC Equations_NumberOfJacobianMatricesGet

  PUBLIC Equations_NumberOfDynamicMatricesGet

  PUBLIC Equations_LinearMatrixGet

  PUBLIC Equations_JacobianMatrixGet

  PUBLIC Equations_DynamicMatrixGet

  PUBLIC Equations_DynamicMatrixGetByType

  PUBLIC Equations_DynamicMatrixTypeGet

  PUBLIC Equations_RhsVectorGet

  PUBLIC Equations_ResidualVectorGet

  PUBLIC Equations_ResidualNumberOfVariablesGet

  PUBLIC Equations_ResidualVariablesGet

  PUBLIC Equations_SourceVectorGet

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
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to create equations for
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

  !>Gets the equality type for equations.
  SUBROUTINE Equations_EqualityTypeGet(equations,equalityType,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER :: equations !<A pointer to the equations to get the equality type for
    INTEGER(INTG), INTENT(OUT) :: equalityType !<On exit, the equality type of the equations. \see EquationsRoutines_EquationEqualityTypes,EquationsRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Equations_EqualityTypeGet",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
    IF(.NOT.equations%equationsFinished) CALL FlagError("Equations has not been finished.",err,error,*999)
    
    equalityType=equations%equalityType
       
    EXITS("Equations_EqualityTypeGet")
    RETURN
999 ERRORSEXITS("Equations_EqualityTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE Equations_EqualityTypeGet
  
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

  !>Gets the equation type for equations.
  SUBROUTINE Equations_EquationTypeGet(equations,equationType,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER :: equations !<A pointer to the equations to get the equation type for
    INTEGER(INTG), INTENT(OUT) :: equationType !<On exit, the equation type of the equations. \see EquationsRoutines_EquationTypes,EquationsRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Equations_EquationTypeGet",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
    IF(.NOT.equations%equationsFinished) CALL FlagError("Equations has not been finished.",err,error,*999)
    
    equationType=equations%equationType
       
    EXITS("Equations_EquationTypeGet")
    RETURN
999 ERRORSEXITS("Equations_EquaitonTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE Equations_EquationTypeGet
  
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
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to initialise the equations for
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
    equationsSet%equations%jacobianCalculationType=EQUATIONS_JACOBIAN_FINITE_DIFFERENCE_CALCULATED
    equationsSet%equations%jacobianFiniteDifferenceStepSize=1.0E-6_DP
       
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
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
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
    
    equations%interpolation%geometricField=>equationsSet%geometry%GEOMETRIC_FIELD
    equations%interpolation%fibreField=>equationsSet%geometry%FIBRE_FIELD
    equations%interpolation%dependentField=>equationsSet%dependent%DEPENDENT_FIELD
    IF(ASSOCIATED(equationsSet%independent)) THEN
      equations%interpolation%independentField=>equationsSet%independent%INDEPENDENT_FIELD
    ELSE
      NULLIFY(equations%interpolation%independentField)
    ENDIF
    IF(ASSOCIATED(equationsSet%materials)) THEN
      equations%interpolation%materialsField=>equationsSet%materials%MATERIALS_FIELD
    ELSE
      NULLIFY(equations%interpolation%materialsField)
    ENDIF
    IF(ASSOCIATED(equationsSet%source)) THEN
      equations%interpolation%sourceField=>equationsSet%source%SOURCE_FIELD
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
      IF(equations%interpolation%independentField%type==FIELD_GEOMETRIC_TYPE.OR. &
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

  !>Gets the linearity type for equations.
  SUBROUTINE Equations_LinearityTypeGet(equations,linearityType,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER :: equations !<A pointer to the equations to get the linearity for
    INTEGER(INTG), INTENT(OUT) :: linearityType !<On exit, the linearity type of the equations. \see EquationsSetConstants_LinearityTypes,EquationsSetConstants
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Equations_LinearityTypeGet",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
    IF(.NOT.equations%equationsFinished) CALL FlagError("Equations has not been finished.",err,error,*999)
    
    linearityType=equations%linearity
       
    EXITS("Equations_LinearityTypeGet")
    RETURN
999 ERRORSEXITS("Equations_LinearityTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE Equations_LinearityTypeGet
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the linearity type for equations.
  SUBROUTINE Equations_LinearityTypeSet(equations,linearityType,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER :: equations !<A pointer to the equations to set the linearity for
    INTEGER(INTG), INTENT(IN) :: linearityType !<The linearity type to set \see EquationsSetConstants_LinearityTypes,EquationsSetConstants
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

  !>Gets the lumping type for equations.
  SUBROUTINE Equations_LumpingTypeGet(equations,lumpingType,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER :: equations !<A pointer to the equations to get the lumping type for
    INTEGER(INTG), INTENT(OUT) :: lumpingType !<On exit, the lumping type of the equations \see EquationsRoutines_LumpingTypes,EquationsRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Equations_LumpingTypeGet",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
    IF(.NOT.equations%equationsFinished) CALL FlagError("Equations has not been finished.",err,error,*999)
      
    lumpingType=equations%lumpingType
       
    EXITS("Equations_LumpingTypeGet")
    RETURN
999 ERRORSEXITS("Equations_LumpingTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE Equations_LumpingTypeGet
  
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

  !>Gets the output type for equations.
  SUBROUTINE Equations_OutputTypeGet(equations,outputType,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER :: equations !<A pointer to the equations to get the output type for
    INTEGER(INTG), INTENT(OUT) :: outputType !<On exit, the output type of the equations. \see EquationsRoutines_OutputTypes,EquationsRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Equations_OutputTypeGet",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
    IF(.NOT.equations%equationsFinished) CALL FlagError("Equations has not been finished.",err,error,*999)
    
    outputType=equations%outputType
       
    EXITS("Equations_OutputTypeGet")
    RETURN
999 ERRORSEXITS("Equations_OutputTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE Equations_OutputTypeGet
  
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

  !>Gets the sparsity type for equations.
  SUBROUTINE Equations_SparsityTypeGet(equations,sparsityType,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER :: equations !<A pointer to the equations to get the sparsity type for
    INTEGER(INTG), INTENT(OUT) :: sparsityType !<On exit, the sparsity type of the equations. \see EquationsRoutines_SparsityTypes,EquationsRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Equations_SparsityTypeGet",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
    IF(.NOT.equations%equationsFinished) CALL FlagError("Equations has not been finished.",err,error,*999)
    
    sparsityType=equations%sparsityType
       
    EXITS("Equations_SparsityTypeGet")
    RETURN
999 ERRORSEXITS("Equations_SparsityTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE Equations_SparsityTypeGet
  
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

  !>Gets the time dependence type for equations.
  SUBROUTINE Equations_TimeDependenceTypeGet(equations,timeDependenceType,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER :: equations !<A pointer to the equations to get the time dependence type for
    INTEGER(INTG), INTENT(OUT) :: timeDependenceType !<On exit, the time dependence type of the equations. \see EquationsRoutines_TimeDependenceTypes,EquationsRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Equations_TimeDependenceTypeGet",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
    IF(.NOT.equations%equationsFinished) CALL FlagError("Equations has not been finished.",err,error,*999)
    
    timeDependenceType=equations%timeDependence
       
    EXITS("Equations_TimeDependenceTypeGet")
    RETURN
999 ERRORSEXITS("Equations_TimeDependenceTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE Equations_TimeDependenceTypeGet
  
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

  !>Gets the field variable for the derived variable type
  SUBROUTINE Equations_DerivedVariableGet(equations,derivedType,fieldVariable,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER, INTENT(IN) :: equations !<A pointer to the equations to get the field variable for.
    INTEGER(INTG), INTENT(IN) :: derivedType !<The derived value type to get the field variable for. \see EquationsSetConstants_DerivedTypes.
    TYPE(FIELD_VARIABLE_TYPE), POINTER, INTENT(INOUT) :: fieldVariable !<On return, the field variable for the derived variable type.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(FIELD_TYPE), POINTER :: derivedField
    INTEGER(INTG) :: fieldVariableType
    TYPE(VARYING_STRING) :: localError

    ENTERS("Equations_DerivedVariableGet",err,error,*999)

    NULLIFY(derivedField)

    !Check pointers
    IF(ASSOCIATED(fieldVariable)) CALL FlagError("Derived field variable is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations are not associated.",err,error,*999)
    IF(derivedType<0.OR.derivedType>EQUATIONS_SET_NUMBER_OF_TENSOR_TYPES)   THEN
      localError="The derived variable type of "//TRIM(NumberToVString(derivedType,"*",err,error))// &
        & " is invalid. It should be >= 1 and <= "//TRIM(NumberToVString(EQUATIONS_SET_NUMBER_OF_TENSOR_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    NULLIFY(equationsSet)
    CALL Equations_EquationsSetGet(equations,equationsSet,err,error,*999)
    IF(.NOT.equationsSet%EQUATIONS_SET_FINISHED) CALL FlagError("Equations set has not been finished.",err,error,*999)
    IF(.NOT.ASSOCIATED(equationsSet%derived)) CALL FlagError("Equations set derived is not associated.",err,error,*999)
    IF(.NOT.equationsSet%derived%derivedFinished) &
      & CALL FlagError("Equations set derived has not been finished.",err,error,*999)
    IF(.NOT.ASSOCIATED(equationsSet%derived%derivedField))  &
      & CALL FlagError("Equations set derived field is not associated.",err,error,*999)
 
    fieldVariableType=equationsSet%derived%variableTypes(derivedType)
    NULLIFY(fieldVariable)
    CALL Field_VariableGet(equationsSet%derived%derivedField,fieldVariableType,fieldVariable,err,error,*999)

    EXITS("Equations_DerivedVariableGet")
    RETURN
999 ERRORSEXITS("Equations_DerivedVariableGet",err,error)
    RETURN 1
  END SUBROUTINE Equations_DerivedVariableGet

  !
  !================================================================================================================================
  !

  !>Get the number of linear matrices in the equations
  SUBROUTINE Equations_NumberOfLinearMatricesGet(equations,numberOfMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER, INTENT(IN) :: equations !<The equations to get the number of linear matrices for
    INTEGER(INTG), INTENT(OUT) :: numberOfMatrices !<On return, the number of linear matrices
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    !Local variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    
    ENTERS("Equations_NumberOfLinearMatricesGet",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)

    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)

    linearMatrices=>vectorMatrices%linearMatrices
    IF(ASSOCIATED(linearMatrices)) THEN
      numberOfMatrices=linearMatrices%numberOfLinearMatrices
    ELSE
      numberOfMatrices=0
    END IF

    EXITS("Equations_NumberOfLinearMatricesGet")
    RETURN
999 ERRORSEXITS("Equations_NumberOfLinearMatricesGet",err,error)
    RETURN 1

  END SUBROUTINE Equations_NumberOfLinearMatricesGet

  !
  !================================================================================================================================
  !

  !>Get the number of Jacobian matrices in the equations
  SUBROUTINE Equations_NumberOfJacobianMatricesGet(equations,numberOfMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER, INTENT(IN) :: equations !<The equations to get the number of Jacobian matrices for
    INTEGER(INTG), INTENT(OUT) :: numberOfMatrices !<On return, the number of Jacobian matrices
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    !Local variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    
    ENTERS("Equations_NumberOfJacobianMatricesGet",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)

    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)

    nonlinearMatrices=>vectorMatrices%nonlinearMatrices
    IF(ASSOCIATED(nonlinearMatrices)) THEN
      numberOfMatrices=nonlinearMatrices%numberOfJacobians
    ELSE
      numberOfMatrices=0
    END IF
 
    EXITS("Equations_NumberOfJacobianMatricesGet")
    RETURN
999 ERRORSEXITS("Equations_NumberOfJacobianMatricesGet",err,error)
    RETURN 1

  END SUBROUTINE Equations_NumberOfJacobianMatricesGet

  !
  !================================================================================================================================
  !

  !>Get the number of dynamic matrices in the equations
  SUBROUTINE Equations_NumberOfDynamicMatricesGet(equations,numberOfMatrices,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER, INTENT(IN) :: equations !<The equations to get the number of dynamic matrices for
    INTEGER(INTG), INTENT(OUT) :: numberOfMatrices !<On return, the number of dynamic matrices
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    !Local variables
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    
    ENTERS("Equations_NumberOfDynamicMatricesGet",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)

    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)

    dynamicMatrices=>vectorMatrices%dynamicMatrices
    IF(ASSOCIATED(dynamicMatrices)) THEN
      numberOfMatrices=dynamicMatrices%numberOfDynamicMatrices
    ELSE
      numberOfMatrices=0
    END IF

    EXITS("Equations_NumberOfDynamicMatricesGet")
    RETURN
999 ERRORSEXITS("Equations_NumberOfDynamicMatricesGet",err,error)
    RETURN 1

  END SUBROUTINE Equations_NumberOfDynamicMatricesGet

  !
  !================================================================================================================================
  !

  !>Get a linear equations matrix from equations
  SUBROUTINE Equations_LinearMatrixGet(equations,matrixIndex,matrix,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER, INTENT(IN) :: equations !<The equations to get the linear matrix for
    INTEGER(INTG), INTENT(IN) :: matrixIndex !<The index of the linear matrix to get
    TYPE(DistributedMatrixType), POINTER, INTENT(INOUT) :: matrix !<On return, the linear matrix requested. Must not be associated on entry.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    !Local variables
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    
    ENTERS("Equations_LinearMatrixGet",err,error,*998)

    IF(ASSOCIATED(matrix)) CALL FlagError("The matrix is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*998)

    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(linearMatrices)
    CALL EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*999)
    NULLIFY(equationsMatrix)
    CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,matrixIndex,equationsMatrix,err,error,*999)

    matrix=>equationsMatrix%matrix
    IF(.NOT.ASSOCIATED(matrix)) CALL FlagError("Equations matrix matrix is not associated.",err,error,*999)

    EXITS("Equations_LinearMatrixGet")
    RETURN
999 NULLIFY(matrix)
998 ERRORSEXITS("Equations_LinearMatrixGet",err,error)
    RETURN 1

  END SUBROUTINE Equations_LinearMatrixGet

  !
  !================================================================================================================================
  !

  !>Get a Jacobian matrix from equations
  SUBROUTINE Equations_JacobianMatrixGet(equations,residualIndex,variableType,matrix,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER, INTENT(IN) :: equations !<The equations to get the Jacobian matrix for
    INTEGER(INTG), INTENT(IN) :: residualIndex !<The index of the residual vector to get the Jacobian matrix for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type that the residual is differentiated with respect to for this Jacobian
    TYPE(DistributedMatrixType), POINTER, INTENT(INOUT) :: matrix !<On return, the requested Jacobian matrix
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    !Local variables
    INTEGER(INTG) :: matrixIndex,variableIndex
    TYPE(EquationsJacobianType), POINTER :: jacobianMatrix
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations

    ENTERS("Equations_JacobianMatrixGet",err,error,*999)

    !Check for pointer associations
    IF(ASSOCIATED(matrix)) CALL FlagError("The matrix is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*998)
    IF(residualIndex/=1) CALL FlagError("Multiple residual vectors are not yet implemented.",err,error,*998)
 
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*998)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*998)
    NULLIFY(vectorMapping)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*998)
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*998)
    NULLIFY(nonlinearMapping)
    CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*998)

    !Find Jacobian matrix index using the nonlinear equations mapping
    matrixIndex=0
    DO variableIndex=1,nonlinearMapping%numberOfResidualVariables
      IF(nonlinearMapping%residualVariables(variableIndex)%ptr%VARIABLE_TYPE==variableType) THEN
        matrixIndex=nonlinearMapping%varToJacobianMap(variableIndex)%jacobianNumber
        EXIT
      END IF
    END DO
    IF(matrixIndex==0) THEN
      CALL FlagError("Equations do not have a Jacobian matrix for residual index "// &
        & TRIM(NumberToVstring(residualIndex,"*",err,error))//" and variable type "// &
        & TRIM(NumberToVstring(variableType,"*",err,error))//".",err,error,*999)
    END IF

    !Now get Jacobian matrix using the matrix index
    NULLIFY(jacobianMatrix)
    CALL EquationsMatricesNonlinear_JacobianMatrixGet(nonlinearMatrices,matrixIndex,jacobianMatrix,err,error,*999)

    matrix=>jacobianMatrix%jacobian
    IF(.NOT.ASSOCIATED(matrix)) CALL FlagError("Jacobian matrix matrix is not associated.",err,error,*999)

    EXITS("Equations_JacobianMatrixGet")
    RETURN
999 NULLIFY(matrix)
998 ERRORSEXITS("Equations_JacobianMatrixGet",err,error)
    RETURN 1

  END SUBROUTINE Equations_JacobianMatrixGet

  !
  !================================================================================================================================
  !

  !>Get a dynamic equations matrix from equations using the dynamic matrix index
  SUBROUTINE Equations_DynamicMatrixGet(equations,matrixIndex,matrix,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER, INTENT(IN) :: equations !<The equations to get the dynamic matrix for
    INTEGER(INTG), INTENT(IN) :: matrixIndex !<The number of the dynamic matrix to get
    TYPE(DistributedMatrixType), POINTER, INTENT(INOUT) :: matrix !<On return, the requested dynamic matrix. Must not be associated on entry.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    !Local variables
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
 
    ENTERS("Equations_DynamicMatrixGet",err,error,*999)

    IF(ASSOCIATED(matrix)) CALL FlagError("The matrix is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*998)

    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(dynamicMatrices)
    CALL EquationsMatricesVector_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*999)
    NULLIFY(equationsMatrix)
    CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,matrixIndex,equationsMatrix,err,error,*999)

    matrix=>equationsMatrix%matrix
    IF(.NOT.ASSOCIATED(matrix)) CALL FlagError("Equations matrix matrix is not associated.",err,error,*999)
 
    EXITS("Equations_DynamicMatrixGet")
    RETURN
999 NULLIFY(matrix)
998 ERRORSEXITS("Equations_DynamicMatrixGet",err,error)
    RETURN 1

  END SUBROUTINE Equations_DynamicMatrixGet

  !
  !================================================================================================================================
  !

  !>Get a dynamic equations matrix from equations using the dynamic matrix type
  SUBROUTINE Equations_DynamicMatrixGetByType(equations,matrixType,matrix,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER, INTENT(IN) :: equations !<The equations to get the dynamic matrix for
    INTEGER(INTG), INTENT(IN) :: matrixType !<The type of the dynamic matrix to get. \see EquationsSetConstants_DynamicMatrixTypes,EquationsSetConstants
    TYPE(DistributedMatrixType), POINTER, INTENT(INOUT) :: matrix !<On return, the requested dynamic matrix
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    !Local variables
    INTEGER(INTG) :: matrixIndex
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(VARYING_STRING) :: localError

    ENTERS("Equations_DynamicMatrixGetByType",err,error,*998)

    !Check all pointer associations
    IF(ASSOCIATED(matrix)) CALL FlagError("The matrix is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*998)

    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
    NULLIFY(dynamicMatrices)
    CALL EquationsMatricesVector_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*999)
    NULLIFY(dynamicMapping)
    CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
 
    !Find matrix index using the equations mapping
    SELECT CASE(matrixType)
    CASE(EQUATIONS_MATRIX_STIFFNESS)
      matrixIndex=dynamicMapping%stiffnessMatrixNumber
    CASE(EQUATIONS_MATRIX_DAMPING)
      matrixIndex=dynamicMapping%dampingMatrixNumber
    CASE(EQUATIONS_MATRIX_MASS)
      matrixIndex=dynamicMapping%massMatrixNumber
    CASE DEFAULT
      localError="The specified dynamic matrix type of "//TRIM(NumberToVstring(matrixType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    IF(matrixIndex==0) THEN
      localError="The equations dynamic matrices do not have a matrix with the specified type of "// &
        & TRIM(NumberToVstring(matrixType,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
     
    NULLIFY(equationsMatrix)
    CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,matrixIndex,equationsMatrix,err,error,*999)

    matrix=>equationsMatrix%matrix
    IF(.NOT.ASSOCIATED(matrix)) CALL FlagError("Equations matrix matrix is not associated.",err,error,*999)

    EXITS("Equations_DynamicMatrixGetByType")
    RETURN
999 NULLIFY(matrix)
998 ERRORSEXITS("Equations_DynamicMatrixGetByType",err,error)
    RETURN 1

  END SUBROUTINE Equations_DynamicMatrixGetByType

  !
  !================================================================================================================================
  !

  !>Get the type of a dynamic matrix, eg. stiffness, damping or mass
  SUBROUTINE Equations_DynamicMatrixTypeGet(equations,matrixIndex,matrixType,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER, INTENT(IN) :: equations !<The equations to get the dynamic matrix for
    INTEGER(INTG), INTENT(IN) :: matrixIndex !<The number of the dynamic matrix to get
    INTEGER(INTG), INTENT(INOUT) :: matrixType !<On return, the type of the dynamic matrix. \see EquationsMatricesRoutines_DynamicMatrixTypes,EquationsMatricesRoutines
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    !Local variables
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(VARYING_STRING) :: localError

    ENTERS("Equations_DynamicMatrixTypeGet",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
    
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
    NULLIFY(dynamicMapping)
    CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)

    IF(matrixIndex<1.OR.matrixIndex>dynamicMapping%numberOfDynamicMatrices) THEN
      localError="The specified matrix index of "//TRIM(NumberToVString(matrixIndex,"*",err,error))// &
        & " is invalid. The matrix index must be >= 1 and <= "// &
        & TRIM(NumberToVstring(dynamicMapping%numberOfDynamicMatrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    IF(matrixIndex==dynamicMapping%stiffnessMatrixNumber) THEN
      matrixType=EQUATIONS_MATRIX_STIFFNESS
    ELSE IF(matrixIndex==dynamicMapping%dampingMatrixNumber) THEN
      matrixType=EQUATIONS_MATRIX_DAMPING
    ELSE IF(matrixIndex==dynamicMapping%massMatrixNumber) THEN
      matrixType=EQUATIONS_MATRIX_MASS
    ELSE
      localError="Could not find the dynamic matrix type for matrix index "//TRIM(NumberToVString(matrixIndex,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    END IF

    EXITS("Equations_DynamicMatrixTypeGet")
    RETURN
999 ERRORSEXITS("Equations_DynamicMatrixTypeGet",err,error)
    RETURN 1

  END SUBROUTINE Equations_DynamicMatrixTypeGet

  !
  !================================================================================================================================
  !

  !>Get the right hand side vector for equations
  SUBROUTINE Equations_RhsVectorGet(equations,vector,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER, INTENT(IN) :: equations !<The equations to get the right hand side vector for
    TYPE(DistributedVectorType), POINTER, INTENT(INOUT) :: vector !<On return, the right hand side vector for the equations. Must not be associated on entry.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    !Local variables
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
 
    ENTERS("Equations_RhsVectorGet",err,error,*998)

    IF(ASSOCIATED(vector)) CALL FlagError("Vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations are not associated.",err,error,*998)

    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*998)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*998)
    NULLIFY(rhsVector)
    CALL EquationsMatricesVector_RHSVectorGet(vectorMatrices,rhsVector,err,error,*998)
    
    vector=>rhsVector%vector
    IF(.NOT.ASSOCIATED(vector)) CALL FlagError("The RHS vector vector is not associated.",err,error,*999)

    EXITS("Equations_RhsVectorGet")
    RETURN
999 NULLIFY(vector)
998 ERRORSEXITS("Equations_RhsVectorGet",err,error)
    RETURN 1

  END SUBROUTINE Equations_RhsVectorGet

  !
  !================================================================================================================================
  !

  !>Get a residual vector for nonlinear equations
  SUBROUTINE Equations_ResidualVectorGet(equations,residualIndex,vector,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER, INTENT(IN) :: equations !<The equations to get the residual vector for
    INTEGER(INTG), INTENT(IN) :: residualIndex !<The index of the residual vector to get
    TYPE(DistributedVectorType), POINTER, INTENT(INOUT) :: vector !<On return, the residual vector for the equations
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    !Local variables
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations

    ENTERS("Equations_ResidualVectorGet",err,error,*998)

    IF(ASSOCIATED(vector)) CALL FlagError("Vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations are not associated.",err,error,*998)
    IF(residualIndex/=1) CALL FlagError("Multiple residual vectors are not yet implemented.",err,error,*998)
   
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*998)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*998)
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*998)
 
    vector=>nonlinearMatrices%residual
    IF(.NOT.ASSOCIATED(vector)) CALL FlagError("The residual vector vector is not associated.",err,error,*999)

    EXITS("Equations_ResidualVectorGet")
    RETURN
999 NULLIFY(vector)
998 ERRORSEXITS("Equations_ResidualVectorGet",err,error)
    RETURN 1

  END SUBROUTINE Equations_ResidualVectorGet

  !
  !================================================================================================================================
  !

  !>Get the number of field variables that contribute to the residual vector
  SUBROUTINE Equations_ResidualNumberOfVariablesGet(equations,residualIndex,numberOfVariables,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER, INTENT(IN) :: equations !<The equations to get the residual vector number of variables for
    INTEGER(INTG), INTENT(IN) :: residualIndex !<The index of the residual vector to get the number of variables for
    INTEGER(INTG), INTENT(OUT) :: numberOfVariables !<On return, the number of variables that contribute to the residual vector
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    !Local variables
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsVectorType), POINTER :: vectorEquations

    ENTERS("Equations_ResidualNumberOfVariablesGet",err,error,*999)

    !Check for pointer associations
    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations are not associated.",err,error,*999)
    IF(residualIndex/=1) CALL FlagError("Multiple residual vectors are not yet implemented.",err,error,*999)
    
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
    NULLIFY(nonlinearMapping)
    CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
    
    numberOfVariables=nonlinearMapping%numberOfResidualVariables

    EXITS("Equations_ResidualNumberOfVariablesGet")
    RETURN
999 ERRORSEXITS("Equations_ResidualNumberOfVariablesGet",err,error)
    RETURN 1

  END SUBROUTINE Equations_ResidualNumberOfVariablesGet

  !
  !================================================================================================================================
  !

  !>Get the field variables that contribute to the residual vector
  SUBROUTINE Equations_ResidualVariablesGet(equations,residualIndex,residualVariables,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER, INTENT(IN) :: equations !<The equations to get the residual vector variables for
    INTEGER(INTG), INTENT(IN) :: residualIndex !<The index of the residual vector to get the variables for
    INTEGER(INTG), INTENT(OUT) :: residualVariables(:) !<residualVariables(varIdx). On return, the field variable type for the varIdx'th residual variable
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    !Local variables
    INTEGER(INTG) :: numberOfVariables,variableIdx
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(VARYING_STRING) :: localError

    ENTERS("Equations_ResidualVariablesGet",err,error,*999)

    !Check for pointer associations
    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations are not associated.",err,error,*999)
    IF(residualIndex/=1) CALL FlagError("Multiple residual vectors are not yet implemented.",err,error,*999)

    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
    NULLIFY(nonlinearMapping)
    CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)

    numberOfVariables=nonlinearMapping%numberOfResidualVariables
    IF(SIZE(residualVariables,1)<numberOfVariables) THEN
      localError="The size of the specified residual variables array is too small. The array must have a size of >= "// &
        & TRIM(numberToVstring(numberOfVariables,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
      
    DO variableIdx=1,numberOfVariables
      IF(.NOT.ASSOCIATED(nonlinearMapping%residualVariables(variableIdx)%ptr)) THEN
        localError="The residual variable is not associated for variable index "// &
          & TRIM(NumberToVString(variableIdx,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      residualVariables(variableIdx)=nonlinearMapping%residualVariables(variableIdx)%ptr%VARIABLE_TYPE
    END DO
 
    EXITS("Equations_ResidualVariablesGet")
    RETURN
999 ERRORSEXITS("Equations_ResidualVariablesGet",err,error)
    RETURN 1

  END SUBROUTINE Equations_ResidualVariablesGet

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

  !>Get the source vector for equations
  SUBROUTINE Equations_SourceVectorGet(equations,vector,err,error,*)

    !Argument variables
    TYPE(EquationsType), POINTER, INTENT(IN) :: equations !<The equations to get the source vector for
    TYPE(DistributedVectorType), POINTER, INTENT(INOUT) :: vector !<On return, the source vector for the equations. Must not be associated on entry.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    !Local variables
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations

    ENTERS("Equations_SourceVectorGet",err,error,*999)

    IF(ASSOCIATED(vector)) CALL FlagError("Vector is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations are not associated.",err,error,*998)

    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*998)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*998)
    NULLIFY(sourceVector)
    CALL EquationsMatricesVector_SourceVectorGet(vectorMatrices,sourceVector,err,error,*998)
    
    vector=>sourceVector%vector
    IF(.NOT.ASSOCIATED(vector)) CALL FlagError("The source vector vector is not associated.",err,error,*999)

    EXITS("Equations_SourceVectorGet")
    RETURN
999 NULLIFY(vector)
998 ERRORSEXITS("Equations_SourceVectorGet",err,error)
    RETURN 1

  END SUBROUTINE Equations_SourceVectorGet

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
