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

  USE BASE_ROUTINES
  USE EQUATIONS_MAPPING_ROUTINES
  USE EQUATIONS_MATRICES_ROUTINES
  USE EQUATIONS_SET_CONSTANTS
  USE FIELD_ROUTINES
  USE FieldAccessRoutines
  USE ISO_VARYING_STRING
  USE KINDS
  USE STRINGS
  USE TYPES

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

  INTERFACE EQUATIONS_CREATE_START
    MODULE PROCEDURE Equations_CreateStart
  END INTERFACE EQUATIONS_CREATE_START

  INTERFACE Equations_CreateStart
    MODULE PROCEDURE EQUATIONS_CREATE_START
  END INTERFACE Equations_CreateStart

  INTERFACE Equations_CreateFinish
    MODULE PROCEDURE EQUATIONS_CREATE_FINISH
  END INTERFACE Equations_CreateFinish

  INTERFACE Equations_LinearityTypeGet
    MODULE PROCEDURE EQUATIONS_LINEARITY_TYPE_GET
  END INTERFACE Equations_LinearityTypeGet
 
  INTERFACE Equations_LinearityTypeSet
    MODULE PROCEDURE EQUATIONS_LINEARITY_TYPE_SET
  END INTERFACE Equations_LinearityTypeSet

  INTERFACE Equations_LumpingTypeGet
    MODULE PROCEDURE EQUATIONS_LUMPING_TYPE_GET
  END INTERFACE Equations_LumpingTypeGet
 
  INTERFACE Equations_LumpingTypeSet
    MODULE PROCEDURE EQUATIONS_LUMPING_TYPE_SET
  END INTERFACE Equations_LumpingTypeSet
 
  INTERFACE Equations_OutputTypeGet
    MODULE PROCEDURE EQUATIONS_OUTPUT_TYPE_GET
  END INTERFACE Equations_OutputTypeGet
 
  INTERFACE Equations_OutputTypeSet
    MODULE PROCEDURE EQUATIONS_OUTPUT_TYPE_SET
  END INTERFACE Equations_OutputTypeSet

  INTERFACE Equations_SparsityTypeGet
    MODULE PROCEDURE EQUATIONS_SPARSITY_TYPE_GET
  END INTERFACE Equations_SparsityTypeGet

  INTERFACE Equations_SparsityTypeSet
    MODULE PROCEDURE EQUATIONS_SPARSITY_TYPE_SET
  END INTERFACE Equations_SparsityTypeSet
 
  INTERFACE Equations_TimeDependenceTypeGet
    MODULE PROCEDURE EQUATIONS_TIME_DEPENDENCE_TYPE_GET
  END INTERFACE Equations_TimeDependenceTypeGet

  INTERFACE Equations_TimeDependenceTypeSet
    MODULE PROCEDURE EQUATIONS_TIME_DEPENDENCE_TYPE_SET
  END INTERFACE Equations_TimeDependenceTypeSet

  PUBLIC EQUATIONS_SCALAR_TYPE,EQUATIONS_VECTOR_TYPE,EQUATIONS_FUNCTIONAL_TYPE

  PUBLIC EQUATIONS_EQUALS_TYPE,EQUATIONS_LESS_THAN_TYPE,EQUATIONS_LESS_THAN_EQUALS_TYPE,EQUATIONS_GREATER_THAN_TYPE, &
    & EQUATIONS_GREATER_THAN_EQUALS_TYPE
  PUBLIC EQUATIONS_NO_OUTPUT,EQUATIONS_TIMING_OUTPUT,EQUATIONS_MATRIX_OUTPUT,EQUATIONS_ELEMENT_MATRIX_OUTPUT

  PUBLIC EQUATIONS_NODAL_MATRIX_OUTPUT

  PUBLIC EQUATIONS_SPARSE_MATRICES,EQUATIONS_FULL_MATRICES

  PUBLIC EQUATIONS_UNLUMPED_MATRICES,EQUATIONS_LUMPED_MATRICES
  
  PUBLIC EQUATIONS_CREATE_START,EQUATIONS_CREATE_FINISH

  PUBLIC Equations_CreateStart,Equations_CreateFinish

  PUBLIC EQUATIONS_DESTROY

  PUBLIC Equations_Destroy

  PUBLIC EQUATIONS_INITIALISE,EQUATIONS_FINALISE

  PUBLIC Equations_Initialise,Equations_Finalise

  PUBLIC EQUATIONS_LINEARITY_TYPE_GET,EQUATIONS_LINEARITY_TYPE_SET

  PUBLIC Equations_LinearityTypeGet,Equations_LinearityTypeSet

  PUBLIC EQUATIONS_LUMPING_TYPE_GET,EQUATIONS_LUMPING_TYPE_SET

  PUBLIC Equations_LumpingTypeGet,Equations_LumpingTypeSet

  PUBLIC EQUATIONS_OUTPUT_TYPE_GET,EQUATIONS_OUTPUT_TYPE_SET

  PUBLIC Equations_OutputTypeGet,Equations_OutputTypeSet

  PUBLIC EQUATIONS_SPARSITY_TYPE_GET,EQUATIONS_SPARSITY_TYPE_SET

  PUBLIC Equations_SparsityTypeGet,Equations_SparsityTypeSet

  PUBLIC EQUATIONS_TIME_DEPENDENCE_TYPE_GET,EQUATIONS_TIME_DEPENDENCE_TYPE_SET

  PUBLIC Equations_TimeDependenceTypeGet,Equations_TimeDependenceTypeSet

  PUBLIC EQUATIONS_SET_EQUATIONS_GET

  PUBLIC EquationsSet_EquationsGet

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
    TYPE(EQUATIONS_TYPE), POINTER :: equations !<A pointer to the equations to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Equations_CreateFinish",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
    IF(equations%EQUATIONS_FINISHED) CALL FlagError("Equations have already been finished.",err,error,*999)
    
    !Initialise the eqations interpolation
    CALL Equations_InterpolationInitialise(equations,err,error,*999)
    !Set the finished flag
    equations%EQUATIONS_FINISHED=.TRUE.
       
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
    TYPE(EQUATIONS_TYPE), POINTER :: equations !<On exit, a pointer to the created equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Equations_CreateStart",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(ASSOCIATED(equationsSet%equations)) CALL FlagError("Equations are already associated for the equations set.",err,error,*999)
    IF(ASSOCIATED(equations)) CALL FlagError("Equations is already associated.",err,error,*999)

    !Initialise the equations
    CALL Equations_Initialise(equationsSet,err,error,*999)
    !Return the pointer
    equations=>equationsSet%equations
       
    EXITS("Equations_CreateStart")
    RETURN
999 ERRORSEXITS("Equations_CreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE Equations_CreateStart

  !
  !================================================================================================================================
  !

  !>Destroys equations
  SUBROUTINE Equations_Destroy(equations,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER :: equations !<A pointer to the equations to destroy
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
    TYPE(EQUATIONS_TYPE), POINTER :: equations !<A pointer to the equations to get the equality type for
    INTEGER(INTG), INTENT(OUT) :: equalityType !<On exit, the equality type of the equations. \see EquationsRoutines_EquationEqualityTypes,EquationsRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Equations_EqualityTypeGet",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
    IF(.NOT.equations%EQUATIONS_FINISHED) CALL FlagError("Equations has not been finished.",err,error,*999)
    
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
    TYPE(EQUATIONS_TYPE), POINTER :: equations !<A pointer to the equations to set the equality for
    INTEGER(INTG), INTENT(IN) :: equalityType !<The equality type to set \see EquationsRoutines_EquationEqualityTypes,EquationsRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Equations_EqualityTypeSet",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
    IF(equations%EQUATIONS_FINISHED) CALL FlagError("Equations has already been finished.",err,error,*999)

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
    TYPE(EQUATIONS_TYPE), POINTER :: equations !<A pointer to the equations to get the equation type for
    INTEGER(INTG), INTENT(OUT) :: equationType !<On exit, the equation type of the equations. \see EquationsRoutines_EquationTypes,EquationsRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Equations_EquationTypeGet",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
    IF(.NOT.equations%EQUATIONS_FINISHED) CALL FlagError("Equations has not been finished.",err,error,*999)
    
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
    TYPE(EQUATIONS_TYPE), POINTER :: equations !<A pointer to the equations to set the equation for
    INTEGER(INTG), INTENT(IN) :: equationType !<The equation type to set \see EquationsRoutines_EquationsTypes,EquationsRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Equations_EquationTypeSet",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
    IF(equations%EQUATIONS_FINISHED) CALL FlagError("Equations has already been finished.",err,error,*999)

    SELECT CASE(equationType)
    CASE(EQUATIONS_SCALAR_TYPE)
      equations%equationsType=EQUATIONS_SCALAR_TYPE
    CASE(EQUATIONS_VECTOR_TYPE)
      equations%equationsType=EQUATIONS_VECTOR_TYPE
    CASE(EQUATIONS_FUNCTIONAL_TYPE)
      equations%equationsType=EQUATIONS_FUNCTIONAL_TYPE
    CASE DEFAULT
      localError="The specified equation type of "//TRIM(NumberToVString(equationType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("Equations_EquationTypeSet")
    RETURN
999 ERRORSEXITS("Equations_EquationTypeSet",err,error)
    RETURN 1
  END SUBROUTINE Equations_EquationTypeSet
  
  !
  !================================================================================================================================
  !

  !>Finalise the equations and deallocate all memory.
  SUBROUTINE Equations_Finalise(equations,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER :: equations !<A pointer to the equations to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Equations_Finalise",err,error,*999)

    IF(ASSOCIATED(equations)) THEN
      CALL Equations_InterpolationFinalise(equations%INTERPOLATION,err,error,*999)
      IF(ASSOCIATED(equations%EQUATIONS_MAPPING)) CALL EQUATIONS_MAPPING_DESTROY(equations%EQUATIONS_MAPPING,err,error,*999)
      IF(ASSOCIATED(equations%EQUATIONS_MATRICES)) CALL EQUATIONS_MATRICES_DESTROY(equations%EQUATIONS_MATRICES,err,error,*999)
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
    equationsSet%equations%EQUATIONS_SET=>EQUATIONS_SET
    equationsSet%equations%LINEARITY=EQUATIONS_LINEAR
    equationsSet%equations%timeDependence=EQUATIONS_STATIC
    equationsSet%equations%OUTPUT_TYPE=EQUATIONS_NO_OUTPUT
    equationsSet%equations%SPARSITY_TYPE=EQUATIONS_SPARSE_MATRICES
    equationsSet%equations%LUMPING_TYPE=EQUATIONS_UNLUMPED_MATRICES
    NULLIFY(equationsSet%equations%INTERPOLATION)
    NULLIFY(equationsSet%equations%EQUATIONS_MAPPING)
    NULLIFY(equationsSet%equations%EQUATIONS_MATRICES)
    equationsSet%equations%EQUATIONS_FINISHED=.FALSE.
       
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
    TYPE(EQUATIONS_INTERPOLATION_TYPE), POINTER :: equationsInterpolation !<A pointer to the equations interpolation to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Equations_InterpolationFinalise",err,error,*999)

    IF(ASSOCIATED(equationsInterpolation)) THEN
      CALL Field_InterpolationParametersFinalise(equationsInterpolation%GEOMETRIC_INTERP_PARAMETERS,err,error,*999)
      CALL Field_InterpolationParametersFinalise(equationsInterpolation%FIBRE_INTERP_PARAMETERS,err,error,*999)
      CALL Field_InterpolationParametersFinalise(equationsInterpolation%DEPENDENT_INTERP_PARAMETERS,err,error,*999)
      CALL Field_InterpolationParametersFinalise(equationsInterpolation%prevDependentInterpParameters,err,error,*999)
      CALL Field_InterpolationParametersFinalise(equationsInterpolation%INDEPENDENT_INTERP_PARAMETERS,err,error,*999)
      CALL Field_InterpolationParametersFinalise(equationsInterpolation%MATERIALS_INTERP_PARAMETERS,err,error,*999)
      CALL Field_InterpolationParametersFinalise(equationsInterpolation%SOURCE_INTERP_PARAMETERS,err,error,*999)
      CALL Field_InterpolatedPointsFinalise(equationsInterpolation%GEOMETRIC_INTERP_POINT,err,error,*999)
      CALL Field_InterpolatedPointsFinalise(equationsInterpolation%DEPENDENT_INTERP_POINT,err,error,*999)
      CALL Field_InterpolatedPointsFinalise(equationsInterpolation%prevDependentInterpPoint,err,error,*999)
      CALL Field_InterpolatedPointsFinalise(equationsInterpolation%INDEPENDENT_INTERP_POINT,err,error,*999)
      CALL Field_InterpolatedPointsFinalise(equationsInterpolation%FIBRE_INTERP_POINT,err,error,*999)
      CALL Field_InterpolatedPointsFinalise(equationsInterpolation%MATERIALS_INTERP_POINT,err,error,*999)
      CALL Field_InterpolatedPointsFinalise(equationsInterpolation%SOURCE_INTERP_POINT,err,error,*999)
      CALL Field_PhysicalPointsFinalise(equationsInterpolation%DEPENDENT_PHYSICAL_POINT,err,error,*999)
      CALL Field_InterpolatedPointsMetricsFinalise(equationsInterpolation%DEPENDENT_INTERP_POINT_METRICS,err,error,*999)
      CALL Field_InterpolatedPointsMetricsFinalise(equationsInterpolation%prevDependentInterpPointMetrics,err,error,*999)
      CALL Field_InterpolatedPointsMetricsFinalise(equationsInterpolation%INDEPENDENT_INTERP_POINT_METRICS,err,error,*999)
      CALL Field_InterpolatedPointsMetricsFinalise(equationsInterpolation%GEOMETRIC_INTERP_POINT_METRICS,err,error,*999)
      CALL Field_InterpolatedPointsMetricsFinalise(equationsInterpolation%FIBRE_INTERP_POINT_METRICS,err,error,*999)
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
    TYPE(EQUATIONS_TYPE), POINTER :: equations !<The pointer to the equations to initialise the interpolation for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("Equations_InterpolationInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated",err,error,*998)
    IF(.NOT.ASSOCIATED(equations%EQUATIONS_SET)) CALL FlagError("Equations equation set is not associated",err,error,*998)THEN
    IF(ASSOCIATED(equations%interpolation)) &
      & CALL FlagError("Interpolation is already associated for these equations.",err,error,*998)

    equationsSet=>equations%EQUATIONS_SET
    ALLOCATE(equations%interpolation,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate equations interpolation",err,error,*999)
    equations%interpolation%equations=>equations
    NULLIFY(equations%interpolation%GEOMETRIC_INTERP_PARAMETERS)
    NULLIFY(equations%interpolation%FIBRE_INTERP_PARAMETERS)
    NULLIFY(equations%interpolation%DEPENDENT_INTERP_PARAMETERS)
    NULLIFY(equations%interpolation%prevDependentInterpParameters)
    NULLIFY(equations%interpolation%INDEPENDENT_INTERP_PARAMETERS)
    NULLIFY(equations%interpolation%MATERIALS_INTERP_PARAMETERS)
    NULLIFY(equations%interpolation%SOURCE_INTERP_PARAMETERS)
    NULLIFY(equations%interpolation%GEOMETRIC_INTERP_POINT)
    NULLIFY(equations%interpolation%FIBRE_INTERP_POINT)
    NULLIFY(equations%interpolation%DEPENDENT_INTERP_POINT)
    NULLIFY(equations%interpolation%prevDependentInterpPoint)
    NULLIFY(equations%interpolation%INDEPENDENT_INTERP_POINT)
    NULLIFY(equations%interpolation%MATERIALS_INTERP_POINT)
    NULLIFY(equations%interpolation%SOURCE_INTERP_POINT)
    NULLIFY(equations%interpolation%DEPENDENT_PHYSICAL_POINT)
    NULLIFY(equations%interpolation%DEPENDENT_INTERP_POINT_METRICS)
    NULLIFY(equations%interpolation%prevDependentInterpPointMetrics)
    NULLIFY(equations%interpolation%INDEPENDENT_INTERP_POINT_METRICS)
    NULLIFY(equations%interpolation%GEOMETRIC_INTERP_POINT_METRICS)
    NULLIFY(equations%interpolation%FIBRE_INTERP_POINT_METRICS)
    
    equations%interpolation%GEOMETRIC_FIELD=>equationsSet%geometry%GEOMETRIC_FIELD
    equations%interpolation%FIBRE_FIELD=>equationsSet%geometry%FIBRE_FIELD
    equations%interpolation%DEPENDENT_FIELD=>equationsSet%dependent%DEPENDENT_FIELD
    IF(ASSOCIATED(equationsSet%INDEPENDENT)) THEN
      equations%interpolation%INDEPENDENT_FIELD=>equationsSet%independent%INDEPENDENT_FIELD
    ELSE
      NULLIFY(equations%interpolation%INDEPENDENT_FIELD)
    ENDIF
    IF(ASSOCIATED(equationsSet%materials)) THEN
      equations%interpolation%MATERIALS_FIELD=>equationsSet%materials%MATERIALS_FIELD
    ELSE
      NULLIFY(equations%interpolation%MATERIALS_FIELD)
    ENDIF
    IF(ASSOCIATED(equationsSet%source)) THEN
      equations%interpolation%SOURCE_FIELD=>equationsSet%source%SOURCE_FIELD
    ELSE
      NULLIFY(equations%interpolation%SOURCE_FIELD)
    ENDIF
    
    CALL Field_InterpolationParametersInitialise(equations%interpolation%GEOMETRIC_FIELD, &
      & equations%interpolation%GEOMETRIC_INTERP_PARAMETERS,err,error,*999)
    CALL Field_InterpolatedPointsInitialise(equations%interpolation%GEOMETRIC_INTERP_PARAMETERS, &
      & equations%interpolation%GEOMETRIC_INTERP_POINT,err,error,*999)
    CALL Field_InterpolatedPointsMetricsInitialise(equations%interpolation%GEOMETRIC_INTERP_POINT, &
      & equations%interpolation%GEOMETRIC_INTERP_POINT_METRICS,err,error,*999)
    CALL Field_InterpolationParametersInitialise(equations%interpolation%DEPENDENT_FIELD, &
      & equations%interpolation%DEPENDENT_INTERP_PARAMETERS,err,error,*999)
    CALL Field_InterpolatedPointsInitialise(equations%interpolation%DEPENDENT_INTERP_PARAMETERS, &
      & equations%interpolation%DEPENDENT_INTERP_POINT,err,error,*999)
    !CALL Field_PhysicalPointsInitialise(equations%interpolation%DEPENDENT_INTERP_POINT, &
    !  & equations%interpolation%GEOMETRIC_INTERP_POINT,equations%interpolation%DEPENDENT_PHYSICAL_POINT, &
    !  & err,error,*999)
    IF(equations%interpolation%DEPENDENT_FIELD%TYPE==FIELD_GEOMETRIC_TYPE.OR. &
      & equations%interpolation%DEPENDENT_FIELD%TYPE==FIELD_FIBRE_TYPE.OR. &
      & equations%interpolation%DEPENDENT_FIELD%TYPE==FIELD_GEOMETRIC_GENERAL_TYPE) THEN
      CALL Field_InterpolatedPointsMetricsInitialise(equations%interpolation%DEPENDENT_INTERP_POINT, &
        & equations%interpolation%DEPENDENT_INTERP_POINT_METRICS,err,error,*999)
    ENDIF
    IF(equations%timeDependence/=EQUATIONS_STATIC) THEN
      CALL Field_InterpolationParametersInitialise(equations%interpolation%dependent_Field, &
        & equations%interpolation%prevDependentInterpParameters,err,error,*999)
      CALL Field_InterpolatedPointsInitialise(equations%interpolation%prevDependentInterpParameters, &
        & equations%interpolation%prevDependentInterpPoint,err,error,*999)
      IF(equations%interpolation%DEPENDENT_FIELD%TYPE==FIELD_GEOMETRIC_TYPE.OR. &
        & equations%interpolation%DEPENDENT_FIELD%TYPE==FIELD_FIBRE_TYPE.OR. &
        & equations%interpolation%DEPENDENT_FIELD%TYPE==FIELD_GEOMETRIC_GENERAL_TYPE) THEN
        CALL Field_InterpolatedPointsMetricsInitialise(equations%interpolation%prevDependentInterpPoint, &
          & equations%interpolation%prevDependentInterpPointMetrics,err,error,*999)
      ENDIF
    ENDIF
    IF(ASSOCIATED(equations%interpolation%FIBRE_FIELD)) THEN
      CALL Field_InterpolationParametersInitialise(equations%interpolation%FIBRE_FIELD, &
        & equations%interpolation%FIBRE_INTERP_PARAMETERS,err,error,*999)
      CALL Field_InterpolatedPointsInitialise(equations%interpolation%FIBRE_INTERP_PARAMETERS,  &
        & equations%interpolation%FIBRE_INTERP_POINT,err,error,*999)
      CALL Field_InterpolatedPointsMetricsInitialise(equations%interpolation%FIBRE_INTERP_POINT,  &
        & equations%interpolation%FIBRE_INTERP_POINT_METRICS,err,error,*999)
    ENDIF
    IF(ASSOCIATED(equations%interpolation%INDEPENDENT_FIELD)) THEN
      CALL Field_InterpolationParametersInitialise(equations%interpolation%INDEPENDENT_FIELD, &
        & equations%interpolation%INDEPENDENT_INTERP_PARAMETERS,err,error,*999)
      CALL Field_InterpolatedPointsInitialise(equations%interpolation%INDEPENDENT_INTERP_PARAMETERS,  &
        & equations%interpolation%INDEPENDENT_INTERP_POINT,err,error,*999)
      IF(equations%interpolation%INDEPENDENT_FIELD%TYPE==FIELD_GEOMETRIC_TYPE.OR. &
        & equations%interpolation%INDEPENDENT_FIELD%TYPE==FIELD_FIBRE_TYPE) THEN
        CALL Field_InterpolatedPointsMetricsInitialise(equations%interpolation%INDEPENDENT_INTERP_POINT,  &
          &  equations%interpolation%INDEPENDENT_INTERP_POINT_METRICS,err,error,*999)
      END IF
    ENDIF
    IF(ASSOCIATED(equations%interpolation%MATERIALS_FIELD)) THEN
      CALL Field_InterpolationParametersInitialise(equations%interpolation%MATERIALS_FIELD, &
        & equations%interpolation%MATERIALS_INTERP_PARAMETERS,err,error,*999)
      CALL Field_InterpolatedPointsInitialise(equations%interpolation%MATERIALS_INTERP_PARAMETERS,  &
        & equations%interpolation%MATERIALS_INTERP_POINT,err,error,*999)
    ENDIF
    IF(ASSOCIATED(equations%interpolation%SOURCE_FIELD)) THEN
      CALL Field_InterpolationParametersInitialise(equations%interpolation%SOURCE_FIELD, &
        & equations%interpolation%SOURCE_INTERP_PARAMETERS,err,error,*999)
      CALL Field_InterpolatedPointsInitialise(equations%interpolation%SOURCE_INTERP_PARAMETERS, &
        & equations%interpolation%SOURCE_INTERP_POINT,err,error,*999)
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
    TYPE(EQUATIONS_TYPE), POINTER :: equations !<A pointer to the equations to get the linearity for
    INTEGER(INTG), INTENT(OUT) :: linearityType !<On exit, the linearity type of the equations. \see EQUATIONS_SET_CONSTANTS_LinearityTypes,EQUATIONS_SET_CONSTANTS
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Equations_LinearityTypeGet",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
    IF(.NOT.equations%EQUATIONS_FINISHED) CALL FlagError("Equations has not been finished.",err,error,*999)
    
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
    TYPE(EQUATIONS_TYPE), POINTER :: equations !<A pointer to the equations to set the linearity for
    INTEGER(INTG), INTENT(IN) :: linearityType !<The linearity type to set \see EQUATIONS_SET_CONSTANTS_LinearityTypes,EQUATIONS_SET_CONSTANTS
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Equations_LinearityTypeSet",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
    IF(equations%EQUATIONS_FINISHED) CALL FlagError("Equations has already been finished.",err,error,*999)

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
    TYPE(EQUATIONS_TYPE), POINTER :: equations !<A pointer to the equations to get the lumping type for
    INTEGER(INTG), INTENT(OUT) :: lumpingType !<On exit, the lumping type of the equations \see EquationsRoutines_LumpingTypes,EquationsRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Equations_LumpingTypeGet",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
    IF(.NOT.equations%EQUATIONS_FINISHED) CALL FlagError("Equations has not been finished.",err,error,*999)
      
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
    TYPE(EQUATIONS_TYPE), POINTER :: equations !<A pointer to the equations to set the lumping for
    INTEGER(INTG), INTENT(IN) :: lumpingType !<The lumping type to set \see EquationsRoutines_LumpingTypes,EquationsRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Equations_LumpingTypeSet",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
    IF(equations%EQUATIONS_FINISHED) CALL FlagError("Equations has already been finished.",err,error,*999)

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
    TYPE(EQUATIONS_TYPE), POINTER :: equations !<A pointer to the equations to get the output type for
    INTEGER(INTG), INTENT(OUT) :: outputType !<On exit, the output type of the equations. \see EquationsRoutines_OutputTypes,EquationsRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Equations_OutputTypeGet",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
    IF(.NOT.equations%EQUATIONS_FINISHED) CALL FlagError("Equations has not been finished.",err,error,*999)
    
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
    TYPE(EQUATIONS_TYPE), POINTER :: equations !<A pointer to the equations to set the output type for
    INTEGER(INTG), INTENT(IN) :: outputType !<The output type to set \see EquationsRoutines_OutputTypes,EquationsRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Equations_OutputTypeSet",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
    IF(equations%EQUATIONS_FINISHED) CALL FlagError("Equations has already been finished.",err,error,*999)

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
    TYPE(EQUATIONS_TYPE), POINTER :: equations !<A pointer to the equations to get the sparsity type for
    INTEGER(INTG), INTENT(OUT) :: sparsityType !<On exit, the sparsity type of the equations. \see EquationsRoutines_SparsityTypes,EquationsRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Equations_SparsityTypeGet",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
    IF(.NOT.equations%EQUATIONS_FINISHED) CALL FlagError("Equations has not been finished.",err,error,*999)
    
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
    TYPE(EQUATIONS_TYPE), POINTER :: equations !<A pointer to the equations to set the sparsity type for
    INTEGER(INTG), INTENT(IN) :: sparsityType !<The sparsity type to set \see EquationsRoutines_SparsityTypes,EquationsRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Equations_SparsityTypeSet",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
    IF(equations%EQUATIONS_FINISHED) CALL FlagError("Equations has already been finished.",err,error,*999)
    
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
    TYPE(EQUATIONS_TYPE), POINTER :: equations !<A pointer to the equations to get the time dependence type for
    INTEGER(INTG), INTENT(OUT) :: timeDependenceType !<On exit, the time dependence type of the equations. \see EquationsRoutines_TimeDependenceTypes,EquationsRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Equations_TimeDependenceTypeGet",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
    IF(.NOT.equations%EQUATIONS_FINISHED) CALL FlagError("Equations has not been finished.",err,error,*999)
    
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
    TYPE(EQUATIONS_TYPE), POINTER :: equations !<A pointer to the equations to set the linearity for
    INTEGER(INTG), INTENT(IN) :: timeDependenceType !<The time dependence type to set \see EquationsRoutines_TimeDependenceTypes,EquationsRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Equations_TimeDependenceTypeSet",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
    IF(equations%EQUATIONS_FINISHED) CALL FlagError("Equations has already been finished.",err,error,*999)

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
    TYPE(EQUATIONS_TYPE), POINTER, INTENT(IN) :: equations !<A pointer to the equations to get the field variable for.
    INTEGER(INTG), INTENT(IN) :: derivedType !<The derived value type to get the field variable for. \see EQUATIONS_SET_CONSTANTS_DerivedTypes.
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
    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations are not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equations%equations_set)) CALL FlagError("Equations equations set is not associated.",err,error,*999)
    IF(ASSOCIATED(fieldVariable)) CALL FlagError("Derived field variable is already associated.",err,error,*999)
    IF(derivedType<0.OR.derivedType>EQUATIONS_SET_NUMBER_OF_DERIVED_TYPES)   THEN
      localError="The derived variable type of "//TRIM(NumberToVString(derivedType,"*",err,error))// &
        & " is invalid. It should be >= 1 and <= "//TRIM(NumberToVString(EQUATIONS_SET_NUMBER_OF_DERIVED_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.equations%equations_Set%EQUATIONS_SET_FINISHED) CALL FlagError("Equations set has not been finished.",err,error,*999)
    IF(.NOT.ASSOCIATED(equations%equations_Set%derived)) CALL FlagError("Equations set derived is not associated.",err,error,*999)
    IF(.NOT.equations%equations_Set%derived%derivedFinished) &
      & CALL FlagError("Equations set derived has not been finished.",err,error,*999)
    IF(.NOT.ASSOCIATED(equationsSet%derived%derivedField))  &
      & CALL FlagError("Equations set derived field is not associated.",err,error,*999)
 
    equationsSet=>equations%equations_set
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
    TYPE(EQUATIONS_TYPE), POINTER, INTENT(IN) :: equations !<The equations to get the number of linear matrices for
    INTEGER(INTG), INTENT(OUT) :: numberOfMatrices !<On return, the number of linear matrices
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    !Local variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: equationsMatrices
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: linearMatrices

    ENTERS("Equations_NumberOfLinearMatricesGet",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equations%equations_Matrices)) CALL FlagError("Equations equations matrices are not associated.", &
      & err,error,*999)
    
    equationsMatrices=>equations%equations_Matrices
    linearMatrices=>equationsMatrices%linear_Matrices
    IF(ASSOCIATED(linearMatrices)) THEN
      numberOfMatrices=linearMatrices%number_of_linear_matrices
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
    TYPE(EQUATIONS_TYPE), POINTER, INTENT(IN) :: equations !<The equations to get the number of Jacobian matrices for
    INTEGER(INTG), INTENT(OUT) :: numberOfMatrices !<On return, the number of Jacobian matrices
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    !Local variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: equationsMatrices
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: nonlinearMatrices

    ENTERS("Equations_NumberOfJacobianMatricesGet",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equations%equations_Matrices)) CALL FlagError("Equations equations matrices are not associated.", &
      & err,error,*999)
    
    equationsMatrices=>equations%equations_Matrices
    nonlinearMatrices=>equationsMatrices%nonlinear_matrices
    IF(ASSOCIATED(nonlinearMatrices)) THEN
      numberOfMatrices=nonlinearMatrices%number_of_jacobians
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
    TYPE(EQUATIONS_TYPE), POINTER, INTENT(IN) :: equations !<The equations to get the number of dynamic matrices for
    INTEGER(INTG), INTENT(OUT) :: numberOfMatrices !<On return, the number of dynamic matrices
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    !Local variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: equationsMatrices
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: dynamicMatrices

    ENTERS("Equations_NumberOfDynamicMatricesGet",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equations%equations_Matrices)) CALL FlagError("Equations equations matrices are not associated.", &
      & err,error,*999)
    
    equationsMatrices=>equations%equations_Matrices
    dynamicMatrices=>equationsMatrices%dynamic_Matrices
    IF(ASSOCIATED(dynamicMatrices)) THEN
      numberOfMatrices=dynamicMatrices%number_Of_Dynamic_Matrices
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
    TYPE(EQUATIONS_TYPE), POINTER, INTENT(IN) :: equations !<The equations to get the linear matrix for
    INTEGER(INTG), INTENT(IN) :: matrixIndex !<The index of the linear matrix to get
    TYPE(DISTRIBUTED_MATRIX_TYPE), POINTER, INTENT(INOUT) :: matrix !<On return, the linear matrix requested. Must not be associated on entry.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    !Local variables
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: equationsMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Equations_LinearMatrixGet",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equations%equations_Matrices)) CALL FlagError("Equations equations matrices are not associated.", &
      & err,error,*999)
    IF(.NOT.ASSOCIATED(equations%equations_Matrices%linear_Matrices)) &
      & CALL FlagError("Equations matrices linear matrices are not associated.",err,error,*999)
    IF(matrixIndex<1.OR.matrixIndex>equations%equations_Matrices%linear_Matrices%number_Of_Linear_Matrices) THEN
      localError="The specified matrix index of "//TRIM(NumberToVString(matrixIndex,"*",err,error))// &
        & " is invalid. The matrix index must be > 0 and <= "// &
        & TRIM(NumberToVstring(equations%equations_Matrices%linear_Matrices%number_Of_Linear_Matrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(ASSOCIATED(matrix)) CALL FlagError("The matrix is already associated.",err,error,*999)

    equationsMatrix=>equations%equations_Matrices%linear_Matrices%matrices(matrixIndex)%ptr
    IF(ASSOCIATED(equationsMatrix)) THEN
      matrix=>equationsMatrix%matrix
    ELSE
      localError="The linear equations matrix for matrix index "//TRIM(NumberToVString(matrixIndex,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)
    END IF

    EXITS("Equations_LinearMatrixGet")
    RETURN
999 ERRORSEXITS("Equations_LinearMatrixGet",err,error)
    RETURN 1

  END SUBROUTINE Equations_LinearMatrixGet

  !
  !================================================================================================================================
  !

  !>Get a Jacobian matrix from equations
  SUBROUTINE Equations_JacobianMatrixGet(equations,residualIndex,variableType,matrix,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER, INTENT(IN) :: equations !<The equations to get the Jacobian matrix for
    INTEGER(INTG), INTENT(IN) :: residualIndex !<The index of the residual vector to get the Jacobian matrix for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type that the residual is differentiated with respect to for this Jacobian
    TYPE(DISTRIBUTED_MATRIX_TYPE), POINTER, INTENT(INOUT) :: matrix !<On return, the requested Jacobian matrix
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    !Local variables
    INTEGER(INTG) :: matrixIndex,variableIndex
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: equationsMapping
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: nonlinearMapping
    TYPE(EQUATIONS_JACOBIAN_TYPE), POINTER :: equationsJacobian

    ENTERS("Equations_JacobianMatrixGet",err,error,*999)

    !Check for pointer associations
    IF(ASSOCIATED(equations)) THEN
      equationsMapping=>equations%equations_mapping
      IF(ASSOCIATED(equationsMapping)) THEN
        nonlinearMapping=>equationsMapping%nonlinear_mapping
        IF(.NOT.ASSOCIATED(nonlinearMapping)) THEN
          CALL FlagError("The equations nonlinear mapping is not associated.",err,error,*999)
        END IF
      ELSE
        CALL FlagError("The equations mapping is not associated.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("The equations are not associated.",err,error,*999)
    END IF
    IF(ASSOCIATED(matrix)) THEN
      CALL FlagError("The matrix is already associated.",err,error,*999)
    END IF

    IF(residualIndex/=1) THEN
      CALL FlagError("Multiple residual vectors are not yet implemented.",err,error,*999)
    END IF

    !Find Jacobian matrix index using the nonlinear equations mapping
    matrixIndex=0
    DO variableIndex=1,nonlinearMapping%number_of_residual_variables
      IF(nonlinearMapping%residual_variables(variableIndex)%ptr%variable_type==variableType) THEN
        matrixIndex=nonlinearMapping%var_to_jacobian_map(variableIndex)%jacobian_number
      END IF
    END DO
    IF(matrixIndex==0) THEN
      CALL FlagError("Equations do not have a Jacobian matrix for residual index "// &
        & TRIM(NumberToVstring(residualIndex,"*",err,error))//" and variable type "// &
        & TRIM(NumberToVstring(variableType,"*",err,error))//".",err,error,*999)
    END IF

    !Now get Jacobian matrix using the matrix index
    equationsJacobian=>nonlinearMapping%jacobian_to_var_map(matrixIndex)%jacobian
    IF(ASSOCIATED(equationsJacobian)) THEN
      matrix=>equationsJacobian%jacobian
    ELSE
      CALL FlagError("The equations Jacobian matrix is not associated.",err,error,*999)
    END IF

    EXITS("Equations_JacobianMatrixGet")
    RETURN
999 ERRORSEXITS("Equations_JacobianMatrixGet",err,error)
    RETURN 1

  END SUBROUTINE Equations_JacobianMatrixGet

  !
  !================================================================================================================================
  !

  !>Get a dynamic equations matrix from equations using the dynamic matrix index
  SUBROUTINE Equations_DynamicMatrixGet(equations,matrixIndex,matrix,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER, INTENT(IN) :: equations !<The equations to get the dynamic matrix for
    INTEGER(INTG), INTENT(IN) :: matrixIndex !<The number of the dynamic matrix to get
    TYPE(DISTRIBUTED_MATRIX_TYPE), POINTER, INTENT(INOUT) :: matrix !<On return, the requested dynamic matrix. Must not be associated on entry.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    !Local variables
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: equationsMatrix
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Equations_DynamicMatrixGet",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equations%equations_Matrices)) CALL FlagError("Equations equations matrices are not associated.", &
      & err,error,*999)
    IF(.NOT.ASSOCIATED(equations%equations_Matrices%dynamic_Matrices)) &
      & CALL FlagError("Equations matrices dynamic matrices are not associated.",err,error,*999)
    IF(matrixIndex<1.OR.matrixIndex>equations%equations_Matrices%dynamic_Matrices%number_Of_Dynamic_Matrices) THEN
      localError="The specified matrix index of "//TRIM(NumberToVString(matrixIndex,"*",err,error))// &
        & " is invalid. The matrix index must be > 0 and <= "// &
        & TRIM(NumberToVstring(equations%equations_Matrices%dynamic_Matrices%number_Of_Dynamic_Matrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(ASSOCIATED(matrix)) CALL FlagError("The matrix is already associated.",err,error,*999)

    equationsMatrix=>equations%equations_Matrices%dynamic_Matrices%matrices(matrixIndex)%ptr
    IF(ASSOCIATED(equationsMatrix)) THEN
      matrix=>equationsMatrix%matrix
    ELSE
      localError="The dynamic equations matrix for matrix index "//TRIM(NumberToVString(matrixIndex,"*",err,error))// &
        & " is not associated."
      CALL FlagError(localError,err,error,*999)
    END IF

    EXITS("Equations_DynamicMatrixGet")
    RETURN
999 ERRORSEXITS("Equations_DynamicMatrixGet",err,error)
    RETURN 1

  END SUBROUTINE Equations_DynamicMatrixGet

  !
  !================================================================================================================================
  !

  !>Get a dynamic equations matrix from equations using the dynamic matrix type
  SUBROUTINE Equations_DynamicMatrixGetByType(equations,matrixType,matrix,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER, INTENT(IN) :: equations !<The equations to get the dynamic matrix for
    INTEGER(INTG), INTENT(IN) :: matrixType !<The type of the dynamic matrix to get. \see EQUATIONS_SET_CONSTANTS_DynamicMatrixTypes,EQUATIONS_SET_CONSTANTS
    TYPE(DISTRIBUTED_MATRIX_TYPE), POINTER, INTENT(INOUT) :: matrix !<On return, the requested dynamic matrix
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    !Local variables
    INTEGER(INTG) :: matrixIndex
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: equationsMatrix
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: equationsMatrices
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: dynamicMatrices
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: equationsMapping
    TYPE(EQUATIONS_MAPPING_DYNAMIC_TYPE), POINTER :: dynamicMapping

    ENTERS("Equations_DynamicMatrixGetByType",err,error,*999)

    !Check all pointer associations
    IF(ASSOCIATED(equations)) THEN
      equationsMatrices=>equations%equations_matrices
      IF(ASSOCIATED(equationsMatrices)) THEN
        dynamicMatrices=>equationsMatrices%dynamic_matrices
        IF(.NOT.ASSOCIATED(dynamicMatrices)) THEN
          CALL FlagError("The equations dynamic matrices are not associated.",err,error,*999)
        END IF
      ELSE
        CALL FlagError("The equations matrices are not associated.",err,error,*999)
      END IF
      equationsMapping=>equations%equations_mapping
      IF(ASSOCIATED(equationsMapping)) THEN
        dynamicMapping=>equationsMapping%DYNAMIC_MAPPING
        IF(.NOT.ASSOCIATED(dynamicMapping)) THEN
          CALL FlagError("The equations dynamic mapping is not associated.",err,error,*999)
        END IF
      ELSE
        CALL FlagError("The equations mapping is not associated.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("The equations are not associated.",err,error,*999)
    END IF
    IF(ASSOCIATED(matrix)) THEN
      CALL FlagError("The matrix is already associated.",err,error,*999)
    END IF

    !Now get the dynamic matrix
    !Find matrix index using the equations mapping
    SELECT CASE(matrixType)
    CASE(EQUATIONS_MATRIX_STIFFNESS)
      matrixIndex=dynamicMapping%stiffness_matrix_number
    CASE(EQUATIONS_MATRIX_DAMPING)
      matrixIndex=dynamicMapping%damping_matrix_number
    CASE(EQUATIONS_MATRIX_MASS)
      matrixIndex=dynamicMapping%mass_matrix_number
    CASE DEFAULT
      CALL FlagError("Invalid dynamic matrix type "//TRIM(NumberToVstring(matrixType,"*",err,error))// &
        & " specified.",err,error,*999)
    END SELECT
    IF(matrixIndex==0) THEN
      CALL FlagError("The equations dynamic matrices do not have a matrix with the specified type of "// &
        & TRIM(NumberToVstring(matrixType,"*",err,error))//".",err,error,*999)
    ELSE
      equationsMatrix=>dynamicMatrices%matrices(matrixIndex)%ptr
      IF(ASSOCIATED(equationsMatrix)) THEN
        matrix=>equationsMatrix%matrix
      ELSE
        CALL FlagError("The equations dynamic matrix for index "// &
          & TRIM(NumberToVstring(matrixIndex,"*",err,error))//" is not associated.",err,error,*999)
      END IF
    END IF

    EXITS("Equations_DynamicMatrixGetByType")
    RETURN
999 ERRORSEXITS("Equations_DynamicMatrixGetByType",err,error)
    RETURN 1

  END SUBROUTINE Equations_DynamicMatrixGetByType

  !
  !================================================================================================================================
  !

  !>Get the type of a dynamic matrix, eg. stiffness, damping or mass
  SUBROUTINE Equations_DynamicMatrixTypeGet(equations,matrixIndex,matrixType,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER, INTENT(IN) :: equations !<The equations to get the dynamic matrix for
    INTEGER(INTG), INTENT(IN) :: matrixIndex !<The number of the dynamic matrix to get
    INTEGER(INTG), INTENT(INOUT) :: matrixType !<On return, the type of the dynamic matrix. \see EQUATIONS_MATRICES_ROUTINES_DynamicMatrixTypes,EQUATIONS_MATRICES_ROUTINES
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    !Local variables
    TYPE(EQUATIONS_MAPPING_DYNAMIC_TYPE), POINTER :: dynamicMapping
    TYPE(VARYING_STRING) :: localError

    ENTERS("Equations_DynamicMatrixTypeGet",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equations%equations_Mapping)) CALL FlagError("Equations equations mapping is not associated.", &
      & err,error,*999)
    IF(.NOT.ASSOCIATED(equations%equations_Mapping%dynamic_Mapping)) &
      & CALL FlagError("Equations mapping dynamic mapping is not associated.",err,error,*999)
    IF(matrixIndex<1.OR.matrixIndex>equations%equations_Mapping%dynamic_Mapping%number_Of_Dynamic_Equations_Matrices) THEN
      localError="The specified matrix index of "//TRIM(NumberToVString(matrixIndex,"*",err,error))// &
        & " is invalid. The matrix index must be > 0 and <= "// &
        & TRIM(NumberToVstring(equations%equations_Mapping%dynamic_Mapping%number_Of_Dynamic_Equations_Matrices,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    dynamicMapping=>equations%equations_Mapping%dynamic_Mapping
    IF(matrixIndex==dynamicMapping%stiffness_matrix_number) THEN
      matrixType=EQUATIONS_MATRIX_STIFFNESS
    ELSE IF(matrixIndex==dynamicMapping%damping_matrix_number) THEN
      matrixType=EQUATIONS_MATRIX_DAMPING
    ELSE IF(matrixIndex==dynamicMapping%mass_matrix_number) THEN
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
    TYPE(EQUATIONS_TYPE), POINTER, INTENT(IN) :: equations !<The equations to get the right hand side vector for
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER, INTENT(INOUT) :: vector !<On return, the right hand side vector for the equations. Must not be associated on entry.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    !Local variables

    ENTERS("Equations_RhsVectorGet",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations are not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equations%equations_Matrices)) &
      & CALL FlagError("Equations equations matrices are not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equations%equations_Matrices%rhs_Vector)) &
      & CALL FlagError("Equations matrices RHS vector is not associated.",err,error,*999)
    IF(ASSOCIATED(vector)) CALL FlagError("Vector is already associated.",err,error,*999)
    
    vector=>equations%equations_Matrices%rhs_Vector%vector
    IF(.NOT.ASSOCIATED(vector)) CALL FlagError("The RHS vector vector is not associated.",err,error,*999)

    EXITS("Equations_RhsVectorGet")
    RETURN
999 ERRORSEXITS("Equations_RhsVectorGet",err,error)
    RETURN 1

  END SUBROUTINE Equations_RhsVectorGet

  !
  !================================================================================================================================
  !

  !>Get a residual vector for nonlinear equations
  SUBROUTINE Equations_ResidualVectorGet(equations,residualIndex,vector,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER, INTENT(IN) :: equations !<The equations to get the residual vector for
    INTEGER(INTG), INTENT(IN) :: residualIndex !<The index of the residual vector to get
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER, INTENT(INOUT) :: vector !<On return, the residual vector for the equations
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    !Local variables

    ENTERS("Equations_ResidualVectorGet",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations are not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equations%equations_Matrices)) &
      & CALL FlagError("Equations equations matrices are not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equations%equations_Matrices%nonlinear_Matrices)) &
      & CALL FlagError("Equations matrices nonlinear matrices is not associated.",err,error,*999)
    IF(residualIndex/=1) CALL FlagError("Multiple residual vectors are not yet implemented.",err,error,*999)
    IF(.NOT.ASSOCIATED(equations%equations_Matrices%nonlinear_Matrices%residual)) &
      & CALL FlagError("Nonlinear matrices residual is not associated.",err,error,*999)
    IF(ASSOCIATED(vector)) CALL FlagError("Vector is already associated.",err,error,*999)
    
    vector=>equations%equations_Matrices%nonlinear_Matrices%residual%vector
    IF(.NOT.ASSOCIATED(vector)) CALL FlagError("The residual vector vector is not associated.",err,error,*999)

    EXITS("Equations_ResidualVectorGet")
    RETURN
999 ERRORSEXITS("Equations_ResidualVectorGet",err,error)
    RETURN 1

  END SUBROUTINE Equations_ResidualVectorGet

  !
  !================================================================================================================================
  !

  !>Get the number of field variables that contribute to the residual vector
  SUBROUTINE Equations_ResidualNumberOfVariablesGet(equations,residualIndex,numberOfVariables,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER, INTENT(IN) :: equations !<The equations to get the residual vector number of variables for
    INTEGER(INTG), INTENT(IN) :: residualIndex !<The index of the residual vector to get the number of variables for
    INTEGER(INTG), INTENT(OUT) :: numberOfVariables !<On return, the number of variables that contribute to the residual vector
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    !Local variables

    ENTERS("Equations_ResidualNumberOfVariablesGet",err,error,*999)

    !Check for pointer associations
    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations are not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equations%equations_Mapping)) &
      & CALL FlagError("Equations equations mapping is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equations%equations_Mapping%nonlinear_mapping)) &
      & CALL FlagError("Equations mapping nonlinear mapping is not associated.",err,error,*999)
    IF(residualIndex/=1) CALL FlagError("Multiple residual vectors are not yet implemented.",err,error,*999)
    
    numberOfVariables=equations%equations_Mapping%nonlinear_mapping%number_of_residual_variables

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
    TYPE(EQUATIONS_TYPE), POINTER, INTENT(IN) :: equations !<The equations to get the residual vector variables for
    INTEGER(INTG), INTENT(IN) :: residualIndex !<The index of the residual vector to get the variables for
    INTEGER(INTG), INTENT(OUT) :: residualVariables(:) !<residualVariables(varIdx). On return, the field variable type for the varIdx'th residual variable
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    !Local variables
    INTEGER(INTG) :: numberOfVariables,variableIdx
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: equationsMapping
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: nonlinearMapping

    ENTERS("Equations_ResidualVariablesGet",err,error,*999)

    !Check for pointer associations
    IF(ASSOCIATED(equations)) THEN
      equationsMapping=>equations%equations_mapping
      IF(ASSOCIATED(equationsMapping)) THEN
        nonlinearMapping=>equationsMapping%nonlinear_mapping
        IF(.NOT.ASSOCIATED(nonlinearMapping)) THEN
          CALL FlagError("The equations nonlinear mapping is not associated.",err,error,*999)
        END IF
      ELSE
        CALL FlagError("The equations mapping is not associated.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("The equations are not associated.",err,error,*999)
    END IF

    IF(residualIndex==1) THEN
      numberOfVariables=nonlinearMapping%number_of_residual_variables
      IF(SIZE(residualVariables,1)>=numberOfVariables) THEN
        DO variableIdx=1,numberOfVariables
          residualVariables(variableIdx)=nonlinearMapping%residual_variables(variableIdx)%ptr%variable_type
        END DO
      ELSE
        CALL FlagError("residualVariables array must have size of at least "// &
          & TRIM(numberToVstring(numberOfVariables,"*",err,error))//".",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Multiple residual vectors are not yet implemented.",err,error,*999)
    END IF

    EXITS("Equations_ResidualVariablesGet")
    RETURN
999 ERRORSEXITS("Equations_ResidualVariablesGet",err,error)
    RETURN 1

  END SUBROUTINE Equations_ResidualVariablesGet

  !
  !================================================================================================================================
  !

  !>Get the source vector for equations
  SUBROUTINE Equations_SourceVectorGet(equations,vector,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER, INTENT(IN) :: equations !<The equations to get the source vector for
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER, INTENT(INOUT) :: vector !<On return, the source vector for the equations. Must not be associated on entry.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    !Local variables

    ENTERS("Equations_SourceVectorGet",err,error,*999)

    IF(.NOT.ASSOCIATED(equations)) CALL FlagError("Equations are not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(equations%equations_Matrices)) &
      & CALL FlagError("Equations equations matrices are not associated.",err,error,*999)THEN
    IF(.NOT.ASSOCIATED(equations%equations_Matrices%source_Vector)) &
      & CALL FlagError("Equations matrices source vector is not associated.",err,error,*999)THEN
    IF(ASSOCIATED(vector)) CALL FlagError("Vector is already associated.",err,error,*999)
    
    vector=>equations%equations_Matricess%source_Vector%vector
    IF(.NOT.ASSOCIATED(vector)) CALL FlagError("The source vector vector is not associated.",err,error,*999)

    EXITS("Equations_SourceVectorGet")
    RETURN
999 ERRORSEXITS("Equations_SourceVectorGet",err,error)
    RETURN 1

  END SUBROUTINE Equations_SourceVectorGet

  !
  !================================================================================================================================
  !

END MODULE EquationsRoutines
