!> \file
!> \author Chris Bradley
!> \brief This module handles all elasticity routines.
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

!> This module handles all elasticity class routines.
MODULE ElasticityRoutines

  USE BaseRoutines
  USE ControlLoopRoutines
  USE ControlLoopAccessRoutines
  USE EquationsSetConstants
  USE EquationsSetAccessRoutines
  USE FINITE_ELASTICITY_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE Kinds
  USE LINEAR_ELASTICITY_ROUTINES
  USE PROBLEM_CONSTANTS
  USE SolverAccessRoutines
  USE Strings
  USE Types

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC Elasticity_EquationsSetSpecificationSet

  PUBLIC Elasticity_FiniteElementCalculate

  PUBLIC Elasticity_FiniteElementJacobianEvaluate,Elasticity_FiniteElementResidualEvaluate

  PUBLIC Elasticity_FiniteElementPreResidualEvaluate,Elasticity_FiniteElementPostResidualEvaluate

  PUBLIC Elasticity_EquationsSetSetup

  PUBLIC Elasticity_EquationsSetSolutionMethodSet

  PUBLIC Elasticity_EquationsSetDerivedVariableCalculate

  PUBLIC Elasticity_TensorInterpolateGaussPoint

  PUBLIC Elasticity_TensorInterpolateXi

  PUBLIC Elasticity_BoundaryConditionsAnalyticCalculate
  
  PUBLIC Elasticity_ProblemSpecificationSet

  PUBLIC Elasticity_ProblemSetup

  PUBLIC Elasticity_PreSolve,Elasticity_PostSolve

  PUBLIC Elasticity_ControlLoopPreLoop,Elasticity_ControlLoopPostLoop

  PUBLIC Elasticity_LoadIncrementApply

CONTAINS

  !
  !================================================================================================================================
  !
  
  !>Sets the problem specification for an elasticity equation set class.
  SUBROUTINE Elasticity_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Elasticity_EquationsSetSpecificationSet",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(SIZE(specification,1)<2) &
      & CALL FlagError("Equations set specification must have at least two entries for an elasticity class equations set.", &
      & err,error,*999)
    
    SELECT CASE(specification(2))
    CASE(EQUATIONS_SET_LINEAR_ELASTICITY_TYPE)
      CALL LinearElasticity_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
    CASE(EQUATIONS_SET_FINITE_ELASTICITY_TYPE)
      CALL FiniteElasticity_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
    CASE DEFAULT
      localError="The second equations set specification of "//TRIM(NumberToVstring(specification(2),"*",err,error))// &
        & " is not valid for an elasticity equations set."
      CALL FlagError(localError,err,error,*999)
    END SELECT
 
    EXITS("Elasticity_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("Elasticity_EquationsSetSpecificationSet",err,error)
    EXITS("Elasticity_EquationsSetSpecificationSet")
    RETURN 1
    
  END SUBROUTINE Elasticity_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matries and rhs vector for the given element number for an elasticity class finite element equation set.
  SUBROUTINE Elasticity_FiniteElementCalculate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calcualate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Elasticity_FiniteElementCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) &
      & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<2) &
      & CALL FlagError("Equations set specification must have at least two entries for an elasticity class equations set.", &
      & err,error,*999)
    
    SELECT CASE(equationsSet%specification(2))
    CASE(EQUATIONS_SET_LINEAR_ELASTICITY_TYPE)
      CALL LINEAR_ELASTICITY_FINITE_ELEMENT_CALCULATE(equationsSet,elementNumber,err,error,*999)
    CASE(EQUATIONS_SET_FINITE_ELASTICITY_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="Equations set type "//TRIM(NumberToVString(equationsSet%specification(2),"*",err,error))// &
        & " is not valid for an elasticity equation set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("Elasticity_FiniteElementCalculate")
    RETURN
999 ERRORSEXITS("Elasticity_FiniteElementCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE Elasticity_FiniteElementCalculate

  !
  !================================================================================================================================
  !

  !>Evaluates the Jacobian for the given element number for an elasticity class finite element equation set.
  SUBROUTINE Elasticity_FiniteElementJacobianEvaluate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calcualate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Elasticity_FiniteElementJacobianEvaluate",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated",err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) &
      & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<2) &
      & CALL FlagError("Equations set specification must have at least two entries for an elasticity class equations set.", &
      & err,error,*999)
     
    SELECT CASE(equationsSet%specification(2))
    CASE(EQUATIONS_SET_LINEAR_ELASTICITY_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_FINITE_ELASTICITY_TYPE)
      CALL FiniteElasticity_FiniteElementJacobianEvaluate(equationsSet,elementNumber,err,error,*999)
    CASE DEFAULT
      localError="Equations set type "//TRIM(NumberToVString(equationsSet%specification(2),"*",err,error))// &
        & " is not valid for an elasticity equation set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("Elasticity_FiniteElementJacobianEvaluate")
    RETURN
999 ERRORSEXITS("Elasticity_FiniteElementJacobianEvaluate",err,error)
    RETURN 1
    
  END SUBROUTINE Elasticity_FiniteElementJacobianEvaluate

  !
  !================================================================================================================================
  !

  !>Evaluates the residual and rhs vector for the given element number for an elasticity class finite element equation set.
  SUBROUTINE Elasticity_FiniteElementResidualEvaluate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calcualate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Elasticity_FiniteElementResidualEvaluate",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated",err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) &
      & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<2) &
      & CALL FlagError("Equations set specification must have at least two entries for an elasticity class equations set.", &
      & err,error,*999)
      
    SELECT CASE(equationsSet%specification(2))
    CASE(EQUATIONS_SET_LINEAR_ELASTICITY_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_FINITE_ELASTICITY_TYPE)
      CALL FiniteElasticity_FiniteElementResidualEvaluate(equationsSet,elementNumber,err,error,*999)
    CASE DEFAULT
      localError="Equations set type "//TRIM(NumberToVString(equationsSet%specification(2),"*",err,error))// &
        & " is not valid for an elasticity equation set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("Elasticity_FiniteElementResidualEvaluate")
    RETURN
999 ERRORSEXITS("Elasticity_FiniteElementResidualEvaluate",err,error)
    RETURN 1
    
  END SUBROUTINE Elasticity_FiniteElementResidualEvaluate

  !
  !================================================================================================================================
  !

  !>Pre-evaluates the residual for an elasticity class finite element equation set.
  SUBROUTINE Elasticity_FiniteElementPreResidualEvaluate(equationsSet,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Elasticity_FiniteElementPreResidualEvaluate",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) &
      & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<2) &
      & CALL FlagError("Equations set specification must have at least two entries for an elasticity class equations set.", &
      & err,error,*999)
    
    SELECT CASE(equationsSet%specification(2))      
    CASE(EQUATIONS_SET_LINEAR_ELASTICITY_TYPE)
      CALL FlagError("Cannot pre-evaluate the residual for a linear equations set.",err,error,*999)
    CASE(EQUATIONS_SET_FINITE_ELASTICITY_TYPE)
      CALL FiniteElasticity_FiniteElementPreResidualEvaluate(equationsSet,err,error,*999)
    CASE DEFAULT
      localError="Equations set type "//TRIM(NumberToVString(equationsSet%specification(2),"*",err,error))// &
          & " is not valid for an elasticity equation set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("Elasticity_FiniteElementPreResidualEvaluate")
    RETURN
999 ERRORSEXITS("Elasticity_FiniteElementPreResidualEvaluate",err,error)
    RETURN 1
    
  END SUBROUTINE Elasticity_FiniteElementPreResidualEvaluate

  !
  !================================================================================================================================
  !

  !>Post-evaluates the residual for an elasticity class finite element equation set.
  SUBROUTINE Elasticity_FiniteElementPostResidualEvaluate(equationsSet,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Elasticity_FiniteElementPostResidualEvaluate",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated",err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) &
      & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<2) &
      & CALL FlagError("Equations set specification must have at least two entries for an elasticity class equations set.", &
      & err,error,*999)
     
    SELECT CASE(equationsSet%specification(2))
    CASE(EQUATIONS_SET_LINEAR_ELASTICITY_TYPE)
      CALL FlagError("Cannot post-evaluate the residual for a linear equations set.",err,error,*999)
    CASE(EQUATIONS_SET_FINITE_ELASTICITY_TYPE)
      CALL FiniteElasticity_FiniteElementPostResidualEvaluate(equationsSet,err,error,*999)
    CASE DEFAULT
      localError="Equations set type "//TRIM(NumberToVString(equationsSet%specification(2),"*",err,error))// &
        & " is not valid for an elasticity equation set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("Elasticity_FiniteElementPostResidualEvaluate")
    RETURN
999 ERRORS("Elasticity_FiniteElementPostResidualEvaluate",err,error)
    EXITS("Elasticity_FiniteElementPostResidualEvaluate")
    RETURN 1
    
  END SUBROUTINE Elasticity_FiniteElementPostResidualEvaluate

  !
  !================================================================================================================================
  !

  !>Sets up the equations set for an elasticity equations set class.
  SUBROUTINE Elasticity_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: equationsSetSetup !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Elasticity_EquationsSetSetup",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) &
      & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<2) &
      & CALL FlagError("Equations set specification must have at least two entries for an elasticity class equations set.", &
      & err,error,*999)
      
    SELECT CASE(equationsSet%specification(2))
    CASE(EQUATIONS_SET_LINEAR_ELASTICITY_TYPE)
      CALL LINEAR_ELASTICITY_EQUATIONS_SET_SETUP(equationsSet,equationsSetSetup,err,error,*999)
    CASE(EQUATIONS_SET_FINITE_ELASTICITY_TYPE)
      CALL FINITE_ELASTICITY_EQUATIONS_SET_SETUP(equationsSet,equationsSetSetup,err,error,*999)
    CASE DEFAULT
      localError="Equation set type "//TRIM(NumberToVString(equationsSet%specification(2),"*",err,error))// &
        & " is not valid for an elasticity equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("Elasticity_EquationsSetSetup")
    RETURN
999 ERRORSEXITS("Elasticity_EquationsSetSetup",err,error)
    RETURN 1
    
  END SUBROUTINE Elasticity_EquationsSetSetup

  !
  !================================================================================================================================
  !
  
  !>Sets/changes the solution method for an elasticity equation set class.
  SUBROUTINE Elasticity_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: solutionMethod !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Elasticity_EquationsSetSolutionMethodSet",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) &
      & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<2) &
      & CALL FlagError("Equations set specification must have at least two entries for an elasticity class equations set.", &
      & err,error,*999)
    
    SELECT CASE(equationsSet%specification(2))
    CASE(EQUATIONS_SET_LINEAR_ELASTICITY_TYPE)
      CALL LinearElasticity_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*999)
    CASE(EQUATIONS_SET_FINITE_ELASTICITY_TYPE)
      CALL FiniteElasticity_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*999)
    CASE DEFAULT
      localError="Equations set equation type of "//TRIM(NumberToVString(equationsSet%specification(2),"*",err,error))// &
        & " is not valid for an elasticity equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("Elasticity_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORSEXITS("Elasticity_EquationsSetSolutionMethodSet",err,error)
    RETURN 1
    
  END SUBROUTINE Elasticity_EquationsSetSolutionMethodSet

  !
  !================================================================================================================================
  !

  !>Calculates a derived value for the elasticity equations set. \see OpenCMISS::Iron::cmfe_EquationsSet_DerivedCalculate
  SUBROUTINE Elasticity_EquationsSetDerivedVariableCalculate(equationsSet,derivedType,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER, INTENT(IN) :: equationsSet !<A pointer to the equations set to calculate the output for
    INTEGER(INTG), INTENT(IN) :: derivedType !<The derived field type to calculate. \see EquationsSetConstants_DerivedTypes.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("Elasticity_EquationsSetDerivedVariableCalculate",err,error,*999)

    CALL EquationsSet_AssertIsFinished(equationsSet,err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) &
      & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<2) &
      & CALL FlagError("Equations set specification must have at least two entries for an elasticity class equations set.", &
      & err,error,*999)
    
    SELECT CASE(equationsSet%specification(2))
    CASE(EQUATIONS_SET_LINEAR_ELASTICITY_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_FINITE_ELASTICITY_TYPE)
      CALL FiniteElasticityEquationsSet_DerivedVariableCalculate(equationsSet,derivedType, &
        & err,error,*999)
    CASE DEFAULT
      CALL FlagError("The second equations set specification of "// &
        & TRIM(NumberToVString(equationsSet%specification(2),"*",err,error))// &
        & " is not valid for an elasticity equations set.",err,error,*999)
    END SELECT

    EXITS("Elasticity_EquationsSetDerivedVariableCalculate")
    RETURN
999 ERRORS("Elasticity_EquationsSetDerivedVariableCalculate",err,error)
    EXITS("Elasticity_EquationsSetDerivedVariableCalculate")
    RETURN 1
    
  END SUBROUTINE Elasticity_EquationsSetDerivedVariableCalculate

  !
  !================================================================================================================================
  !

  !>Evaluate a tensor at a given element Gauss point.
  SUBROUTINE Elasticity_TensorInterpolateGaussPoint(equationsSet,tensorEvaluateType,gaussPointNumber,userElementNumber,values, &
    & err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER, INTENT(IN) :: equationsSet !<A pointer to the equations set to interpolate the tensor for.
    INTEGER(INTG), INTENT(IN) :: tensorEvaluateType !<The type of tensor to evaluate.
    INTEGER(INTG), INTENT(IN) :: gaussPointNumber !<The Gauss point number of the field to interpolate.
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number of the field to interpolate.
    REAL(DP), INTENT(OUT) :: values(:,:) !<On exit, the interpolated tensor values.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Elasticity_TensorInterpolateGaussPoint",err,error,*999)

    CALL EquationsSet_AssertIsFinished(equationsSet,err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<2) CALL FlagError("Equations set specification must have at least two entries.", &
      & err,error,*999)

    SELECT CASE(equationsSet%specification(2))
    CASE(EQUATIONS_SET_LINEAR_ELASTICITY_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_FINITE_ELASTICITY_TYPE)
      CALL FiniteElasticity_TensorInterpolateGaussPoint(equationsSet,tensorEvaluateType,gaussPointNumber,userElementNumber, &
        & values,err,error,*999)
    CASE DEFAULT
      localError="The second equations set specification of "// &
        & TRIM(NumberToVstring(equationsSet%specification(2),"*",err,error))// &
        & " is not valid for an elasticity equation set."
      CALL FlagError(localError,err,error,*999)      
    END SELECT

    EXITS("Elasticity_TensorInterpolateGaussPoint")
    RETURN
999 ERRORSEXITS("Elasticity_TensorInterpolateGaussPoint",err,error)
    RETURN 1
    
  END SUBROUTINE Elasticity_TensorInterpolateGaussPoint

  !
  !================================================================================================================================
  !

  !>Evaluate a tensor at a given element xi location.
  SUBROUTINE Elasticity_TensorInterpolateXi(equationsSet,tensorEvaluateType,userElementNumber,xi,values,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER, INTENT(IN) :: equationsSet !<A pointer to the equations set to interpolate the tensor for.
    INTEGER(INTG), INTENT(IN) :: tensorEvaluateType !<The type of tensor to evaluate.
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number of the field to interpolate.
    REAL(DP), INTENT(IN) :: xi(:) !<xi(xiIdx). The element xi to interpolate the field at.
    REAL(DP), INTENT(OUT) :: values(:,:) !<On exit, the interpolated tensor values.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("Elasticity_TensorInterpolateXi",err,error,*999)

    CALL EquationsSet_AssertIsFinished(equationsSet,err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) &
      & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<2) &
      & CALL FlagError("Equations set specification must have at least two entries for an elasticity class equations set.", &
      & err,error,*999)
 
    SELECT CASE(equationsSet%specification(2))
    CASE(EQUATIONS_SET_LINEAR_ELASTICITY_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_FINITE_ELASTICITY_TYPE)
      CALL FiniteElasticity_TensorInterpolateXi(equationsSet,tensorEvaluateType,userElementNumber,xi,values,err,error,*999)
    CASE DEFAULT
      CALL FlagError("The second equations set specification of "// &
        & TRIM(NumberToVstring(equationsSet%specification(2),"*",err,error))// &
        & " is not valid for an elasticity equation set.",err,error,*999)
    END SELECT

    EXITS("Elasticity_TensorInterpolateXi")
    RETURN
999 ERRORSEXITS("Elasticity_TensorInterpolateXi",err,error)
    RETURN 1
    
  END SUBROUTINE Elasticity_TensorInterpolateXi

  !
  !================================================================================================================================
  !

  !>Sets the analytic boundary conditions for an elasticity equation set class.
  SUBROUTINE Elasticity_BoundaryConditionsAnalyticCalculate(equationsSet,boundaryConditions,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to set the solution method for
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: boundaryConditions !<A pointer to the boundary conditionsn to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Elasticity_BoundaryConditionsAnalyticCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) &
      & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<2) &
      & CALL FlagError("Equations set specification must have at least two entries for an elasticity class equations set.", &
      & err,error,*999)
     
    SELECT CASE(equationsSet%specification(2))
    CASE(EQUATIONS_SET_LINEAR_ELASTICITY_TYPE)
      CALL LinearElasticity_BoundaryConditionsAnalyticCalculate(equationsSet,boundaryConditions,err,error,*999)
    CASE(EQUATIONS_SET_FINITE_ELASTICITY_TYPE)
      CALL FiniteElasticity_BoundaryConditionsAnalyticCalculate(equationsSet,boundaryConditions,err,error,*999)
    CASE DEFAULT
      localError="Equations set equation type of "//TRIM(NumberToVString(equationsSet%specification(2),"*",err,error))// &
        & " is not valid for an elasticity equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Elasticity_BoundaryConditionsAnalyticCalculate")
    RETURN
999 ERRORS("Elasticity_BoundaryConditionsAnalyticCalculate",err,error)
    EXITS("Elasticity_BoundaryConditionsAnalyticCalculate")
    RETURN 1
    
  END SUBROUTINE Elasticity_BoundaryConditionsAnalyticCalculate

  !
  !================================================================================================================================
  !

  !>Sets/changes the problem type and subtype for an elasticity problem class.
  SUBROUTINE Elasticity_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification array
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: problemType
    TYPE(VARYING_STRING) :: localError

    ENTERS("Elasticity_ProblemSpecificationSet",err,error,*999)

    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
    IF(SIZE(problemSpecification,1)<2) CALL FlagError("Elasticity problem specification requires a type.",err,error,*999)
    
    problemType=problemSpecification(2)
    SELECT CASE(problemType)
    CASE(PROBLEM_LINEAR_ELASTICITY_TYPE)
      CALL LinearElasticity_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
    CASE(PROBLEM_FINITE_ELASTICITY_TYPE)
      CALL FiniteElasticity_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
    CASE(PROBLEM_LINEAR_ELASTICITY_CONTACT_TYPE)
      CALL FLAG_ERROR("Not implemented yet.",err,error,*999)
    CASE(PROBLEM_FINITE_ELASTICITY_CONTACT_TYPE)
      CALL FiniteElasticity_ContactProblemSpecificationSet(problem,problemSpecification,err,error,*999)
    CASE DEFAULT
      localError="The second problem specification of "//TRIM(NumberToVstring(problemType,"*",err,error))// &
        & " is not valid for an elasticity problem."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Elasticity_ProblemSpecificationSet")
    RETURN
999 ERRORS("Elasticity_ProblemSpecificationSet",err,error)
    EXITS("Elasticity_ProblemSpecificationSet")
    RETURN 1
    
  END SUBROUTINE Elasticity_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the problem for an elasticity problem class.
  SUBROUTINE Elasticity_ProblemSetup(problem,problemSetup,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: problemSetup !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Elasticity_ProblemSetup",err,error,*999)

    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(problem%specification)) &
      & CALL FlagError("Problem specification is not allocated.",err,error,*999)
    IF(SIZE(problem%specification,1)<2) &
      & CALL FlagError("Problem specification must have at least two entries for an elasticity problem.",err,error,*999)
      
    SELECT CASE(problem%specification(2))
    CASE(PROBLEM_LINEAR_ELASTICITY_TYPE)
      CALL LINEAR_ELASTICITY_PROBLEM_SETUP(problem,problemSetup,err,error,*999)
    CASE(PROBLEM_FINITE_ELASTICITY_TYPE)
      CALL FINITE_ELASTICITY_PROBLEM_SETUP(problem,problemSetup,err,error,*999)
    CASE(PROBLEM_LINEAR_ELASTICITY_CONTACT_TYPE)
      CALL FlagError("Not implemented yet.",err,error,*999)
    CASE(PROBLEM_FINITE_ELASTICITY_CONTACT_TYPE)
      CALL FiniteElasticity_ContactProblemSetup(problem,problemSetup,err,error,*999)
    CASE DEFAULT
      localError="Problem type "//TRIM(NumberToVString(problem%specification(2),"*",err,error))// &
        & " is not valid for an elasticity problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("Elasticity_ProblemSetup")
    RETURN
999 ERRORSEXITS("Elasticity_ProblemSetup",err,error)
    RETURN 1
    
  END SUBROUTINE Elasticity_ProblemSetup

  !
  !================================================================================================================================
  !
  
  !>Performs pre-solve actions for an elasticity problem class.
  SUBROUTINE Elasticity_PreSolve(solver,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Elasticity_PreSolve",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    IF(.NOT.ALLOCATED(problem%specification)) &
      & CALL FlagError("Problem specification is not allocated.",err,error,*999)
    IF(SIZE(problem%specification,1)<2) &
      & CALL FlagError("Problem specification must have at least two entries for an elasticity problem.",err,error,*999)
    
    SELECT CASE(problem%specification(2))
    CASE(EQUATIONS_SET_LINEAR_ELASTICITY_TYPE)
      !Do Nothing
    CASE(EQUATIONS_SET_FINITE_ELASTICITY_TYPE)
      CALL FiniteElasticity_PreSolve(solver,err,error,*999)
    CASE(PROBLEM_LINEAR_ELASTICITY_CONTACT_TYPE)
      !Do Nothing
    CASE(PROBLEM_FINITE_ELASTICITY_CONTACT_TYPE)
      !Do Nothing
    CASE DEFAULT
      localError="Problem type "//TRIM(NumberToVString(problem%specification(2),"*",err,error))// &
        & " is not valid for an elasticity problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("Elasticity_PreSolve")
    RETURN
999 ERRORSEXITS("Elasticity_PreSolve",err,error)
    RETURN 1
    
  END SUBROUTINE Elasticity_PreSolve

  !
  !================================================================================================================================
  !
  
  !>Sets up the output type for an elasticity problem class.
  SUBROUTINE Elasticity_PostSolve(solver,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Elasticity_PostSolve",err,error,*999)

    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    IF(.NOT.ALLOCATED(problem%specification)) &
      & CALL FlagError("Problem specification is not allocated.",err,error,*999)
    IF(SIZE(controlLoop%PROBLEM%specification,1)<2) &
      & CALL FlagError("Problem specification must have at least two entries for an elasticity problem.",err,error,*999)
    
    SELECT CASE(problem%specification(2))
    CASE(EQUATIONS_SET_LINEAR_ELASTICITY_TYPE)
      !Do Nothing
    CASE(EQUATIONS_SET_FINITE_ELASTICITY_TYPE)
      CALL FiniteElasticity_PostSolve(solver,err,error,*999)
    CASE(PROBLEM_LINEAR_ELASTICITY_CONTACT_TYPE)
      !Do Nothing
    CASE(PROBLEM_FINITE_ELASTICITY_CONTACT_TYPE)
      !Do Nothing
    CASE DEFAULT
      localError="Problem type "//TRIM(NumberToVString(controlLoop%PROBLEM%specification(2),"*",err,error))// &
        & " is not valid for an elasticity problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("Elasticity_PostSolve")
    RETURN
999 ERRORSEXITS("Elasticity_PostSolve",err,error)
    RETURN 1
    
  END SUBROUTINE Elasticity_PostSolve

  !
  !================================================================================================================================
  !

  !>Executes before each loop of a control loop, ie before each time step for a time loop
  SUBROUTINE Elasticity_ControlLoopPreLoop(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("Elasticity_ControlLoopPreLoop",err,error,*999)

    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    
    SELECT CASE(controlLoop%loopType)
    CASE(CONTROL_TIME_LOOP_TYPE)
      IF(.NOT.ALLOCATED(problem%specification)) &
        & CALL FlagError("Problem specification is not allocated.",err,error,*999)
      IF(SIZE(controlLoop%PROBLEM%specification,1)<2) &
        & CALL FlagError("Problem specification must have at least two entries for an elasticity problem.",err,error,*999)
        
      SELECT CASE(problem%specification(2))
      CASE(EQUATIONS_SET_LINEAR_ELASTICITY_TYPE)
        !Do nothing for now
      CASE(EQUATIONS_SET_FINITE_ELASTICITY_TYPE)
        CALL FiniteElasticity_ControlTimeLoopPreLoop(controlLoop,err,error,*999)
      CASE DEFAULT
        localError="Problem type "//TRIM(NumberToVString(controlLoop%PROBLEM%specification(2),"*",err,error))// &
          & " is not valid for an elasticity problem class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      !Do nothing
    END SELECT

    EXITS("Elasticity_ControlLoopPreLoop")
    RETURN
999 ERRORSEXITS("Elasticity_ControlLoopPreLoop",err,error)
    RETURN 1
    
  END SUBROUTINE Elasticity_ControlLoopPreLoop
  
  !
  !================================================================================================================================
  !

  !>Executes after each loop of a control loop
  SUBROUTINE Elasticity_ControlLoopPostLoop(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("Elasticity_ControlLoopPostLoop",err,error,*999)

    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    IF(.NOT.ALLOCATED(problem%specification)) CALL FlagError("Problem specification is not allocated.",err,error,*999)
    IF(SIZE(problem%specification,1)<2) &
      & CALL FlagError("Problem specification must have at least two entries for an elasticity problem.",err,error,*999)
    
    SELECT CASE(problem%specification(2))      
    CASE(PROBLEM_LINEAR_ELASTICITY_TYPE,PROBLEM_LINEAR_ELASTICITY_CONTACT_TYPE)
      !Do nothing
    CASE(PROBLEM_FINITE_ELASTICITY_TYPE,PROBLEM_FINITE_ELASTICITY_CONTACT_TYPE)
      CALL FiniteElasticity_PostLoop(controlLoop,err,error,*999)
    CASE DEFAULT
      localError="Problem type "//TRIM(NumberToVString(problem%specification(2),"*",err,error))// &
        & " is not valid for a elasticity problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Elasticity_ControlLoopPostLoop")
    RETURN
999 ERRORSEXITS("Elasticity_ControlLoopPostLoop",err,error)
    RETURN 1
    
  END SUBROUTINE Elasticity_ControlLoopPostLoop

  !
  !================================================================================================================================
  !

  !> Apply load increments for equations sets
  SUBROUTINE Elasticity_LoadIncrementApply(equationsSet,iterationNumber,maximumNumberOfIterations,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    INTEGER(INTG), INTENT(IN) :: iterationNumber !<The current load increment iteration index
    INTEGER(INTG), INTENT(IN) :: maximumNumberOfIterations !<Final index for load increment loop
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("Elasticity_LoadIncrementApply",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) &
      & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<2) &
      & CALL FlagError("Equations set specification must have at least two entries for an elasticity class equations set.", &
      & err,error,*999)
    
    SELECT CASE(equationsSet%specification(2))
    CASE(EQUATIONS_SET_FINITE_ELASTICITY_TYPE)
      CALL FINITE_ELASTICITY_LOAD_INCREMENT_APPLY(equationsSet,iterationNumber,maximumNumberOfIterations,err,error,*999)
    CASE DEFAULT
      !Do nothing
    END SELECT

    EXITS("Elasticity_LoadIncrementApply")
    RETURN
999 ERRORSEXITS("Elasticity_LoadIncrementApply",err,error)
    RETURN 1

  END SUBROUTINE Elasticity_LoadIncrementApply

  !
  !================================================================================================================================
  !

END MODULE ElasticityRoutines

