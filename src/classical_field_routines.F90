!> \file
!> \author Chris Bradley
!> \brief This module handles all classical field routines.
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

!> This module handles all classical field class routines.
MODULE ClassicalFieldRoutines

  USE AdvectionEquationsRoutines
  USE AdvectionDiffusionEquationsRoutines
  USE BaseRoutines
  USE ControlLoopAccessRoutines
  USE DiffusionEquationsRoutines
  USE EquationsSetAccessRoutines
  USE HelmholtzEquationsRoutines
  USE ISO_VARYING_STRING
  USE Kinds
  USE HamiltonJacobiRoutines
  USE LaplaceEquationsRoutines
  USE PoissonEquationsRoutines
  USE ProblemAccessRoutines
  USE ReactionDiffusionEquationsRoutines
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

  PUBLIC ClassicalField_AnalyticFunctionsEvaluate
  
  PUBLIC ClassicalField_PostLoop

  PUBLIC ClassicalField_FiniteElementJacobianEvaluate,ClassicalField_FiniteElementResidualEvaluate
  
  PUBLIC ClassicalField_EquationsSetSpecificationSet

  PUBLIC ClassicalField_FiniteElementCalculate
  
  PUBLIC ClassicalField_EquationsSetSetup

  PUBLIC ClassicalField_EquationsSetSolutionMethodSet

  PUBLIC ClassicalField_ProblemSpecificationSet

  PUBLIC ClassicalField_ProblemSetup

  PUBLIC ClassicalField_PreSolve,ClassicalField_PostSolve

  PUBLIC ClassicalField_BoundaryConditionsAnalyticCalculate

CONTAINS

  !
  !================================================================================================================================
  !

  !>Evaluate the analytic solution for a classical field equations set.
  SUBROUTINE ClassicalField_AnalyticFunctionsEvaluate(equationsSet,equationsType,analyticFunctionType,position,tangents, &
    & normal,time,variableType,globalDerivative,computerNumber,analyticParameters,materialsParameters,value,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to evaluate the analytic for
    INTEGER(INTG), INTENT(IN) :: equationsType !<The type of equation to evaluate
    INTEGER(INTG), INTENT(IN) :: analyticFunctionType !<The type of analytic function to evaluate
    REAL(DP), INTENT(IN) :: position(:) !<position(dimention_idx). The geometric position to evaluate at
    REAL(DP), INTENT(IN) :: tangents(:,:) !<tangents(dimention_idx,xi_idx). The geometric tangents at the point to evaluate at.
    REAL(DP), INTENT(IN) :: normal(:) !<normal(dimension_idx). The normal vector at the point to evaluate at.
    REAL(DP), INTENT(IN) :: time !<The time to evaluate at
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to evaluate at
    INTEGER(INTG), INTENT(IN) :: globalDerivative !<The global derivative direction to evaluate at
    INTEGER(INTG), INTENT(IN) :: computerNumber !<The dependent field component number to evaluate
    REAL(DP), INTENT(IN) :: analyticParameters(:) !<A pointer to any analytic field parameters
    REAL(DP), INTENT(IN) :: materialsParameters(:) !<A pointer to any materials field parameters
    REAL(DP), INTENT(OUT) :: value !<On return, the analtyic function value.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("ClassicalField_AnalyticFunctionsEvaluate",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    
    SELECT CASE(equationsType)
    CASE(EQUATIONS_SET_LAPLACE_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_HJ_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_POISSON_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_HELMHOLTZ_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_WAVE_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_DIFFUSION_EQUATION_TYPE)
      CALL Diffusion_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,position, &
        & tangents,normal,time,variableType,globalDerivative,computerNumber,analyticParameters, &
        & materialsParameters,value,err,error,*999)
    CASE(EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_ADVECTION_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_BIHARMONIC_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="Equations set equation type "//TRIM(NumberToVString(equationsType,"*",err,error))// &
        & " is not valid for a classical field equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("ClassicalField_AnalyticFunctionsEvaluate")
    RETURN
999 ERRORSEXITS("ClassicalField_AnalyticFunctionsEvaluate",err,error)
    RETURN 1
    
  END SUBROUTINE ClassicalField_AnalyticFunctionsEvaluate

  !
  !================================================================================================================================
  !

  !>Executes after each loop of a control loop for bioelectric problems, i.e., after each time step for a time loop
  SUBROUTINE ClassicalField_PostLoop(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("ClassicalField_PostLoop",err,error,*999)

    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*999)

    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    SELECT CASE(controlLoop%loopType)
    CASE(CONTROL_TIME_LOOP_TYPE)
      SELECT CASE(problem%specification(2))
      CASE(PROBLEM_DIFFUSION_EQUATION_TYPE)
        CALL Diffusion_PostLoop(controlLoop,err,error,*999)
      CASE(PROBLEM_REACTION_DIFFUSION_EQUATION_TYPE)
        CALL ReactionDiffusion_PostLoop(controlLoop,err,error,*999)
      CASE DEFAULT
        localError="The second problem specification of "// &
          & TRIM(NumberToVString(problem%specification(2),"*",err,error))// &
          & " is not valid for a classical field problem."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      !do nothing
    END SELECT

    EXITS("ClassicalField_PostLoop")
    RETURN
999 ERRORSEXITS("ClassicalField_PostLoop",err,error)
    RETURN 1
    
  END SUBROUTINE ClassicalField_PostLoop

  !
  !================================================================================================================================
  !

  !>Sets the equations set specification for a classical field equation set.
  SUBROUTINE ClassicalField_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("ClassicalField_EquationsSetSpecificationSet",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated",err,error,*999)
    IF(SIZE(specification,1)<2) THEN
      localError="The size of the specified specification array of "// &
        & TRIM(NumberToVString(SIZE(specification,1),"*",err,error))//" is invalid. The size should be >= 2."
      CALL FlagError(localError,err,error,*999)
    ENDIF
      
    SELECT CASE(specification(2))
    CASE(EQUATIONS_SET_LAPLACE_EQUATION_TYPE)
      CALL Laplace_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
    CASE(EQUATIONS_SET_POISSON_EQUATION_TYPE)
      CALL Poisson_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
    CASE(EQUATIONS_SET_HELMHOLTZ_EQUATION_TYPE)
      CALL Helmholtz_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
    CASE(EQUATIONS_SET_WAVE_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_DIFFUSION_EQUATION_TYPE)
      CALL Diffusion_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
    CASE(EQUATIONS_SET_ADVECTION_EQUATION_TYPE)
      CALL Advection_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
    CASE(EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TYPE)
      CALL AdvectionDiffusion_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
    CASE(EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE)
      CALL ReactionDiffusion_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
    CASE(EQUATIONS_SET_BIHARMONIC_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_HJ_EQUATION_TYPE)
      CALL HamiltonJacobi_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
    CASE DEFAULT
      localError="The second equations set specification of "//TRIM(NumberToVstring(specification(2),"*",err,error))// &
        & " is not valid for a classical field equations set."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("ClassicalField_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("ClassicalField_EquationsSetSpecificationSet",err,error)
    EXITS("ClassicalField_EquationsSetSpecificationSet")
    RETURN 1
    
  END SUBROUTINE ClassicalField_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matries and rhs vector for the given element number for a clasical field class finite element equation set.
  SUBROUTINE ClassicalField_FiniteElementCalculate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calcualate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("ClassicalField_FiniteElementCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated",err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) &
      & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<2) &
      & CALL FlagError("Equations set specification must have at least two entries for a classical field class equations set.", &
      & err,error,*999)
    
    SELECT CASE(equationsSet%specification(2))
    CASE(EQUATIONS_SET_LAPLACE_EQUATION_TYPE)
      CALL Laplace_FiniteElementCalculate(equationsSet,elementNumber,err,error,*999)
    CASE(EQUATIONS_SET_HJ_EQUATION_TYPE)
      CALL HamiltonJacobi_FiniteElementCalculate(equationsSet,elementNumber,err,error,*999)
    CASE(EQUATIONS_SET_POISSON_EQUATION_TYPE)
      CALL Poisson_FiniteElementCalculate(equationsSet,elementNumber,err,error,*999)
    CASE(EQUATIONS_SET_HELMHOLTZ_EQUATION_TYPE)
      CALL Helmholtz_FiniteElementCalculate(equationsSet,elementNumber,err,error,*999)
    CASE(EQUATIONS_SET_WAVE_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_DIFFUSION_EQUATION_TYPE)
      CALL Diffusion_FiniteElementCalculate(equationsSet,elementNumber,err,error,*999)
    CASE(EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TYPE)
      CALL AdvectionDiffusion_FiniteElementCalculate(equationsSet,elementNumber,err,error,*999)
    CASE(EQUATIONS_SET_ADVECTION_EQUATION_TYPE)
      CALL Advection_FiniteElementCalculate(equationsSet,elementNumber,err,error,*999)
    CASE(EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE)
      CALL ReactionDiffusion_FiniteElementCalculate(equationsSet,elementNumber,err,error,*999)
    CASE(EQUATIONS_SET_BIHARMONIC_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="Equations set type "//TRIM(NumberToVString(equationsSet%specification(2),"*",err,error))// &
        & " is not valid for a classical field equation set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("ClassicalField_FiniteElementCalculate")
    RETURN
999 ERRORSEXITS("ClassicalField_FiniteElementCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE ClassicalField_FiniteElementCalculate

  !
  !================================================================================================================================
  !

  !>Evaluates the element Jacobian matrix for the given element number for a clasical field class finite element equation set.
  SUBROUTINE ClassicalField_FiniteElementJacobianEvaluate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to evaluate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("ClassicalField_FiniteElementJacobianEvaluate",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated",err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) &
      & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<2) &
      & CALL FlagError("Equations set specification must have at least two entries for a classical field class equations set.", &
      & err,error,*999)
    
    SELECT CASE(equationsSet%specification(2))
    CASE(EQUATIONS_SET_LAPLACE_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_HJ_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_POISSON_EQUATION_TYPE)
      CALL Poisson_FiniteElementJacobianEvaluate(equationsSet,elementNumber,err,error,*999)
    CASE(EQUATIONS_SET_HELMHOLTZ_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_WAVE_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_DIFFUSION_EQUATION_TYPE)
      CALL Diffusion_FiniteElementJacobianEvaluate(equationsSet,elementNumber,err,error,*999)
    CASE(EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_ADVECTION_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_BIHARMONIC_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="Equations set type "//TRIM(NumberToVString(equationsSet%specification(2),"*",err,error))// &
        & " is not valid for a classical field equation set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("ClassicalField_FiniteElementJacobianEvaluate")
    RETURN
999 ERRORS("ClassicalField_FiniteElementJacobianEvaluate",err,error)
    EXITS("ClassicalField_FiniteElementJacobianEvaluate")
    RETURN 1
    
  END SUBROUTINE ClassicalField_FiniteElementJacobianEvaluate

  !
  !================================================================================================================================
  !

  !>Evaluates the element residual and rhs vectors for the given element number for a clasical field class finite element equation set.
  SUBROUTINE ClassicalField_FiniteElementResidualEvaluate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("ClassicalField_FiniteElementResidualEvaluate",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated",err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) &
      & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<2) &
      & CALL FlagError("Equations set specification must have at least two entries for a classical field class equations set.", &
      & err,error,*999)
    
    SELECT CASE(equationsSet%specification(2))
    CASE(EQUATIONS_SET_LAPLACE_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_HJ_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_POISSON_EQUATION_TYPE)
      CALL Poisson_FiniteElementResidualEvaluate(equationsSet,elementNumber,err,error,*999)
    CASE(EQUATIONS_SET_HELMHOLTZ_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_WAVE_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_DIFFUSION_EQUATION_TYPE)
      CALL Diffusion_FiniteElementResidualEvaluate(equationsSet,elementNumber,err,error,*999)
    CASE(EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_ADVECTION_EQUATION_TYPE)
      CALL Advection_FiniteElementCalculate(equationsSet,elementNumber,err,error,*999)
    CASE(EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_BIHARMONIC_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="Equations set type "//TRIM(NumberToVString(equationsSet%specification(2),"*",err,error))// &
        & " is not valid for a classical field equation set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("ClassicalField_FiniteElementResidualEvaluate")
    RETURN
999 ERRORS("ClassicalField_FiniteElementResidualEvaluate",err,error)
    EXITS("ClassicalField_FiniteElementResidualEvaluate")
    RETURN 1
    
  END SUBROUTINE ClassicalField_FiniteElementResidualEvaluate

  !
  !================================================================================================================================
  !

  !>Sets up the equations set for a classical field equations set class.
  SUBROUTINE ClassicalField_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    TYPE(EquationsSetSetupType), INTENT(INOUT) :: equationsSetSetup !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("ClassicalField_EquationsSetSetup",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated",err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) &
      & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<2) &
      & CALL FlagError("Equations set specification must have at least two entries for a classical field class equations set.", &
      & err,error,*999)
    
    SELECT CASE(equationsSet%specification(2))
    CASE(EQUATIONS_SET_LAPLACE_EQUATION_TYPE)
      CALL Laplace_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*999)
    CASE(EQUATIONS_SET_HJ_EQUATION_TYPE)
      CALL HamiltonJacobi_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*999)
    CASE(EQUATIONS_SET_POISSON_EQUATION_TYPE)
      CALL Poisson_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*999)
    CASE(EQUATIONS_SET_HELMHOLTZ_EQUATION_TYPE)
      CALL Helmholtz_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*999)
    CASE(EQUATIONS_SET_WAVE_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_DIFFUSION_EQUATION_TYPE)
      CALL Diffusion_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*999)
    CASE(EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TYPE)
      CALL AdvectionDiffusion_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*999)
    CASE(EQUATIONS_SET_ADVECTION_EQUATION_TYPE)
      CALL Advection_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*999)
    CASE(EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE)
      CALL ReactionDiffusion_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*999)
    CASE(EQUATIONS_SET_BIHARMONIC_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="Equation set type "//TRIM(NumberToVString(equationsSet%specification(2),"*",err,error))// &
        & " is not valid for a classical field equation set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
        
    EXITS("ClassicalField_EquationsSetSetup")
    RETURN
999 ERRORSEXITS("ClassicalField_EquationsSetSetup",err,error)
    RETURN 1
    
  END SUBROUTINE ClassicalField_EquationsSetSetup
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a classical field equation set class.
  SUBROUTINE ClassicalField_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: solutionMethod !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("ClassicalField_EquationsSetSolutionMethodSet",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated",err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) &
      & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<2) &
      & CALL FlagError("Equations set specification must have at least two entries for a classical field class equations set.", &
      & err,error,*999)
    
     SELECT CASE(equationsSet%specification(2))
     CASE(EQUATIONS_SET_LAPLACE_EQUATION_TYPE)
       CALL Laplace_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*999)
     CASE(EQUATIONS_SET_HJ_EQUATION_TYPE)
       CALL HamiltonJacobi_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*999)
     CASE(EQUATIONS_SET_POISSON_EQUATION_TYPE)
       CALL Poisson_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*999)
     CASE(EQUATIONS_SET_HELMHOLTZ_EQUATION_TYPE)
       CALL Helmholtz_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*999)
     CASE(EQUATIONS_SET_WAVE_EQUATION_TYPE)
       CALL FlagError("Not implemented.",err,error,*999)
     CASE(EQUATIONS_SET_DIFFUSION_EQUATION_TYPE)
       CALL Diffusion_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*999)
     CASE(EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TYPE)
       CALL AdvectionDiffusion_EquationsSetSolnMethodSet(equationsSet,solutionMethod,err,error,*999)
     CASE(EQUATIONS_SET_ADVECTION_EQUATION_TYPE)
       CALL Advection_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*999)
     CASE(EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE)
       CALL ReactionDiffusion_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*999)
     CASE(EQUATIONS_SET_BIHARMONIC_EQUATION_TYPE)
       CALL FlagError("Not implemented.",err,error,*999)
     CASE DEFAULT
       localError="Equations set equation type of "//TRIM(NumberToVString(equationsSet%specification(2),"*",err,error))// &
         & " is not valid for a classical field equations set class."
       CALL FlagError(localError,err,error,*999)
     END SELECT
       
    EXITS("ClassicalField_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORS("ClassicalField_EquationsSetSolutionMethodSet",err,error)
    EXITS("ClassicalField_EquationsSetSolutionMethodSet")
    RETURN 1
    
  END SUBROUTINE ClassicalField_EquationsSetSolutionMethodSet

  !
  !================================================================================================================================
  !

  !>Sets the analytic boundary conditions for a classical field equation set class.
  SUBROUTINE ClassicalField_BoundaryConditionsAnalyticCalculate(equationsSet,boundaryConditions,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the solution method for
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditionsn to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: esSpecification(2)
    TYPE(VARYING_STRING) :: localError

    ENTERS("ClassicalField_BoundaryConditionsAnalyticCalculate",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,2,esSpecification,err,error,*999)
    
    SELECT CASE(esSpecification(2))
    CASE(EQUATIONS_SET_LAPLACE_EQUATION_TYPE)
      CALL Laplace_BoundaryConditionsAnalyticCalculate(equationsSet,boundaryConditions,err,error,*999)
    CASE(EQUATIONS_SET_HJ_EQUATION_TYPE)
      CALL HamiltonJacobi_BoundaryConditionsAnalyticCalculate(equationsSet,boundaryConditions,err,error,*999)
    CASE(EQUATIONS_SET_POISSON_EQUATION_TYPE)
      CALL Poisson_BoundaryConditionsAnalyticCalculate(equationsSet,boundaryConditions,.FALSE.,err,error,*999)
    CASE(EQUATIONS_SET_HELMHOLTZ_EQUATION_TYPE)
      CALL Helmholtz_BoundaryConditionsAnalyticCalculate(equationsSet,boundaryConditions,err,error,*999)
    CASE(EQUATIONS_SET_WAVE_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_DIFFUSION_EQUATION_TYPE)
      CALL Diffusion_BoundaryConditionAnalyticCalculate(equationsSet,boundaryConditions,err,error,*999)
    CASE(EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TYPE)
      CALL AdvectionDiffusion_BoundaryConditionsAnalyticCalculate(equationsSet,boundaryConditions,err,error,*999)
    CASE(EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_BIHARMONIC_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="Equations set equation type of "//TRIM(NumberToVString(esSpecification(2),"*",err,error))// &
        & " is not valid for a classical field equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("ClassicalField_BoundaryConditionsAnalyticCalculate")
    RETURN
999 ERRORS("ClassicalField_BoundaryConditionsAnalyticCalculate",err,error)
    EXITS("ClassicalField_BoundaryConditionsAnalyticCalculate")
    RETURN 1
    
  END SUBROUTINE ClassicalField_BoundaryConditionsAnalyticCalculate

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a classical field problem class.
  SUBROUTINE ClassicalField_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to set the specification for.
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: problemType
    TYPE(VARYING_STRING) :: localError

    ENTERS("ClassicalField_ProblemSpecificationSet",err,error,*999)

    IF(.NOT.ASSOCIATED(problem))  CALL FlagError("Problem is not associated.",err,error,*999)
    IF(SIZE(problemSpecification,1)<2) THEN
      localError="The size of the specified specification array of "// &
        & TRIM(NumberToVString(SIZE(problemSpecification,1),"*",err,error))//" is invalid. The size should be >= 2."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    problemType=problemSpecification(2)
    SELECT CASE(problemType)
    CASE(PROBLEM_LAPLACE_EQUATION_TYPE)
      CALL Laplace_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
    CASE(PROBLEM_HJ_EQUATION_TYPE)
      CALL HamiltonJacobi_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
    CASE(PROBLEM_POISSON_EQUATION_TYPE)
      CALL Poisson_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
    CASE(PROBLEM_HELMHOLTZ_EQUATION_TYPE)
      CALL Helmholtz_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
    CASE(PROBLEM_WAVE_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(PROBLEM_DIFFUSION_EQUATION_TYPE)
      CALL Diffusion_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
    CASE(PROBLEM_ADVECTION_EQUATION_TYPE)
      CALL Advection_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
    CASE(PROBLEM_ADVECTION_DIFFUSION_EQUATION_TYPE)
      CALL AdvectionDiffusion_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
    CASE(PROBLEM_REACTION_DIFFUSION_EQUATION_TYPE)
      CALL ReactionDiffusion_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
    CASE(PROBLEM_BIHARMONIC_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The second problem specification of "//TRIM(NumberToVstring(problemType,"*",err,error))// &
        & " is not valid for a classical field problem."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("ClassicalField_ProblemSpecificationSet")
    RETURN
999 ERRORS("ClassicalField_ProblemSpecificationSet",err,error)
    EXITS("ClassicalField_ProblemSpecificationSet")
    RETURN 1
    
  END SUBROUTINE ClassicalField_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the problem for a classical field problem class.
  SUBROUTINE ClassicalField_ProblemSetup(problem,problemSetup,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem
    TYPE(ProblemSetupType), INTENT(INOUT) :: problemSetup !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("ClassicalField_ProblemSetup",err,error,*999)

    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(problem%SPECIFICATION)) CALL FlagError("Problem specification is not allocated.",err,error,*999)
    IF(SIZE(problem%specification,1)<2) &
      & CALL FlagError("Problem specification must have at least two entries for a classical field problem.",err,error,*999)
    
    SELECT CASE(problem%specification(2))
    CASE(PROBLEM_LAPLACE_EQUATION_TYPE)
      CALL Laplace_ProblemSetup(problem,problemSetup,err,error,*999)
    CASE(PROBLEM_POISSON_EQUATION_TYPE)
      CALL Poisson_ProblemSetup(problem,problemSetup,err,error,*999)
    CASE(PROBLEM_HELMHOLTZ_EQUATION_TYPE)
      CALL Helmholtz_ProblemSetup(problem,problemSetup,err,error,*999)
    CASE(PROBLEM_WAVE_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(PROBLEM_DIFFUSION_EQUATION_TYPE)
      CALL Diffusion_ProblemSetup(problem,problemSetup,err,error,*999)
    CASE(PROBLEM_ADVECTION_EQUATION_TYPE)
      CALL Advection_ProblemSetup(problem,problemSetup,err,error,*999)
    CASE(PROBLEM_ADVECTION_DIFFUSION_EQUATION_TYPE)
      CALL AdvectionDiffusion_ProblemSetup(problem,problemSetup,err,error,*999)
    CASE(PROBLEM_REACTION_DIFFUSION_EQUATION_TYPE)
      CALL ReactionDiffusion_ProblemSetup(problem,problemSetup,err,error,*999)
    CASE(PROBLEM_BIHARMONIC_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(PROBLEM_HJ_EQUATION_TYPE)
      CALL HamiltonJacobi_ProblemSetup(problem,problemSetup,err,error,*999)
    CASE DEFAULT
      localError="Problem type "//TRIM(NumberToVString(problem%specification(2),"*",err,error))// &
        & " is not valid for a classical field problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
        
    EXITS("ClassicalField_ProblemSetup")
    RETURN
999 ERRORSEXITS("ClassicalField_ProblemSetup",err,error)
    RETURN 1
    
  END SUBROUTINE ClassicalField_ProblemSetup

  !
  !================================================================================================================================
  !

  !>Performs any pre-solve operations for a classical field problem class.
  SUBROUTINE ClassicalField_PreSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("ClassicalField_PreSolve",err,error,*999)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
    
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    IF(.NOT.ALLOCATED(problem%specification)) CALL FlagError("Problem specification is not allocated.",err,error,*999)
    IF(SIZE(problem%specification,1)<2) &
      & CALL FlagError("Problem specification must have at least two entries for a classical field problem.",err,error,*999)
    
    SELECT CASE(problem%specification(2))      
    CASE(PROBLEM_LAPLACE_EQUATION_TYPE)
      !CALL Laplace_PreSolve(solver,err,error,*999)
    CASE(PROBLEM_HJ_EQUATION_TYPE)
      !CALL HJ_PreSolve(solver,err,error,*999)
    CASE(PROBLEM_POISSON_EQUATION_TYPE)
      CALL Poisson_PreSolve(solver,err,error,*999)
    CASE(PROBLEM_HELMHOLTZ_EQUATION_TYPE)
      !CALL Helmholtz_PreSolve(solver,err,error,*999)
    CASE(PROBLEM_WAVE_EQUATION_TYPE)
      !CALL Wave_PreSolve(solver,err,error,*999)
    CASE(PROBLEM_DIFFUSION_EQUATION_TYPE)
      CALL Diffusion_PreSolve(solver,err,error,*999)
    CASE(PROBLEM_ADVECTION_EQUATION_TYPE)
      CALL Advection_PreSolve(solver,err,error,*999)
    CASE(PROBLEM_ADVECTION_DIFFUSION_EQUATION_TYPE)
      CALL AdvectionDiffusion_PreSolve(solver,err,error,*999)
    CASE(PROBLEM_REACTION_DIFFUSION_EQUATION_TYPE)
      CALL ReactionDiffusion_PreSolve(solver,err,error,*999)
    CASE(PROBLEM_BIHARMONIC_EQUATION_TYPE)
      !CALL Biharmonic_PreSolve(solver,err,error,*999)
    CASE DEFAULT
      localError="Problem type "//TRIM(NumberToVString(problem%specification(2),"*",err,error))// &
        & " is not valid for a classical field problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("ClassicalField_PreSolve")
    RETURN
999 ERRORSEXITS("ClassicalField_PreSolve",err,error)
    RETURN 1
    
  END SUBROUTINE ClassicalField_PreSolve

  !
  !================================================================================================================================
  !

  !>Performs post solve operations a classical field problem class.
  SUBROUTINE ClassicalField_PostSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("ClassicalField_PostSolve",err,error,*999)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
    
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    IF(.NOT.ALLOCATED(problem%specification)) CALL FlagError("Problem specification is not allocated.",err,error,*999)
    IF(SIZE(problem%specification,1)<2) &
      & CALL FlagError("Problem specification must have at least two entries for a classical field problem.",err,error,*999)
    
    SELECT CASE(problem%specification(2))
    CASE(PROBLEM_LAPLACE_EQUATION_TYPE)
      !CALL Laplace_PostSolve(solver,err,error,*999)
    CASE(PROBLEM_HJ_EQUATION_TYPE)
      !CALL HJ_PostSolve(solver,err,error,*999)
    CASE(PROBLEM_POISSON_EQUATION_TYPE)
      CALL Poisson_PostSolve(solver,err,error,*999)
    CASE(PROBLEM_HELMHOLTZ_EQUATION_TYPE)
      !CALL Helmholtz_PostSolve,err,error,*999)
    CASE(PROBLEM_WAVE_EQUATION_TYPE)
      !CALL Wave_PostSolve(solver,err,error,*999)
    CASE(PROBLEM_DIFFUSION_EQUATION_TYPE)
      CALL Diffusion_PostSolve(solver,err,error,*999)
    CASE(PROBLEM_ADVECTION_EQUATION_TYPE)
      !CALL Advection_PostSolve(solver,err,error,*999)
    CASE(PROBLEM_ADVECTION_DIFFUSION_EQUATION_TYPE)
      CALL AdvectionDiffusion_PostSolve(solver,err,error,*999)
    CASE(PROBLEM_REACTION_DIFFUSION_EQUATION_TYPE)
      CALL ReactionDiffusion_PostSolve(solver,err,error,*999)
    CASE(PROBLEM_BIHARMONIC_EQUATION_TYPE)
      !CALL Biharmonic_PostSolve(solver,err,error,*999)
    CASE DEFAULT
      localError="Problem type "//TRIM(NumberToVString(problem%specification(2),"*",err,error))// &
        & " is not valid for a classical field problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("ClassicalField_PostSolve")
    RETURN
999 ERRORSEXITS("ClassicalField_PostSolve",err,error)
    RETURN 1
    
  END SUBROUTINE ClassicalField_PostSolve

  !
  !================================================================================================================================
  !

END MODULE ClassicalFieldRoutines

