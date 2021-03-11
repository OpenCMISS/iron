!> \file
!> $Id: multi_physics_routines.F90 177 2009-04-20 
!> \authors Christian Michler, Jack Lee
!> \brief This module handles all multi physics routines.
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
!> Contributor(s): Chris Bradley, Christian Michler, Jack Lee
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

!> This module handles all multi physics class routines.
MODULE MultiPhysicsRoutines

  USE BaseRoutines
  USE DiffusionAdvectionDiffusionRoutines
  USE DiffusionDiffusionRoutines
  USE FiniteElasticityDarcyRoutines
  USE FiniteElasticityFluidPressureRoutines
  USE FSIRoutines
  USE BioelectricFiniteElasticityRoutines
  USE EquationsSetAccessRoutines
  USE ISO_VARYING_STRING
  USE Kinds
  USE MultiCompartmentTransportRoutines
  USE NavierStokesEquationsRoutines
  USE ProblemAccessRoutines
  USE Strings
  USE Types

#include "macros.h"  


  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC MultiPhysics_FiniteElementJacobianEvaluate

  PUBLIC MultiPhysics_FiniteElementResidualEvaluate
  
  PUBLIC MultiPhysics_EquationsSetSpecificationSet

  PUBLIC MultiPhysics_FiniteElementCalculate

  PUBLIC MultiPhysics_EquationsSetSetup

  PUBLIC MultiPhysics_EquationsSetSolnMethodSet

  PUBLIC MultiPhysics_ProblemSpecificationSet

  PUBLIC MultiPhysics_ProblemSetup

  PUBLIC MultiPhysics_PostSolve

  PUBLIC MultiPhysics_PreSolve

  PUBLIC MultiPhysics_PreLoop

  PUBLIC MultiPhysics_PostLoop
  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a multi physics equation set class.
  SUBROUTINE MultiPhysics_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("MultiPhysics_EquationsSetSpecificationSet",err,error,*999)

    !Note that in general, this routine is never used as most multi-physics problems
    !use standard equations sets and couples them, rather than having a special
    !multi-physics problem equations set

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated",err,error,*999)
    IF(SIZE(specification,1)<2) THEN
      localError="The size of the specified specification array of "// &
        & TRIM(NumberToVString(SIZE(specification,1),"*",err,error))//" is invalid. The size should be >= 2."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(specification(2))
    CASE(EQUATIONS_SET_FINITE_ELASTICITY_DARCY_TYPE)
      CALL FinElasticityFluidPressure_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
    CASE(EQUATIONS_SET_FINITE_ELASTICITY_STOKES_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_FINITE_ELASTICITY_NAVIER_STOKES_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_DIFFUSION_DIFFUSION_TYPE)
      CALL DiffusionDiffusion_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
    CASE(EQUATIONS_SET_DIFFUSION_ADVECTION_DIFFUSION_TYPE)
      CALL DiffusionAdvectionDiffusion_EquationsSetSpecSet(equationsSet,specification,err,error,*999)
    CASE DEFAULT
      localError="The second equations set specification of "//TRIM(NumberToVstring(specification(2),"*",err,error))// &
        & " is not valid for a multi physics equations set."
      CALL FlagError(localError,err,error,*999)
    END SELECT
 
    EXITS("MultiPhysics_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("MultiPhysics_EquationsSetSpecificationSet",err,error)
    EXITS("MultiPhysics_EquationsSetSpecificationSet")
    RETURN 1
    
  END SUBROUTINE MultiPhysics_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matries and rhs vector for the given element number for a multi physics class finite element equation set.
  SUBROUTINE MultiPhysics_FiniteElementCalculate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calcualate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: esSpecification(2)
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("MultiPhysics_FiniteElementCalculate",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,2,esSpecification,err,error,*999)

    SELECT CASE(esSpecification(2))
    CASE(EQUATIONS_SET_FINITE_ELASTICITY_DARCY_TYPE)
      CALL FinElasticityFluidPressure_FiniteElementCalculate(equationsSet,elementNumber,err,error,*999)
    CASE(EQUATIONS_SET_FINITE_ELASTICITY_STOKES_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_FINITE_ELASTICITY_NAVIER_STOKES_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_DIFFUSION_DIFFUSION_TYPE)
      CALL DiffusionDiffusion_FiniteElementCalculate(equationsSet,elementNumber,err,error,*999)
    CASE(EQUATIONS_SET_DIFFUSION_ADVECTION_DIFFUSION_TYPE)
      CALL DiffusionAdvectionDiffusion_FiniteElementCalculate(equationsSet,elementNumber,err,error,*999)
    CASE DEFAULT
      localError="Equations set type "//TRIM(NumberToVString(esSpecification(2),"*",err,error))// &
        & " is not valid for a multi physics equation set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("MultiPhysics_FiniteElementCalculate")
    RETURN
999 ERRORSEXITS("MultiPhysics_FiniteElementCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE MultiPhysics_FiniteElementCalculate

  !
  !================================================================================================================================
  !

  !>Evaluates the element Jacobian matrix for the given element number for a multi physics class finite element equation set.
  SUBROUTINE MultiPhysics_FiniteElementJacobianEvaluate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to evaluate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: esSpecification(2)
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("MultiPhysics_FiniteElementJacobianEvaluate",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,2,esSpecification,err,error,*999)

    SELECT CASE(esSpecification(2))
    CASE(EQUATIONS_SET_FINITE_ELASTICITY_DARCY_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_FINITE_ELASTICITY_STOKES_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_FINITE_ELASTICITY_NAVIER_STOKES_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_DIFFUSION_DIFFUSION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_DIFFUSION_ADVECTION_DIFFUSION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="Equations set type "//TRIM(NumberToVString(esSpecification(2),"*",err,error))// &
        & " is not valid for a multi physics equation set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("MultiPhysics_FiniteElementJacobianEvaluate")
    RETURN
999 ERRORSEXITS("MultiPhysics_FiniteElementJacobianEvaluate",err,error)
    RETURN 1
    
  END SUBROUTINE MultiPhysics_FiniteElementJacobianEvaluate

  !
  !================================================================================================================================
  !

  !>Evaluates the element residual and rhs vectors for the given element number for a multi physics class finite element equation set.
  SUBROUTINE MultiPhysics_FiniteElementResidualEvaluate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: esSpecification(2)
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("MultiPhysics_FiniteElementResidualEvaluate",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,2,esSpecification,err,error,*999)

    SELECT CASE(esSpecification(2))
    CASE(EQUATIONS_SET_FINITE_ELASTICITY_DARCY_TYPE)
      !CALL ELASTICITY_DARCY_FINITE_ELEMENT_RESIDUAL_EVALUATE(equationsSet,elementNumber,err,error,*999)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_FINITE_ELASTICITY_STOKES_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_FINITE_ELASTICITY_NAVIER_STOKES_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_DIFFUSION_DIFFUSION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_DIFFUSION_ADVECTION_DIFFUSION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="Equations set type "//TRIM(NumberToVString(esSpecification(2),"*",err,error))// &
        & " is not valid for a multi physics equation set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("MultiPhysics_FiniteElementResidualEvaluate")
    RETURN
999 ERRORSEXITS("MultiPhysics_FiniteElementResidualEvaluate",err,error)
    RETURN 1
    
  END SUBROUTINE MultiPhysics_FiniteElementResidualEvaluate

  !
  !================================================================================================================================
  !

  !>Sets up the equations set for a multi physics equations set class.
  SUBROUTINE MultiPhysics_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    TYPE(EquationsSetSetupType), INTENT(INOUT) :: equationsSetSetup !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: esSpecification(2)
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("MultiPhysics_EquationsSetSetup",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,2,esSpecification,err,error,*999)

    SELECT CASE(esSpecification(2))
    CASE(EQUATIONS_SET_FINITE_ELASTICITY_DARCY_TYPE)
      CALL FiniteElasticityDarcy_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*999)
    CASE(EQUATIONS_SET_FINITE_ELASTICITY_STOKES_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_FINITE_ELASTICITY_NAVIER_STOKES_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_DIFFUSION_DIFFUSION_TYPE)
      CALL DiffusionDiffusion_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*999)
    CASE(EQUATIONS_SET_DIFFUSION_ADVECTION_DIFFUSION_TYPE)
      CALL DiffusionAdvectionDiffusion_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*999)
    CASE DEFAULT
      localError="Equation set type "//TRIM(NumberToVString(esSpecification(2),"*",err,error))// &
        & " is not valid for a multi physics equation set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
      
    EXITS("MultiPhysics_EquationsSetSetup")
    RETURN
999 ERRORSEXITS("MultiPhysics_EquationsSetSetup",err,error)
    RETURN 1
    
  END SUBROUTINE MultiPhysics_EquationsSetSetup
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a multi physics equation set class.
  SUBROUTINE MultiPhysics_EquationsSetSolnMethodSet(equationsSet,solutionMethod,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: solutionMethod !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: esSpecification(2)
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("MultiPhysics_EquationsSetSolnMethodSet",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,2,esSpecification,err,error,*999)
    
    SELECT CASE(esSpecification(2))
    CASE(EQUATIONS_SET_FINITE_ELASTICITY_DARCY_TYPE)
      CALL FinElasticityFluidPressure_EquationsSetSolnMethodSet(equationsSet,solutionMethod,err,error,*999)
    CASE(EQUATIONS_SET_FINITE_ELASTICITY_STOKES_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_FINITE_ELASTICITY_NAVIER_STOKES_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_DIFFUSION_DIFFUSION_TYPE)
      CALL DiffusionDiffusion_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*999)
    CASE(EQUATIONS_SET_DIFFUSION_ADVECTION_DIFFUSION_TYPE)
      CALL DiffusionAdvectionDiffusion_EquationsSetSolnMethodSet(equationsSet,solutionMethod,err,error,*999)
    CASE DEFAULT
      localError="Equations set equation type of "//TRIM(NumberToVString(esSpecification(2),"*",err,error))// &
        & " is not valid for a multi physics equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("MultiPhysics_EquationsSetSolnMethodSet")
    RETURN
999 ERRORS("MultiPhysics_EquationsSetSolnMethodSet",err,error)
    EXITS("MultiPhysics_EquationsSetSolnMethodSet")
    RETURN 1
    
  END SUBROUTINE MultiPhysics_EquationsSetSolnMethodSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the problem type and subtype for a multi physics problem class.
  SUBROUTINE MultiPhysics_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemType

    ENTERS("MultiPhysics_ProblemSpecificationSet",err,error,*999)

    IF(.NOT.ASSOCIATED(problem))  CALL FlagError("Problem is not associated.",err,error,*999)
    IF(SIZE(problemSpecification,1)<2) THEN
      localError="The size of the specified specification array of "// &
        & TRIM(NumberToVString(SIZE(problemSpecification,1),"*",err,error))//" is invalid. The size should be >= 2."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    problemType=problemSpecification(2)
    SELECT CASE(problemType)
    CASE(PROBLEM_FINITE_ELASTICITY_DARCY_TYPE)
      CALL FiniteElasticityDarcy_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
    CASE(PROBLEM_FINITE_ELASTICITY_FLUID_PRESSURE_TYPE)
      CALL FinElasticityFluidPressure_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
    CASE(PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE)
      CALL BioelectricFiniteElasticity_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
    CASE(PROBLEM_FINITE_ELASTICITY_STOKES_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_TYPE)
      CALL FSI_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
    CASE(PROBLEM_DIFFUSION_DIFFUSION_TYPE)
      CALL DiffusionDiffusion_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
    CASE(PROBLEM_DIFFUSION_ADVECTION_DIFFUSION_TYPE)
      CALL DiffusionAdvectionDiffusion_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
    CASE(PROBLEM_MULTI_COMPARTMENT_TRANSPORT_TYPE)
      CALL MultiCompartmentTransport_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
    CASE DEFAULT
      localError="The second problem specification of "//TRIM(NumberToVstring(problemType,"*",err,error))// &
        & " is not valid for a multi physics problem."
      CALL FlagError(localError,err,error,*999)
    END SELECT
 
    EXITS("MultiPhysics_ProblemSpecificationSet")
    RETURN
999 ERRORS("MultiPhysics_ProblemSpecificationSet",err,error)
    EXITS("MultiPhysics_ProblemSpecificationSet")
    RETURN 1
    
  END SUBROUTINE MultiPhysics_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the problem for a multi physics problem class.
  SUBROUTINE MultiPhysics_ProblemSetup(problem,problemSetup,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem
    TYPE(ProblemSetupType), INTENT(INOUT) :: problemSetup !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(2)
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("MultiPhysics_ProblemSetup",err,error,*999)

    CALL Problem_SpecificationGet(problem,2,pSpecification,err,error,*999)
    
    SELECT CASE(pSpecification(2))
    CASE(PROBLEM_FINITE_ELASTICITY_DARCY_TYPE)
      CALL FiniteElasticityDarcy_ProblemSetup(problem,problemSetup,err,error,*999)
    CASE(PROBLEM_FINITE_ELASTICITY_FLUID_PRESSURE_TYPE)
      CALL FinElasticityFluidPressure_ProblemSetup(problem,problemSetup,err,error,*999)
    CASE(PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE)
      CALL BioelectricFiniteElasticity_ProblemSetup(problem,problemSetup,err,error,*999) 
    CASE(PROBLEM_FINITE_ELASTICITY_STOKES_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_TYPE)
      CALL FSI_ProblemSetup(problem,problemSetup,err,error,*999)
    CASE(PROBLEM_DIFFUSION_DIFFUSION_TYPE)
      CALL DiffusionDiffusion_ProblemSetup(problem,problemSetup,err,error,*999)
    CASE(PROBLEM_DIFFUSION_ADVECTION_DIFFUSION_TYPE)
      CALL DiffusionAdvectionDiffusion_ProblemSetup(problem,problemSetup,err,error,*999)
    CASE(PROBLEM_MULTI_COMPARTMENT_TRANSPORT_TYPE)
      CALL MultiCompartmentTransport_ProblemSetup(problem,problemSetup,err,error,*999)
    CASE DEFAULT
      localError="Problem type "//TRIM(NumberToVString(pSpecification(2),"*",err,error))// &
        & " is not valid for a multi physics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("MultiPhysics_ProblemSetup")
    RETURN
999 ERRORSEXITS("MultiPhysics_ProblemSetup",err,error)
    RETURN 1
    
  END SUBROUTINE MultiPhysics_ProblemSetup

  !
  !================================================================================================================================
  !

  !>Sets up the output type for a multi physics problem class.
  SUBROUTINE MultiPhysics_PostSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(2)
    TYPE(ControlLoopType), POINTER :: controlLoop 
    TYPE(ProblemType), POINTER :: problem 
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("MultiPhysics_PostSolve",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,2,pSpecification,err,error,*999)

    SELECT CASE(pSpecification(2))
    CASE(PROBLEM_FINITE_ELASTICITY_DARCY_TYPE)
      CALL FiniteElasticityDarcy_PostSolve(solver,err,error,*999)
    CASE(PROBLEM_FINITE_ELASTICITY_FLUID_PRESSURE_TYPE)
      CALL FinElasticityFluidPressure_PostSolve(solver,err,error,*999)
    CASE(PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE)
      CALL BioelectricFiniteElasticity_PostSolve(solver,err,error,*999)
    CASE(PROBLEM_FINITE_ELASTICITY_STOKES_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_TYPE)
      CALL FSI_PostSolve(solver,err,error,*999)
    CASE(PROBLEM_DIFFUSION_DIFFUSION_TYPE)
      CALL DiffusionDiffusion_PostSolve(solver,err,error,*999)
    CASE(PROBLEM_DIFFUSION_ADVECTION_DIFFUSION_TYPE)
      !CALL DiffusionAdvectionDiffusion_PostSolve(solver,err,error,*999)
    CASE(PROBLEM_MULTI_COMPARTMENT_TRANSPORT_TYPE)
      CALL MultiCompartmentTransport_PostSolve(solver,err,error,*999)
    CASE DEFAULT
      localError="Problem type "//TRIM(NumberToVString(pSpecification(2),"*",err,error))// &
        & " is not valid for a multi physics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
      
    EXITS("MultiPhysics_PostSolve")
    RETURN
999 ERRORSEXITS("MultiPhysics_PostSolve",err,error)
    RETURN 1
    
  END SUBROUTINE MultiPhysics_PostSolve

  !
  !================================================================================================================================
  !

  !>Sets up the output type for a multi physics problem class.
  SUBROUTINE MultiPhysics_PreSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(2)
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("MultiPhysics_PreSolve",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,2,pSpecification,err,error,*999)
    
    SELECT CASE(pSpecification(2))
    CASE(PROBLEM_FINITE_ELASTICITY_DARCY_TYPE)
      CALL FiniteElasticityDarcy_PreSolve(solver,err,error,*999)
    CASE(PROBLEM_FINITE_ELASTICITY_FLUID_PRESSURE_TYPE)
      CALL FinElasticityFluidPressure_PreSolve(solver,err,error,*999)
    CASE(PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE)
      CALL BioelectricFiniteElasticity_PreSolve(solver,err,error,*999)
    CASE(PROBLEM_FINITE_ELASTICITY_STOKES_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_TYPE)
      CALL FSI_PreSolve(solver,err,error,*999)
    CASE(PROBLEM_DIFFUSION_DIFFUSION_TYPE)
      CALL DiffusionDiffusion_PreSolve(solver,err,error,*999)
    CASE(PROBLEM_DIFFUSION_ADVECTION_DIFFUSION_TYPE)
      CALL DiffusionAdvectionDiffusion_PreSolve(solver,err,error,*999)
    CASE(PROBLEM_MULTI_COMPARTMENT_TRANSPORT_TYPE)
      CALL MultiCompartmentTransport_PreSolve(solver,err,error,*999)
    CASE DEFAULT
      localError="Problem type "//TRIM(NumberToVString(pSpecification(2),"*",err,error))// &
        & " is not valid for a multi physics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
      
    EXITS("MultiPhysics_PreSolve")
    RETURN
999 ERRORSEXITS("MultiPhysics_PreSolve",err,error)
    RETURN 1
    
  END SUBROUTINE MultiPhysics_PreSolve

  !
  !================================================================================================================================
  !

  !>Executes before each loop of a control loop, ie before each time step for a time loop
  SUBROUTINE MultiPhysics_PreLoop(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(2)
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("MultiPhysics_PreLoop",err,error,*999)

    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,2,pSpecification,err,error,*999)

    SELECT CASE(pSpecification(2))
    CASE(PROBLEM_FINITE_ELASTICITY_DARCY_TYPE)
      CALL FiniteElasticityDarcy_PreLoop(controlLoop,err,error,*999)
    CASE(PROBLEM_FINITE_ELASTICITY_FLUID_PRESSURE_TYPE)
      !do nothing
    CASE(PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE)
      CALL BioelectricFiniteElasticity_PreLoop(controlLoop,err,error,*999)
    CASE(PROBLEM_FINITE_ELASTICITY_STOKES_TYPE)
      !do nothing
    CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_TYPE)
      !TODO Store previous data?
    CASE(PROBLEM_DIFFUSION_DIFFUSION_TYPE)
      !do nothing
    CASE(PROBLEM_DIFFUSION_ADVECTION_DIFFUSION_TYPE)
      !do nothing
    CASE(PROBLEM_MULTI_COMPARTMENT_TRANSPORT_TYPE)
      !do nothing
    CASE DEFAULT
      localError="Problem type "//TRIM(NumberToVString(pSpecification(2),"*",err,error))// &
        & " is not valid for a multi physics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("MultiPhysics_PreLoop")
    RETURN
999 ERRORSEXITS("MultiPhysics_PreLoop",err,error)
    RETURN 1
    
  END SUBROUTINE MultiPhysics_PreLoop

  !
  !================================================================================================================================
  !

  !>Executes after each loop of a control loop, ie after each time step for a time loop
  SUBROUTINE MultiPhysics_PostLoop(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(2)
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("MultiPhysics_PostLoop",err,error,*999)

    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,2,pSpecification,err,error,*999)

    SELECT CASE(pSpecification(2))
    CASE(PROBLEM_FINITE_ELASTICITY_DARCY_TYPE)
      CALL FiniteElasticityDarcy_PostLoop(controlLoop,err,error,*999)
    CASE(PROBLEM_FINITE_ELASTICITY_FLUID_PRESSURE_TYPE)
      !do nothing
    CASE(PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE)
      CALL BioelectricFiniteElasticity_PostLoop(controlLoop,err,error,*999)
    CASE(PROBLEM_FINITE_ELASTICITY_STOKES_TYPE)
      !do nothing
    CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_TYPE)
      CALL FSI_PostLoop(controlLoop,err,error,*999)
    CASE(PROBLEM_DIFFUSION_DIFFUSION_TYPE)
      !do nothing
    CASE(PROBLEM_DIFFUSION_ADVECTION_DIFFUSION_TYPE)
      !do nothing
    CASE(PROBLEM_MULTI_COMPARTMENT_TRANSPORT_TYPE)
      !do nothing
    CASE DEFAULT
      localError="Problem type "//TRIM(NumberToVString(pSpecification(2),"*",err,error))// &
        & " is not valid for a multi physics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("MultiPhysics_PostLoop")
    RETURN
999 ERRORSEXITS("MultiPhysics_PostLoop",err,error)
    RETURN 1
    
  END SUBROUTINE MultiPhysics_PostLoop

  !
  !================================================================================================================================
  !

END MODULE MultiPhysicsRoutines

