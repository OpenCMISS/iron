!> \file
!> \author Chris Bradley
!> \brief This module handles all bioelectric routines.
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

!> This module handles all bioelectric class routines.
MODULE BioelectricRoutines

  USE BaseRoutines
  USE BiodomainEquationsRoutines
  USE ControlLoopAccessRoutines
  USE EquationsSetAccessRoutines
  USE ISO_VARYING_STRING
  USE Kinds
  USE MonodomainEquationsRoutines
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

  PUBLIC Bioelectric_PostLoop

  PUBLIC Bioelectric_EquationsSetSpecificationSet

  PUBLIC Bioelectric_FiniteElementCalculate

  PUBLIC Bioelectric_EquationsSetSetup

  PUBLIC Bioelectric_EquationsSetSolutionMethodSet

  PUBLIC Bioelectric_ProblemSpecificationSet

  PUBLIC Bioelectric_ProblemSetup
  
  PUBLIC Bioelectric_PreSolve,Bioelectric_PostSolve
  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Executes after each loop of a control loop for bioelectric problems, i.e., after each time step for a time loop
  SUBROUTINE Bioelectric_PostLoop(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: loopType,pSpecification(2)
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("Bioelectric_PostLoop",err,error,*999)

    CALL ControlLoop_LoopTypeGet(controlLoop,loopType,err,error,*999)
    SELECT CASE(loopType)
    CASE(CONTROL_TIME_LOOP_TYPE)
      NULLIFY(problem)
      CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
      CALL Problem_SpecificationGet(problem,2,pSpecification,err,error,*999)
      SELECT CASE(pSpecification(2))
      CASE(PROBLEM_MONODOMAIN_EQUATION_TYPE,PROBLEM_BIDOMAIN_EQUATION_TYPE)
        CALL Biodomain_PostLoop(controlLoop,err,error,*999)
      CASE(PROBLEM_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE)
        CALL Monodomain_PostLoop(controlLoop,err,error,*999)
      CASE DEFAULT
        localError="Problem type "//TRIM(NumberToVString(pSpecification(2),"*",err,error))// &
          & " is not valid for a bioelectric problem class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      !do nothing
    END SELECT
 
    EXITS("Bioelectric_PostLoop")
    RETURN
999 ERRORSEXITS("Bioelectric_PostLoop",err,error)
    RETURN 1
    
  END SUBROUTINE Bioelectric_PostLoop

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a bioelectric equation set class.
  SUBROUTINE Bioelectric_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Bioelectric_EquationsSetSpecificationSet",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(SIZE(specification,1)<2) THEN
      localError="The size of the specified specification array of "// &
        & TRIM(NumberToVString(SIZE(specification,1),"*",err,error))//" is invalid. The size should be >= 2."
      CALL FlagError(localError,err,error,*999)
    END IF
    
    SELECT CASE(specification(2))
    CASE(EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE)
      CALL Biodomain_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
    CASE(EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE)
      CALL Biodomain_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
    CASE DEFAULT
      localError="The second equations set specification of "//TRIM(NumberToVstring(specification(2),"*",err,error))// &
        & " is not valid for a bioelectric equations set."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Bioelectric_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("Bioelectric_EquationsSetSpecificationSet",err,error)
    EXITS("Bioelectric_EquationsSetSpecificationSet")
    RETURN 1
    
  END SUBROUTINE Bioelectric_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matries and rhs vector for the given element number for a bioelectric class finite element equation set.
  SUBROUTINE Bioelectric_FiniteElementCalculate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calcualate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: esSpecification(2)
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Bioelectric_FiniteElementCalculate",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,2,esSpecification,err,error,*999)
    
    SELECT CASE(equationsSet%SPECIFICATION(2))
    CASE(EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE)
      CALL Biodomain_FiniteElementCalculate(equationsSet,elementNumber,err,error,*999)
    CASE(EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE)
      CALL Biodomain_FiniteElementCalculate(equationsSet,elementNumber,err,error,*999)
    CASE DEFAULT
      localError="Equations set type "//TRIM(NumberToVString(esSpecification(2),"*",err,error))// &
        & " is not valid for a bioelectric equation set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("Bioelectric_FiniteElementCalculate")
    RETURN
999 ERRORSEXITS("Bioelectric_FiniteElementCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE Bioelectric_FiniteElementCalculate

  !
  !================================================================================================================================
  !

  !>Sets up the equations set for a bioelectric equations set class.
  SUBROUTINE Bioelectric_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    TYPE(EquationsSetSetupType), INTENT(INOUT) :: equationsSetSetup !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: esSpecification(2)
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Bioelectric_EquationsSetSetup",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,2,esSpecification,err,error,*999)
    
    SELECT CASE(esSpecification(2))
    CASE(EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE)
      CALL Biodomain_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*999)
    CASE(EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE)
      CALL Biodomain_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*999)
    CASE DEFAULT
      localError="Equation set type "//TRIM(NumberToVString(esSpecification(2),"*",err,error))// &
        & " is not valid for a bioelectric equation set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("Bioelectric_EquationsSetSetup")
    RETURN
999 ERRORSEXITS("Bioelectric_EquationsSetSetup",err,error)
    RETURN 1
    
  END SUBROUTINE Bioelectric_EquationsSetSetup
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a bioelectric equation set class.
  SUBROUTINE Bioelectric_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: solutionMethod !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: esSpecification(2)
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Bioelectric_EquationsSetSolutionMethodSet",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,2,esSpecification,err,error,*999)
    
    SELECT CASE(esSpecification(2))
    CASE(EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE)
      CALL Biodomain_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*999)
    CASE(EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE)
      CALL Biodomain_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*999)
    CASE DEFAULT
      localError="Equations set equation type of "//TRIM(NumberToVString(esSpecification(2),"*",err,error))// &
        & " is not valid for a bioelectric equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("Bioelectric_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORSEXITS("Bioelectric_EquationsSetSolutionMethodSet",err,error)
    RETURN 1
    
  END SUBROUTINE Bioelectric_EquationsSetSolutionMethodSet

  !
  !================================================================================================================================
  !

  !>Perform pre-solve actions for the bioelectrics problem class.
  SUBROUTINE Bioelectric_PreSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(2)
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError

    ENTERS("Bioelectric_PreSolve",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,2,pSpecification,err,error,*999)
    SELECT CASE(pSpecification(2))
    CASE(PROBLEM_MONODOMAIN_EQUATION_TYPE,PROBLEM_BIDOMAIN_EQUATION_TYPE)
      CALL Biodomain_PreSolve(solver,err,error,*999)
    CASE(PROBLEM_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE)
      CALL Monodomain_PreSolve(solver,err,error,*999)
    CASE DEFAULT
      localError="Problem type "//TRIM(NumberToVString(pSpecification(2),"*",err,error))// &
        & " is not valid for a bioelectrics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
      
    EXITS("Bioelectric_PreSolve")
    RETURN
999 ERRORSEXITS("Bioelectric_PreSolve",err,error)
    RETURN 1
    
  END SUBROUTINE Bioelectric_PreSolve

  !
  !================================================================================================================================
  !

  !>Performs post solve actions for a bioelectrics problem class.
  SUBROUTINE Bioelectric_PostSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(2)
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError

    ENTERS("Bioelectric_PostSolve",err,error,*999)
  
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,2,pSpecification,err,error,*999)
    
    SELECT CASE(pSpecification(2))
    CASE(PROBLEM_MONODOMAIN_EQUATION_TYPE,PROBLEM_BIDOMAIN_EQUATION_TYPE,PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE)
      CALL Biodomain_PostSolve(solver,err,error,*999)
    CASE(PROBLEM_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE)
      CALL Monodomain_PostSolve(solver,err,error,*999)
    CASE DEFAULT
      localError="Problem type "//TRIM(NumberToVString(pSpecification(2),"*",err,error))// &
        & " is not valid for a bioelectrics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("Bioelectric_PostSolve")
    RETURN
999 ERRORSEXITS("Bioelectric_PostSolve",err,error)
    RETURN 1
    
  END SUBROUTINE Bioelectric_PostSolve

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a bioelectric problem class.
  SUBROUTINE Bioelectric_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to set the specification for.
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: problemType
    TYPE(VARYING_STRING) :: localError

    ENTERS("Bioelectric_ProblemSpecificationSet",err,error,*999)

    IF(SIZE(problemSpecification,1)<2) THEN
      localError="The size of the specified problem specification array of "// &
        & TRIM(NumberToVString(SIZE(problemSpecification,1),"*",err,error))//" is invalid. The size should be >= 2."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    problemType=problemSpecification(2)
    SELECT CASE(problemType)
    CASE(PROBLEM_MONODOMAIN_EQUATION_TYPE,PROBLEM_BIDOMAIN_EQUATION_TYPE)
      CALL Biodomain_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
    CASE(PROBLEM_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE)
      CALL Monodomain_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
    CASE DEFAULT
      localError="The second problem specification of "//TRIM(NumberToVstring(problemType,"*",err,error))// &
        & " is not valid for a bioelectric problem."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Bioelectric_ProblemSpecificationSet")
    RETURN
999 ERRORS("Bioelectric_ProblemSpecificationSet",err,error)
    EXITS("Bioelectric_ProblemSpecificationSet")
    RETURN 1
    
  END SUBROUTINE Bioelectric_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the problem for a bioelectric problem class.
  SUBROUTINE Bioelectric_ProblemSetup(problem,problemSetup,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem
    TYPE(ProblemSetupType), INTENT(INOUT) :: problemSetup !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(2)
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Bioelectric_ProblemSetup",err,error,*999)

    CALL Problem_SpecificationGet(problem,2,pSpecification,err,error,*999)
    
    SELECT CASE(pSpecification(2))
    CASE(PROBLEM_MONODOMAIN_EQUATION_TYPE)
      CALL Biodomain_ProblemSetup(problem,problemSetup,err,error,*999)
    CASE(PROBLEM_BIDOMAIN_EQUATION_TYPE)
      CALL Biodomain_ProblemSetup(problem,problemSetup,err,error,*999)
    CASE(PROBLEM_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE)
      CALL Monodomain_ProblemSetup(problem,problemSetup,err,error,*999)
    CASE DEFAULT
      localError="Problem type "//TRIM(NumberToVString(pSpecification(2),"*",err,error))// &
        & " is not valid for a bioelectric problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("Bioelectric_ProblemSetup")
    RETURN
999 ERRORSEXITS("Bioelectric_ProblemSetup",err,error)
    RETURN 1
    
  END SUBROUTINE Bioelectric_ProblemSetup

  !
  !================================================================================================================================
  !

END MODULE BioelectricRoutines

