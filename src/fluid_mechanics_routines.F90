!> \file
!> \author Sebastian Krittian
!> \brief This module handles all fluid mechanics routines.
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
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> Kingdom. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
!>
!> Contributor(s): Sebastian Krittian, David Ladd, Chris Bradley
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

!> This module handles all fluid mechanics class routines.
MODULE FluidMechanicsRoutines

  USE AdvectionEquationsRoutines
  USE BaseRoutines
  USE BurgersEquationsRoutines
  USE CharacteristicEquationsRoutines
  USE ControlLoopRoutines
  USE ControlLoopAccessRoutines
  USE DarcyEquationsRoutines
  USE DarcyPressureEquationsRoutines
  USE EquationsSetAccessRoutines
  USE InputOutput
  USE ISO_VARYING_STRING
  USE Kinds
  USE NavierStokesEquationsRoutines
  USE PoiseuilleEquationsRoutines
  USE ProblemAccessRoutines
  USE StokesEquationsRoutines
  USE SolverAccessRoutines
  USE Strings
  USE StreeEquationsRoutines
  USE Types

#include "macros.h"  


  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC FluidMechanics_AnalyticFunctionsEvaluate

  PUBLIC FluidMechanics_BoundaryConditionsAnalyticCalculate

  PUBLIC FluidMechanics_PreLoop,FluidMechanics_PostLoop
  
  PUBLIC FluidMechanics_EquationsSetSolutionMethodSet

  PUBLIC FluidMechanics_EquationsSetSetup

  PUBLIC FluidMechanics_EquationsSetSpecificationSet

  PUBLIC FluidMechanics_FiniteElementCalculate

  PUBLIC FluidMechanics_FiniteElementJacobianEvaluate,FluidMechanics_FiniteElementResidualEvaluate
  
  PUBLIC FluidMechanics_FiniteElementPreResidualEvaluate

  PUBLIC FluidMechanics_NodalJacobianEvaluate,FluidMechanics_NodalResidualEvaluate
  
  PUBLIC FluidMechanics_ProblemSetup
   
  PUBLIC FluidMechanics_PostSolve,FluidMechanics_PreSolve

  PUBLIC FluidMechanics_ProblemSpecificationSet
 
CONTAINS

  !
  !================================================================================================================================
  !

  !>Evaluate the analytic solution for a fluid mechanics equations set.
  SUBROUTINE FluidMechanics_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,position,tangents, &
    & normal,time,variableType,globalDerivative,componentNumber,analyticParameters,materialsParameters,value,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to evaluate the analytic for
    INTEGER(INTG), INTENT(IN) :: analyticFunctionType !<The type of analytic function to evaluate
    REAL(DP), INTENT(IN) :: position(:) !<position(dimentionIdx). The geometric position to evaluate at
    REAL(DP), INTENT(IN) :: tangents(:,:) !<tangents(dimentionIdx,xiIdx). The geometric tangents at the point to evaluate at.
    REAL(DP), INTENT(IN) :: normal(:) !<normal(dimensionIdx). The normal vector at the point to evaluate at.
    REAL(DP), INTENT(IN) :: time !<The time to evaluate at
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to evaluate at
    INTEGER(INTG), INTENT(IN) :: globalDerivative !<The global derivative direction to evaluate at
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The dependent field component number to evaluate
    REAL(DP), INTENT(IN) :: analyticParameters(:) !<A pointer to any analytic field parameters
    REAL(DP), INTENT(IN) :: materialsParameters(:) !<A pointer to any materials field parameters
    REAL(DP), INTENT(OUT) :: value !<On return, the analtyic function value.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("FluidMechanics_AnalyticFunctionsEvaluate",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) &
      & CALL FlagError("Equations set specification has not been allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<2) &
      & CALL FlagError("Equations set specification must have at least two entries for a fluid mechanics equations set.", &
      & err,error,*999)
    
    SELECT CASE(equationsSet%specification(2))
    CASE(EQUATIONS_SET_STOKES_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_DARCY_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_DARCY_PRESSURE_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_POISEUILLE_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_BURGERS_EQUATION_TYPE)
      CALL Burgers_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,position, &
        & tangents,normal,time,variableType,globalDerivative,componentNumber,analyticParameters, &
        & materialsParameters,value,err,error,*999)
    CASE(EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The second equations set specification of "// &
        & TRIM(NumberToVString(equationsSet%specification(2),"*",err,error))// &
        & " is not valid for a fluid mechanics equations set."
      CALL FlagError(localError,err,error,*999)
    END SELECT
        
    EXITS("FluidMechanics_AnalyticFunctionsEvaluate")
    RETURN
999 ERRORSEXITS("FluidMechanics_AnalyticFunctionsEvaluate",err,error)
    RETURN 1
    
  END SUBROUTINE FluidMechanics_AnalyticFunctionsEvaluate

   !
  !================================================================================================================================
  !

  !>Sets the problem specification for a fluid mechanics equation set class.
  SUBROUTINE FluidMechanics_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("FluidMechanics_EquationsSetSpecificationSet",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(SIZE(specification,1)<2) &
       & CALL FlagError("Equations set specification must have at least two entries for a fluid mechanics equations set.", &
       & err,error,*999)
      
    SELECT CASE(specification(2))
    CASE(EQUATIONS_SET_STOKES_EQUATION_TYPE)
      CALL Stokes_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
    CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE)
      CALL NavierStokes_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
    CASE(EQUATIONS_SET_DARCY_EQUATION_TYPE)
      CALL Darcy_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
    CASE(EQUATIONS_SET_DARCY_PRESSURE_EQUATION_TYPE)
      CALL DarcyPressure_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
    CASE(EQUATIONS_SET_POISEUILLE_EQUATION_TYPE)
      CALL Poiseuille_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
    CASE(EQUATIONS_SET_BURGERS_EQUATION_TYPE)
      CALL Burgers_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
    CASE(EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE)
      CALL Characteristic_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
    CASE(EQUATIONS_SET_STREE_EQUATION_TYPE)
      CALL Stree_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
    CASE DEFAULT
      localError="The second equations set specification of "//TRIM(NumberToVstring(specification(2),"*",err,error))// &
        & " is not valid for a fluid mechanics equations set."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("FluidMechanics_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("FluidMechanics_EquationsSetSpecificationSet",err,error)
    EXITS("FluidMechanics_EquationsSetSpecificationSet")
    RETURN 1
    
  END SUBROUTINE FluidMechanics_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matries and rhs vector for the given element number for a fluid mechanics class finite element equation set.
  SUBROUTINE FluidMechanics_FiniteElementCalculate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calcualate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("FluidMechanics_FiniteElementCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated",err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) &
       & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<2) &
      & CALL FlagError("Equations set specification must have at least two entries for a fluid mechanics class equations set.", &
      & err,error,*999)
      
    SELECT CASE(equationsSet%specification(2))
    CASE(EQUATIONS_SET_STOKES_EQUATION_TYPE)
      CALL Stokes_FiniteElementCalculate(equationsSet,elementNumber,err,error,*999)
    CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE)
      CALL FlagError("There are no finite element matrices to be calculated for Navier-Stokes equation.",err,error,*999)
    CASE(EQUATIONS_SET_DARCY_EQUATION_TYPE)
      CALL Darcy_FiniteElementCalculate(equationsSet,elementNumber,err,error,*999)
    CASE(EQUATIONS_SET_DARCY_PRESSURE_EQUATION_TYPE)
      CALL FlagError("There is no element stiffness matrix to be calculated for Darcy pressure.",err,error,*999)
    CASE(EQUATIONS_SET_POISEUILLE_EQUATION_TYPE)
      CALL Poiseuille_FiniteElementCalculate(equationsSet,elementNumber,err,error,*999)
    CASE(EQUATIONS_SET_BURGERS_EQUATION_TYPE)
      CALL FlagError("There are no finite element matrices to be calculated for Burgers equation.",err,error,*999)
    CASE(EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE)
      CALL FlagError("There are no finite element matrices to be calculated for Characteristic equations.",err,error,*999)
    CASE(EQUATIONS_SET_STREE_EQUATION_TYPE)
      CALL Stree_FiniteElementCalculate(equationsSet,elementNumber,err,error,*999)
    CASE DEFAULT
      localError="Equations set type "//TRIM(NumberToVString(equationsSet%specification(2),"*",err,error))// &
        & " is not valid for a fluid mechanics equation set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("FluidMechanics_FiniteElementCalculate")
    RETURN
999 ERRORSEXITS("FluidMechanics_FiniteElementCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE FluidMechanics_FiniteElementCalculate

  !
  !================================================================================================================================
  !

  !>Evaluates the element Jacobian matrix for the given element number for a fluid mechanics class finite element equation set.
  SUBROUTINE FluidMechanics_FiniteElementJacobianEvaluate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to evaluate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("FluidMechanics_FiniteElementJacobianEvaluate",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated",err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) &
      & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<2) &
      & CALL FlagError("Equations set specification must have at least two entries for a fluid mechanics class equations set.", &
      & err,error,*999)
      
    SELECT CASE(equationsSet%specification(2))
    CASE(EQUATIONS_SET_STOKES_EQUATION_TYPE)
      CALL FlagError("There is no Jacobian to be evaluated for Stokes.",err,error,*999)
    CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE)
      CALL NavierStokes_FiniteElementJacobianEvaluate(equationsSet,elementNumber,err,error,*999)
    CASE(EQUATIONS_SET_DARCY_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_DARCY_PRESSURE_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_POISEUILLE_EQUATION_TYPE)
      CALL FlagError("There is no Jacobian to be evaluated for Poiseuille.",err,error,*999)
    CASE(EQUATIONS_SET_BURGERS_EQUATION_TYPE)
      CALL Burgers_FiniteElementJacobianEvaluate(equationsSet,elementNumber,err,error,*999)
    CASE(EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="Equations set type "//TRIM(NumberToVString(equationsSet%specification(2),"*",err,error))// &
        & " is not valid for a fluid mechanics equation set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("FluidMechanics_FiniteElementJacobianEvaluate")
    RETURN
999 ERRORS("FluidMechanics_FiniteElementJacobianEvaluate",err,error)
    EXITS("FluidMechanics_FiniteElementJacobianEvaluate")
    RETURN 1
    
  END SUBROUTINE FluidMechanics_FiniteElementJacobianEvaluate

  !
  !================================================================================================================================
  !

  !>Evaluates the element residual and rhs vectors for the given element number for a fluid mechanics class finite element equation set.
  SUBROUTINE FluidMechanics_FiniteElementResidualEvaluate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("FluidMechanics_FiniteElementResidualEvaluate",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated",err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) &
      & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<2) &
      & CALL FlagError("Equations set specification must have at least two entries for a fluid mechanics class equations set.", &
      & err,error,*999)
     
    SELECT CASE(equationsSet%specification(2))
    CASE(EQUATIONS_SET_STOKES_EQUATION_TYPE)
      CALL FlagError("There is no residual to be evaluated for Stokes.",err,error,*999)
    CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE)
      CALL NavierStokes_FiniteElementResidualEvaluate(equationsSet,elementNumber,err,error,*999)
    CASE(EQUATIONS_SET_DARCY_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_DARCY_PRESSURE_EQUATION_TYPE)
      CALL DarcyPressure_FiniteElementResidualEvaluate(equationsSet,elementNumber,err,error,*999)
    CASE(EQUATIONS_SET_POISEUILLE_EQUATION_TYPE)
      CALL FlagError("There is no residual to be evaluated for Poiseuille.",err,error,*999)
    CASE(EQUATIONS_SET_BURGERS_EQUATION_TYPE)
      CALL Burgers_FiniteElementResidualEvaluate(equationsSet,elementNumber,err,error,*999)
    CASE(EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="Equations set type "//TRIM(NumberToVString(equationsSet%specification(2),"*",err,error))// &
        & " is not valid for a fluid mechanics equation set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("FluidMechanics_FiniteElementResidualEvaluate")
    RETURN
999 ERRORS("FluidMechanics_FiniteElementResidualEvaluate",err,error)
    EXITS("FluidMechanics_FiniteElementResidualEvaluate")
    RETURN 1
    
  END SUBROUTINE FluidMechanics_FiniteElementResidualEvaluate

  !
  !================================================================================================================================
  !

  !>Evaluates the nodal Jacobian matrix for the given node number for a fluid mechanics class nodal equation set.
  SUBROUTINE FluidMechanics_NodalJacobianEvaluate(equationsSet,nodeNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: nodeNumber !<The element number to evaluate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("FluidMechanics_NodalJacobianEvaluate",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated",err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) &
      & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<2) &
      & CALL FlagError("Equations set specification must have at least two entries for a fluid mechanics class equations set.", &
      & err,error,*999)
      
    SELECT CASE(equationsSet%specification(2))
    CASE(EQUATIONS_SET_STOKES_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_DARCY_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_DARCY_PRESSURE_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_POISEUILLE_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_BURGERS_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE)
      CALL Characteristic_NodalJacobianEvaluate(equationsSet,nodeNumber,err,error,*999)
    CASE DEFAULT
      localError="Equations set type "//TRIM(NumberToVString(equationsSet%specification(2),"*",err,error))// &
        & " is not valid for a fluid mechanics equation set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("FluidMechanics_NodalJacobianEvaluate")
    RETURN
999 ERRORSEXITS("FluidMechanics_NodalJacobianEvaluate",err,error)
    RETURN 1
    
  END SUBROUTINE FluidMechanics_NodalJacobianEvaluate

  !
  !================================================================================================================================
  !

  !>Evaluates the nodal residual and rhs vectors for the given node number for a fluid mechanics class nodal equation set.
  SUBROUTINE FluidMechanics_NodalResidualEvaluate(equationsSet,nodeNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: nodeNumber !<The element number to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("FluidMechanics_NodalResidualEvaluate",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated",err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) &
      & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<2) &
      & CALL FlagError("Equations set specification must have at least two entries for a fluid mechanics class equations set.", &
      & err,error,*999)
      
    SELECT CASE(equationsSet%specification(2))
    CASE(EQUATIONS_SET_STOKES_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_DARCY_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_DARCY_PRESSURE_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_POISEUILLE_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_BURGERS_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE)
      CALL Characteristic_NodalResidualEvaluate(equationsSet,nodeNumber,err,error,*999)
    CASE DEFAULT
      localError="Equations set type "//TRIM(NumberToVString(equationsSet%specification(2),"*",err,error))// &
        & " is not valid for a fluid mechanics equation set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("FluidMechanics_NodalResidualEvaluate")
    RETURN
999 ERRORSEXITS("FluidMechanics_NodalResidualEvaluate",err,error)
    RETURN 1
    
  END SUBROUTINE FluidMechanics_NodalResidualEvaluate

  !
  !================================================================================================================================
  !

  !>Sets up the equations set for a fluid mechanics equations set class.
  SUBROUTINE FluidMechanics_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    TYPE(EquationsSetSetupType), INTENT(INOUT) :: equationsSetSetup !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("FluidMechanics_EquationsSetSetup",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) &
      & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<2) &
      & CALL FlagError("Equations set specification must have at least two entries for a fluid mechanics class equations set.", &
      & err,error,*999)
      
    SELECT CASE(equationsSet%specification(2))
    CASE(EQUATIONS_SET_STOKES_EQUATION_TYPE)
      CALL Stokes_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*999)
    CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE)
      CALL NavierStokes_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*999)
    CASE(EQUATIONS_SET_DARCY_EQUATION_TYPE)
      CALL Darcy_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*999)
    CASE(EQUATIONS_SET_DARCY_PRESSURE_EQUATION_TYPE)
      CALL DarcyPressure_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*999)
    CASE(EQUATIONS_SET_POISEUILLE_EQUATION_TYPE)
      CALL Poiseuille_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*999)
    CASE(EQUATIONS_SET_BURGERS_EQUATION_TYPE)
      CALL Burgers_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*999)
    CASE(EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE)
      CALL Characteristic_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*999)
    CASE(EQUATIONS_SET_STREE_EQUATION_TYPE)
      CALL Stree_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*999)
    CASE DEFAULT
      localError="Equation set type "//TRIM(NumberToVString(equationsSet%specification(2),"*",err,error))// &
        & " is not valid for a fluid mechanics equation set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("FluidMechanics_EquationsSetSetup")
    RETURN
999 ERRORSEXITS("FluidMechanics_EquationsSetSetup",err,error)
    RETURN 1
    
  END SUBROUTINE FluidMechanics_EquationsSetSetup
  

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a fluid mechanics equation set class.
  SUBROUTINE FluidMechanics_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: solutionMethod !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("FluidMechanics_EquationsSetSolutionMethodSet",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) &
      & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<2) &
      & CALL FlagError("Equations set specification must have at least two entries for a fluid mechanics class equations set.", &
      & err,error,*999)
      
    SELECT CASE(equationsSet%specification(2))
    CASE(EQUATIONS_SET_STOKES_EQUATION_TYPE)
      CALL Stokes_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*999)
    CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE)
      CALL NavierStokes_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*999)
    CASE(EQUATIONS_SET_DARCY_EQUATION_TYPE)
      CALL DarcyPressure_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*999)
    CASE(EQUATIONS_SET_DARCY_PRESSURE_EQUATION_TYPE)
      CALL DarcyPressure_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*999)
    CASE(EQUATIONS_SET_POISEUILLE_EQUATION_TYPE)
      CALL Poiseuille_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*999)
    CASE(EQUATIONS_SET_BURGERS_EQUATION_TYPE)
      CALL Burgers_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*999)
    CASE(EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE)
      CALL Characteristic_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*999)
    CASE DEFAULT
      localError="Equations set equation type of "//TRIM(NumberToVString(equationsSet%specification(2),"*",err,error))// &
        & " is not valid for a fluid mechancis equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("FluidMechanics_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORS("FluidMechanics_EquationsSetSolutionMethodSet",err,error)
    EXITS("FluidMechanics_EquationsSetSolutionMethodSet")
    RETURN 1
    
  END SUBROUTINE FluidMechanics_EquationsSetSolutionMethodSet

  !
  !================================================================================================================================
  !

  !>Sets the analytic boundary conditions for a fluid mechanics equation set class.
  SUBROUTINE FluidMechanics_BoundaryConditionsAnalyticCalculate(equationsSet,boundaryConditions,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the solution method for
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditionsn to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("FluidMechanics_BoundaryConditionsAnalyticCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated",err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) &
      & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<2) &
      & CALL FlagError("Equations set specification must have at least two entries for a fluid mechanics class equations set.", &
      & err,error,*999)
      
    SELECT CASE(equationsSet%specification(2))
    CASE(EQUATIONS_SET_BURGERS_EQUATION_TYPE)
      CALL Burgers_BoundaryConditionsAnalyticCalculate(equationsSet,boundaryConditions,err,error,*999)
    CASE(EQUATIONS_SET_STOKES_EQUATION_TYPE)
      CALL Stokes_BoundaryConditionsAnalyticCalculate(equationsSet,boundaryConditions,err,error,*999)
    CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE)
      CALL NavierStokes_BoundaryConditionsAnalyticCalculate(equationsSet,boundaryConditions,err,error,*999)
    CASE(EQUATIONS_SET_DARCY_EQUATION_TYPE)
      CALL Darcy_BoundaryConditionsAnalyticCalculate(equationsSet,boundaryConditions,err,error,*999)
    CASE(EQUATIONS_SET_DARCY_PRESSURE_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_POISEUILLE_EQUATION_TYPE)
      CALL Poiseuille_BoundaryConditionsAnalyticCalculate(equationsSet,boundaryConditions,err,error,*999)
    CASE(EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="Equations set equation type of "//TRIM(NumberToVString(equationsSet%specification(2),"*",err,error))// &
        & " is not valid for a fluid mechanics equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("FluidMechanics_BoundaryConditionsAnalyticCalculate")
    RETURN
999 ERRORS("FluidMechanics_BoundaryConditionsAnalyticCalculate",err,error)
    EXITS("FluidMechanics_BoundaryConditionsAnalyticCalculate")
    RETURN 1
    
  END SUBROUTINE FluidMechanics_BoundaryConditionsAnalyticCalculate

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a fluid mechanics problem class.
  SUBROUTINE FluidMechanics_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to set the specification for.
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemType

    ENTERS("FluidMechanics_ProblemSpecificationSet",err,error,*999)

    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
    IF(SIZE(problemSpecification,1)<2) &
      & CALL FlagError("Fluid mechanics problem specification must have a type set.",err,error,*999)

    problemType=problemSpecification(2)
    SELECT CASE(problemType)
    CASE(PROBLEM_STOKES_EQUATION_TYPE)
      CALL Stokes_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
    CASE(PROBLEM_NAVIER_STOKES_EQUATION_TYPE)
      CALL NavierStokes_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
    CASE(PROBLEM_DARCY_EQUATION_TYPE)
      CALL Darcy_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
    CASE(PROBLEM_POISEUILLE_EQUATION_TYPE)
      CALL Poiseuille_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
    CASE(PROBLEM_BURGERS_EQUATION_TYPE)
      CALL Burgers_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
    CASE DEFAULT
      localError="The second problem specification of "//TRIM(NumberToVstring(problemType,"*",err,error))// &
        & " is not valid for a fluid mechanics problem."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("FluidMechanics_ProblemSpecificationSet")
    RETURN
999 ERRORS("FluidMechanics_ProblemSpecificationSet",err,error)
    EXITS("FluidMechanics_ProblemSpecificationSet")
    RETURN 1
    
  END SUBROUTINE FluidMechanics_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the problem for a fluid mechanics problem class.
  SUBROUTINE FluidMechanics_ProblemSetup(problem,problemSetup,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem
    TYPE(ProblemSetupType), INTENT(INOUT) :: problemSetup !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("FluidMechanics_ProblemSetup",err,error,*999)

    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(problem%specification)) &
      & CALL FlagError("Problem specification is not allocated.",err,error,*999)
    IF(SIZE(problem%specification,1)<2) &
      & CALL FlagError("Problem specification must have at least two entries for a fluid mechanics problem.",err,error,*999)
      
    SELECT CASE(problem%specification(2))
    CASE(PROBLEM_STOKES_EQUATION_TYPE)
      CALL Stokes_ProblemSetup(problem,problemSetup,err,error,*999)
    CASE(PROBLEM_NAVIER_STOKES_EQUATION_TYPE)
      CALL NavierStokes_ProblemSetup(problem,problemSetup,err,error,*999)
    CASE(PROBLEM_DARCY_EQUATION_TYPE)
      CALL Darcy_ProblemSetup(problem,problemSetup,err,error,*999)
    CASE(PROBLEM_POISEUILLE_EQUATION_TYPE)
      CALL Poiseuille_ProblemSetup(problem,problemSetup,err,error,*999)
    CASE(PROBLEM_BURGERS_EQUATION_TYPE)
      CALL Burgers_ProblemSetup(problem,problemSetup,err,error,*999)
    CASE DEFAULT
      localError="Problem type "//TRIM(NumberToVString(problem%specification(2),"*",err,error))// &
        & " is not valid for a fluid mechanics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("FluidMechanics_ProblemSetup")
    RETURN
999 ERRORSEXITS("FluidMechanics_ProblemSetup",err,error)
    RETURN 1
    
  END SUBROUTINE FluidMechanics_ProblemSetup

  !
  !================================================================================================================================
  !

  !>Sets up the output type for a fluid mechanics problem class.
  SUBROUTINE FluidMechanics_PostSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("FluidMechanics_PostSolve",err,error,*999)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)    
    IF(.NOT.ALLOCATED(problem%specification)) &
      & CALL FlagError("Problem specification is not allocated.",err,error,*999)
    IF(SIZE(problem%specification,1)<2) &
      & CALL FlagError("Problem specification must have at least two entries for a fluid mechanics problem.",err,error,*999)
      
    SELECT CASE(problem%specification(2))
    CASE(PROBLEM_STOKES_EQUATION_TYPE)
      CALL Stokes_PostSolve(solver,err,error,*999)
    CASE(PROBLEM_NAVIER_STOKES_EQUATION_TYPE)
      CALL NavierStokes_PostSolve(solver,err,error,*999)
    CASE(PROBLEM_DARCY_EQUATION_TYPE)
      CALL Darcy_PostSolve(solver,err,error,*999)
    CASE(PROBLEM_POISEUILLE_EQUATION_TYPE)
      CALL Poiseuille_PostSolve(solver,err,error,*999)
    CASE(PROBLEM_BURGERS_EQUATION_TYPE)
      CALL Burgers_PostSolve(solver,err,error,*999)
    CASE DEFAULT
      localError="Problem type "//TRIM(NumberToVString(problem%specification(2),"*",err,error))// &
        & " is not valid for a fluid mechanics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
      
    EXITS("FluidMechanics_PostSolve")
    RETURN
999 ERRORSEXITS("FluidMechanics_PostSolve",err,error)
    RETURN 1
    
  END SUBROUTINE FluidMechanics_PostSolve

  !
  !================================================================================================================================


  !>Sets up the output type for a fluid mechanics problem class.
  SUBROUTINE FluidMechanics_PreSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("FluidMechanics_PreSolve",err,error,*999)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)    
    IF(.NOT.ALLOCATED(problem%specification)) &
      & CALL FlagError("Problem specification is not allocated.",err,error,*999)
    IF(SIZE(problem%specification,1)<2) &
      & CALL FlagError("Problem specification must have at least two entries for a fluid mechanics problem.",err,error,*999)
    
    SELECT CASE(problem%specification(2))
    CASE(PROBLEM_STOKES_EQUATION_TYPE)
      CALL Stokes_PreSolve(solver,err,error,*999)
    CASE(PROBLEM_NAVIER_STOKES_EQUATION_TYPE)
      CALL NavierStokes_PreSolve(SOLVER,err,error,*999)
    CASE(PROBLEM_DARCY_EQUATION_TYPE)
      CALL Darcy_PreSolve(solver,err,error,*999)
    CASE(PROBLEM_POISEUILLE_EQUATION_TYPE)
      CALL Poiseuille_PreSolve(solver,err,error,*999)
    CASE(PROBLEM_BURGERS_EQUATION_TYPE)
      CALL Burgers_PreSolve(SOLVER,err,error,*999)
    CASE DEFAULT
      localError="Problem type "//TRIM(NumberToVString(problem%specification(2),"*",err,error))// &
        & " is not valid for a fluid mechanics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("FluidMechanics_PreSolve")
    RETURN
999 ERRORSEXITS("FluidMechanics_PreSolve",err,error)
    RETURN 1
    
  END SUBROUTINE FluidMechanics_PreSolve

  !
  !================================================================================================================================
  !

  !>Executes before each loop of a control loop, ie before each time step for a time loop
  SUBROUTINE FluidMechanics_PreLoop(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("FluidMechanics_PreLoop",err,error,*999)

    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("ControlLoop is not associated.",err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    
    SELECT CASE(controlLoop%loopType)
    CASE(CONTROL_TIME_LOOP_TYPE)
      IF(.NOT.ALLOCATED(problem%specification)) &
        & CALL FlagError("Problem specification is not allocated.",err,error,*999)
      IF(SIZE(problem%specification,1)<2) &
        & CALL FlagError("Problem specification must have at least two entries for a fluid mechanics problem.",err,error,*999)
      
      SELECT CASE(problem%specification(2))
      CASE(PROBLEM_STOKES_EQUATION_TYPE)
        !do nothing
      CASE(PROBLEM_NAVIER_STOKES_EQUATION_TYPE)
        !do nothing
      CASE(PROBLEM_DARCY_EQUATION_TYPE)
        CALL Darcy_PreLoop(controlLoop,err,error,*999)
      CASE(PROBLEM_POISEUILLE_EQUATION_TYPE)
        !do nothing
      CASE(PROBLEM_BURGERS_EQUATION_TYPE)
        !do nothing
      CASE DEFAULT
        localError="Problem type "//TRIM(NumberToVString(problem%specification(2),"*",err,error))// &
          & " is not valid for a fluid mechanics problem class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      !do nothing
    END SELECT

    EXITS("FluidMechanics_PreLoop")
    RETURN
999 ERRORSEXITS("FluidMechanics_PreLoop",err,error)
    RETURN 1
    
  END SUBROUTINE FluidMechanics_PreLoop

  !
  !================================================================================================================================
  !

  !>Executes after each loop of a control loop, ie after each time step for a time loop
  SUBROUTINE FluidMechanics_PostLoop(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("FluidMechanics_PostLoop",err,error,*999)

    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("ControlLoop is not associated.",err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    
    SELECT CASE(controlLoop%loopType)
    CASE(CONTROL_TIME_LOOP_TYPE)
      IF(.NOT.ALLOCATED(problem%specification)) &
        & CALL FlagError("Problem specification is not allocated.",err,error,*999)
      IF(SIZE(problem%specification,1)<2) &
        & CALL FlagError("Problem specification must have at least two entries for a fluid mechanics problem.",err,error,*999)
      
      SELECT CASE(problem%specification(2))
      CASE(PROBLEM_STOKES_EQUATION_TYPE)
        !do nothing
      CASE(PROBLEM_NAVIER_STOKES_EQUATION_TYPE)
        CALL NavierStokes_PostLoop(controlLoop,err,error,*999)
      CASE(PROBLEM_DARCY_EQUATION_TYPE)
        !do nothing
      CASE(PROBLEM_POISEUILLE_EQUATION_TYPE)
        !do nothing
      CASE(PROBLEM_BURGERS_EQUATION_TYPE)
        !do nothing
      CASE DEFAULT
        localError="The second problem specification of "// &
          & TRIM(NumberToVString(problem%specification(2),"*",err,error))// &
          & " is not valid for a fluid mechanics problem."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      !do nothing
    END SELECT

    EXITS("FluidMechanics_PostLoop")
    RETURN
999 ERRORSEXITS("FluidMechanics_PostLoop",err,error)
    RETURN 1
    
  END SUBROUTINE FluidMechanics_PostLoop

  !
  !================================================================================================================================
  !

  !>Pre-residual steps for an fluid mechanics class finite element equation set.
  SUBROUTINE FluidMechanics_FiniteElementPreResidualEvaluate(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer the equations set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("FluidMechanics_FiniteElementPreResidualEvaluate",err,error,*999)
    
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated",err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) &
      & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<2) &
      & CALL FlagError("Equations set specification must have at least two entries for a fluid mechanics class equations set.", &
      & err,error,*999)
    
    SELECT CASE(equationsSet%specification(2))
    CASE(EQUATIONS_SET_STOKES_EQUATION_TYPE)
      ! Do nothing
    CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE)
      CALL NavierStokes_FiniteElementPreResidualEvaluate(equationsSet,err,error,*999)
    CASE(EQUATIONS_SET_DARCY_EQUATION_TYPE)
      ! Do nothing
    CASE(EQUATIONS_SET_DARCY_PRESSURE_EQUATION_TYPE)
      ! Do nothing
    CASE(EQUATIONS_SET_POISEUILLE_EQUATION_TYPE)
      ! Do nothing
    CASE(EQUATIONS_SET_BURGERS_EQUATION_TYPE)
      ! Do nothing
    CASE(EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE)
      ! Do nothing
    CASE DEFAULT
      localError="The second equations set specificaiton of "// &
        & TRIM(NumberToVString(equationsSet%specification(2),"*",err,error))// &
        & " is not valid for a fluid mechanics equation set."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("FluidMechanics_FiniteElementPreResidualEvaluate")
    RETURN
999 ERRORS("FluidMechanics_FiniteElementPreResidualEvaluate",err,error)
    EXITS("FluidMechanics_FiniteElementPreResidualEvaluate")
    RETURN 1
    
  END SUBROUTINE FluidMechanics_FiniteElementPreResidualEvaluate

  !
  !================================================================================================================================
  !

END MODULE FluidMechanicsRoutines

