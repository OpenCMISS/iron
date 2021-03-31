!> \file
!> \author Chris Bradley
!> \brief This module contains all boundary condition access method routines.
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

!> \addtogroup OpenCMISS_BoundaryConditions OpenCMISS::Iron::BoundaryConditions
!> This module contains all boundary condition access method routines.
MODULE BoundaryConditionAccessRoutines
  
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

  !> \addtogroup BoundaryConditionsRoutines_Constants OpenCMISS::Iron::BoundaryConditions::Constants
  !> \brief Boundary conditions constants.
  !>@{
  !> \addtogroup BoundaryConditionsRoutines_RowTypes OpenCMISS::Iron::BoundaryConditions::Constants::RowTypes
  !> \brief Boundary conditions type for rows.
  !>@{
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FREE_ROW=0 !<The row is free. \see BoundaryConditionsRoutines_RowTypes,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_DIRICHLET_ROW=1 !<The row has a Dirichlet boundary condition. \see BoundaryConditionsRoutines_RowTypes,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_POINT_NEUMANN_ROW=2 !<The row has a point Neumann boundary condition. \see BoundaryConditionsRoutines_RowTypes,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_INTEGRATED_NEUMANN_ROW=3 !<The row has an integrated Neumann boundary condition. \see BoundaryConditionsRoutines_RowTypes,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_NEUMANN_ROW=3 !<The row has a Neumann boundary condition. \see BoundaryConditionsRoutines_RowTypes,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_ROBIN_ROW=4 !<The row has a Robin boundary condition. \see BoundaryConditionsRoutines_RowTypes,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_CAUCHY_ROW=5 !<The row has a Cauchy boundary condition. \see BoundaryConditionsRoutines_RowTypes,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_CONSTRAINED_ROW=6 !<The row has a constraint boundary condition. \see BoundaryConditionsRoutines_RowTypes,BoundaryConditionsRoutines
  !>@}
  
  !> \addtogroup BoundaryConditionsRoutines_DOFTypes OpenCMISS::Iron:::BoundaryConditions::Constants::DOFTypes
  !> \brief DOF type for boundary conditions.
  !>@{
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_DOF_FREE=0 !<The dof is free. \see BoundaryConditionsRoutines_DOFTypes,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_DOF_FIXED=1 !<The dof is fixed as a boundary condition. \see BoundaryConditionsRoutines_DOFTypes,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_DOF_INCREMENTED=2 !<The dof is set as a mixed boundary condition. \see BoundaryConditionsRoutines_DOFTypes,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_DOF_MIXED=3 !<The dof is set as a mixed boundary condition. \see BoundaryConditionsRoutines_DOFTypes,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_DOF_CONSTRAINED=4 !<The dof is constrained to be a linear combination of other DOFs. \see BoundaryConditionsRoutines_DOFTypes,BoundaryConditionsRoutines
  !>@}
  
  !> \addtogroup BoundaryConditionsRoutines_BoundaryConditionsTypes OpenCMISS::Iron:::BoundaryConditions::Constants::BoundaryConditionTypes
  !> \brief Boundary conditions types. These may be specific to a particular equation type and the solver routines should not need to use these.
  !>@{
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_NONE=0 !<The DOF has no boundary condition. \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED=1 !<The DOF is fixed as a boundary condition. \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED_INLET=2 !<The dof is fixed as a boundary condition. \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED_OUTLET=3 !<The dof is fixed as a boundary condition. \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED_WALL=4 !<The dof is fixed as a boundary condition. \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_MOVED_WALL=5 !<The dof is fixed as a boundary condition. \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FREE_WALL=6 !<The dof is fixed as a boundary condition. \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_NEUMANN_POINT=8 !<The dof is set to a Neumann point boundary condition. \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_NEUMANN_INTEGRATED=9 !<The dof is set to a Neumann integrated boundary condition. \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_DIRICHLET=10 !<The dof is set to a Dirichlet boundary condition. \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_CAUCHY=11 !<The dof is set to a Cauchy boundary condition. \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_ROBIN=12 !<The dof is set to a Robin boundary condition. \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED_INCREMENTED=13 !<The dof is a fixed boundary condition, to be used with load increment loop. \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_PRESSURE=14 !<The dof is a surface pressure boundary condition. \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_PRESSURE_INCREMENTED=15 !<The dof is a surface pressure boundary condition, to be used with load increment loop. \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED=17 !<The dof is fixed as a boundary condition, to be used with load increment loop. \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_CORRECTION_MASS_INCREASE=18 !<The dof is fixed as a boundary condition, to be used with load increment loop. \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_IMPERMEABLE_WALL=19 !<The dof is set such that (via penalty formulation): velocity * normal = 0. \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_NEUMANN_INTEGRATED_ONLY=20 !<A Neumann integrated boundary condition, and no point values will be integrated over a face or line that includes this dof. \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_LINEAR_CONSTRAINT=21 !<The dof is constrained to be a linear combination of other DOFs. \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED=22!<A Neumann point boundary condition that is incremented inside a load increment control loop. \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED_FITTED=23 !<The dof is fixed as a boundary condition to be updated from fitting data \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED_NONREFLECTING=24 !<The dof is fixed and set to a non-reflecting type for 1D wave propagation problems. \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED_CELLML=25 !<The dof is fixed and set to values specified based on the coupled CellML solution at the dof. \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED_STREE=26 !<The dof is fixed and set to values specified based on the transmission line theory at the dof. \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_COUPLING_FLOW=27 !<The dof is fixed and set to values specified based on a coupled flow rate at the dof. \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_COUPLING_STRESS=28 !<The dof is fixed and set to values specified based on a coupled stress at the dof. \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED_PRESSURE=29 !<The dof is a fixed pressure boundary condition. \see BoundaryConditionsRoutines_BoundaryConditionsTypes,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: MAX_BOUNDARY_CONDITION_NUMBER=29 !The maximum boundary condition type identifier, used for allocating an array with an entry for each type
  !>@}
  
  !> \addtogroup BoundaryConditionsRoutines_SparsityTypes OpenCMISS::Iron::BoundaryConditions::Constants::SparsityTypes
  !> \brief Storage type for matrices used by boundary conditions.
  !>@{
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_SPARSE_MATRICES=1 !<The matrices are stored as sparse matrices.
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FULL_MATRICES=2 !<The matrices are stored as full matrices.
  !>@}
  
  !>@}

  !Module types

  !Module variables

  !Interfaces

  PUBLIC BOUNDARY_CONDITION_FREE_ROW,BOUNDARY_CONDITION_DIRICHLET_ROW,BOUNDARY_CONDITION_POINT_NEUMANN_ROW, &
    & BOUNDARY_CONDITION_INTEGRATED_NEUMANN_ROW,BOUNDARY_CONDITION_NEUMANN_ROW,BOUNDARY_CONDITION_ROBIN_ROW, &
    & BOUNDARY_CONDITION_CAUCHY_ROW,BOUNDARY_CONDITION_CONSTRAINED_ROW
  
  PUBLIC BOUNDARY_CONDITION_DOF_FREE,BOUNDARY_CONDITION_DOF_FIXED,BOUNDARY_CONDITION_DOF_INCREMENTED, &
    & BOUNDARY_CONDITION_DOF_MIXED,BOUNDARY_CONDITION_DOF_CONSTRAINED

  PUBLIC BOUNDARY_CONDITION_NONE,BOUNDARY_CONDITION_FIXED,BOUNDARY_CONDITION_FIXED_INLET, &
    & BOUNDARY_CONDITION_FIXED_OUTLET,BOUNDARY_CONDITION_FIXED_WALL,BOUNDARY_CONDITION_MOVED_WALL,BOUNDARY_CONDITION_FREE_WALL, &
    & BOUNDARY_CONDITION_NEUMANN_POINT,BOUNDARY_CONDITION_NEUMANN_INTEGRATED,BOUNDARY_CONDITION_DIRICHLET, &
    & BOUNDARY_CONDITION_CAUCHY,BOUNDARY_CONDITION_ROBIN,BOUNDARY_CONDITION_FIXED_INCREMENTED,BOUNDARY_CONDITION_PRESSURE, &
    & BOUNDARY_CONDITION_PRESSURE_INCREMENTED,BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED, &
    & BOUNDARY_CONDITION_CORRECTION_MASS_INCREASE,BOUNDARY_CONDITION_IMPERMEABLE_WALL,BOUNDARY_CONDITION_NEUMANN_INTEGRATED_ONLY, &
    & BOUNDARY_CONDITION_LINEAR_CONSTRAINT,BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED,BOUNDARY_CONDITION_FIXED_FITTED, &
    & BOUNDARY_CONDITION_FIXED_NONREFLECTING,BOUNDARY_CONDITION_FIXED_CELLML,BOUNDARY_CONDITION_FIXED_STREE, &
    & BOUNDARY_CONDITION_COUPLING_FLOW,BOUNDARY_CONDITION_COUPLING_STRESS,BOUNDARY_CONDITION_FIXED_PRESSURE

  PUBLIC MAX_BOUNDARY_CONDITION_NUMBER

  PUBLIC BOUNDARY_CONDITION_SPARSE_MATRICES,BOUNDARY_CONDITION_FULL_MATRICES

  PUBLIC BoundaryConditions_AssertIsFinished,BoundaryConditions_AssertNotFinished

  PUBLIC BoundaryConditions_NumberOfRowVariablesGet

  PUBLIC BoundaryConditions_NumberOfVariablesGet

  PUBLIC BoundaryConditions_SolverEquationsGet

  PUBLIC BoundaryConditions_RowVariableIndexGet

  PUBLIC BoundaryConditions_VariableIndexGet

  PUBLIC BoundaryConditionsDirichlet_DirichletDOFIndexGet
  
  PUBLIC BoundaryConditionsDirichlet_DynamicSparsityIndicesGet
  
  PUBLIC BoundaryConditionsDirichlet_LinearSparsityIndicesGet

  PUBLIC BoundaryConditionsNeumann_NeumannDOFIndexGet
  
  PUBLIC BoundaryConditionsNeumann_PointDOFMappingGet
  
  PUBLIC BoundaryConditionsNeumann_PointDOFValuesGet

  PUBLIC BoundaryConditionsPressureInc_PressureIncDOFIndexGet
  
  PUBLIC BoundaryConditionsRowVariable_BoundaryConditionsGet

  PUBLIC BoundaryConditionsRowVariable_LHSVariableExists

  PUBLIC BoundaryConditionsRowVariable_LHSVariableGet

  PUBLIC BoundaryConditionsRowVariable_NumberOfDOFsGet

  PUBLIC BoundaryConditionsRowVariable_RHSVariableExists

  PUBLIC BoundaryConditionsRowVariable_RowConditionTypeGet

  PUBLIC BoundaryConditionsRowVariable_TotalNumberOfDOFsGet

  PUBLIC BoundaryConditionsVariable_BoundaryConditionsGet

  PUBLIC BoundaryConditionsVariable_ConditionTypeGet

  PUBLIC BoundaryConditionsVariable_DirichletConditionsExists

  PUBLIC BoundaryConditionsVariable_DirichletConditionsGet

  PUBLIC BoundaryConditionsVariable_DOFConstraintsExists
  
  PUBLIC BoundaryConditionsVariable_DOFConstraintsGet

  PUBLIC BoundaryConditionsVariable_DOFCountGet

  PUBLIC BoundaryConditionsVariable_DOFTypeGet

  PUBLIC BoundaryConditionsVariable_FieldVariableGet

  PUBLIC BoundaryConditionsVariable_FieldVariableTypeGet

  PUBLIC BoundaryConditionsVariable_NeumannConditionsExists

  PUBLIC BoundaryConditionsVariable_NeumannConditionsGet

  PUBLIC BoundaryConditionsVariable_NumberOfDirichletConditionsGet

  PUBLIC BoundaryConditionsVariable_PressureIncConditionsExists

  PUBLIC BoundaryConditionsVariable_PressureIncConditionsGet

  PUBLIC BoundaryConditionsVariableDOFConstraints_ConstraintExists
  
  PUBLIC BoundaryConditionsVariableDOFConstraints_ConstraintGet
  
  PUBLIC BoundaryConditionsVariableDOFConstraints_DOFCouplingExists
 
  PUBLIC BoundaryConditionsVariableDOFConstraints_DOFCouplingGet
 
  PUBLIC BoundaryConditionsVariableDOFConstraints_NumberOfConstraintsGet
 
  PUBLIC BoundaryConditionsVariableDOFConstraints_NumberOfDOFsGet
 
CONTAINS

  !
  !=================================================================================================================================
  !

  !>Assert that a boundary conditions has been finished
  SUBROUTINE BoundaryConditions_AssertIsFinished(boundaryConditions,err,error,*)

    !Argument Variables
    TYPE(BoundaryConditionsType), POINTER, INTENT(IN) :: boundaryConditions !<The boundary conditions to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("BoundaryConditions_AssertIsFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions is not associated.",err,error,*999)
#endif    

    IF(.NOT.boundaryConditions%boundaryConditionsFinished) &
      & CALL FlagError("Boundary conditions has not been finished.",err,error,*999)
    
    EXITS("BoundaryConditions_AssertIsFinished")
    RETURN
999 ERRORSEXITS("BoundaryConditions_AssertIsFinished",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_AssertIsFinished

  !
  !=================================================================================================================================
  !

  !>Assert that a boundary conditions has not been finished
  SUBROUTINE BoundaryConditions_AssertNotFinished(boundaryConditions,err,error,*)

    !Argument Variables
    TYPE(BoundaryConditionsType), POINTER, INTENT(IN) :: boundaryConditions !<The boundary conditions to assert the finished status for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("BoundaryConditions_AssertNotFinished",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions is not associated.",err,error,*999)
#endif
    
    IF(boundaryConditions%boundaryConditionsFinished) &
      & CALL FlagError("Boundary conditions has already been finished.",err,error,*999)
    
    EXITS("BoundaryConditions_AssertNotFinished")
    RETURN
999 ERRORSEXITS("BoundaryConditions_AssertNotFinished",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_AssertNotFinished
  
  !
  !================================================================================================================================
  !

  !>Gets the number of boundary conditions row variables for boundary conditions.
  SUBROUTINE BoundaryConditions_NumberOfRowVariablesGet(boundaryConditions,numberOfRowVariables,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to get the number of boundary conditions row variables for
    INTEGER(INTG), INTENT(OUT) :: numberOfRowVariables !<On exit, the number of boundary conditions row variables in the boundary conditions
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("BoundaryConditions_NumberOfRowVariablesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions is not associated.",err,error,*999)
#endif    

    numberOfRowVariables=boundaryConditions%numberOfBoundaryConditionsRowVariables

    EXITS("BoundaryConditions_NumberOfRowVariablesGet")
    RETURN
999 ERRORSEXITS("BoundaryConditions_NumberOfRowVariablesGet",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_NumberOfRowVariablesGet

  !
  !================================================================================================================================
  !

  !>Gets the number of boundary conditions variables for boundary conditions.
  SUBROUTINE BoundaryConditions_NumberOfVariablesGet(boundaryConditions,numberOfVariables,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to get the number of boundary conditions variables for
    INTEGER(INTG), INTENT(OUT) :: numberOfVariables !<On exit, the number of boundary conditions variables in the boundary conditions
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("BoundaryConditions_NumberOfVariablesGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions is not associated.",err,error,*999)
#endif    

    numberOfVariables=boundaryConditions%numberOfBoundaryConditionsVariables

    EXITS("BoundaryConditions_NumberOfVariablesGet")
    RETURN
999 ERRORSEXITS("BoundaryConditions_NumberOfVariablesGet",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_NumberOfVariablesGet

  !
  !================================================================================================================================
  !

  !>Gets the solver equations for boundary conditions.
  SUBROUTINE BoundaryConditions_SolverEquationsGet(boundaryConditions,solverEquations,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to get the solver equations for
    TYPE(SolverEquationsType), POINTER :: solverEquations !<On exit, a pointer to the solver equations in the specified boundary conditions. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("BoundaryConditions_SolverEquationsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions is not associated.",err,error,*999)
#endif    

    solverEquations=>boundaryConditions%solverEquations

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(solverEquations)) &
      & CALL FlagError("Solver equations is not associated for the boundary conditions.",err,error,*999)
#endif    
       
    EXITS("BoundaryConditions_SolverEquationsGet")
    RETURN
999 NULLIFY(solverEquations)
998 ERRORSEXITS("BoundaryConditions_SolverEquationsGet",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_SolverEquationsGet

  !
  !================================================================================================================================
  !

  !>Gets the specificed boundary conditions row variable for a row variable index of boundary conditions.
  SUBROUTINE BoundaryConditions_RowVariableIndexGet(boundaryConditions,rowVariableIdx,boundaryConditionsRowVariable,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to get the boundary conditions row variable for
    INTEGER(INTG), INTENT(IN) :: rowVariableIdx !<The row variable index to get the boundary conditions row variable for
    TYPE(BoundaryConditionsRowVariableType), POINTER :: boundaryConditionsRowVariable !<On exit, a pointer to the specificed boundary conditions row variable in the specified boundary conditions. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif
 
    ENTERS("BoundaryConditions_RowVariableIndexGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(boundaryConditionsRowVariable)) &
      & CALL FlagError("Boundary conditions row variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions is not associated.",err,error,*999)
    IF(rowVariableIdx<1.OR.rowVariableIdx>boundaryConditions%numberOfBoundaryConditionsRowVariables) THEN
      localError="The specificed boundary conditions row variable index of "// &
        & TRIM(NumberToVString(rowVariableIdx,"*",err,error))//" is invalid. The row variable index should be >= 1 and <= "// &
        & TRIM(NumberToVString(boundaryConditions%numberOfBoundaryConditionsRowVariables,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(boundaryConditions%boundaryConditionsRowVariables)) &
      & CALL FlagError("Boundary conditions row variables is not allocated for the boundary conditions.",err,error,*999)
#endif    

    boundaryConditionsRowVariable=>boundaryConditions%boundaryConditionsRowVariables(rowVariableIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(boundaryConditionsRowVariable)) THEN
      localError="The boundary conditions row variable is not associated for row variable index "// &
        & TRIM(NumberToVString(rowVariableIdx,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
    EXITS("BoundaryConditions_RowVariableIndexGet")
    RETURN
999 NULLIFY(boundaryConditionsRowVariable)
998 ERRORSEXITS("BoundaryConditions_RowVariableIndexGet",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_RowVariableIndexGet

  !
  !================================================================================================================================
  !

  !>Gets the specificed boundary conditions variable for a variable index of boundary conditions.
  SUBROUTINE BoundaryConditions_VariableIndexGet(boundaryConditions,variableIdx,boundaryConditionsVariable,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to get the boundary conditions variable for
    INTEGER(INTG), INTENT(IN) :: variableIdx !<The variable index to get the boundary conditions variable for
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable !<On exit, a pointer to the specificed boundary conditions variable in the specified boundary conditions. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif
 
    ENTERS("BoundaryConditions_VariableIndexGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(boundaryConditionsVariable)) CALL FlagError("Boundary conditions variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions is not associated.",err,error,*999)
    IF(variableIdx<1.OR.variableIdx>boundaryConditions%numberOfBoundaryConditionsVariables) THEN
      localError="The specificed boundary conditions variable index of "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
        & " is invalid. The variable index should be >= 1 and <= "// &
        & TRIM(NumberToVString(boundaryConditions%numberOfBoundaryConditionsVariables,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(boundaryConditions%boundaryConditionsVariables)) &
      & CALL FlagError("Boundary conditions variables is not allocated for the boundary conditions.",err,error,*999)
#endif    

    boundaryConditionsVariable=>boundaryConditions%boundaryConditionsVariables(variableIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(boundaryConditionsVariable)) THEN
      localError="The boundary conditions variable is not associated for variable index "// &
        & TRIM(NumberToVString(variableIdx,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
    EXITS("BoundaryConditions_VariableIndexGet")
    RETURN
999 NULLIFY(boundaryConditionsVariable)
998 ERRORSEXITS("BoundaryConditions_VariableIndexGet",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_VariableIndexGet

  !
  !================================================================================================================================
  !

  !>Gets the Dirichlet DOF index for the Dirichlet boundary conditions.
  SUBROUTINE BoundaryConditionsDirichlet_DirichletDOFIndexGet(dirichletBoundaryConditions,dirichletIdx,dirichletDOFIndex, &
    & err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsDirichletType), POINTER :: dirichletBoundaryConditions !<A pointer to the boundary conditions variable to get the Dirichlet DOF index for
    INTEGER(INTG), INTENT(IN) :: dirichletIdx !<The Dirichlet index to get the Dirichlet DOF index for
    INTEGER(INTG), INTENT(OUT) :: dirichletDOFIndex !<On exit, the Dirichlet DOF index.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif
 
    ENTERS("BoundaryConditionsDirichlet_DirichletDOFIndexGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(dirichletBoundaryConditions)) &
      & CALL FlagError("Dirichlet boundary conditions is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(dirichletBoundaryConditions%dirichletDOFIndices)) &
      & CALL FlagError("The Dirichlet DOF indices are not allocated for the Dirichlet boundary conditions.",err,error,*999)
    IF(dirichletIdx<1.OR.dirichletIdx>SIZE(dirichletBoundaryConditions%dirichletDOFIndices,1)) THEN
      localError="The specified Dirichlet index of "//TRIM(NumberToVString(dirichletIdx,"*",err,error))// &
        & " is invalid. The Dirichlet index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(dirichletBoundaryConditions%dirichletDOFIndices,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    dirichletDOFIndex=dirichletBoundaryConditions%dirichletDOFIndices(dirichletIdx)
      
    EXITS("BoundaryConditionsDirichlet_DirichletDOFIndexGet")
    RETURN
999 ERRORS("BoundaryConditionsDirichlet_DirichletDOFIndexGet",err,error)
    EXITS("BoundaryConditionsDirichlet_DirichletDOFIndexGet")
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsDirichlet_DirichletDOFIndexGet

  !
  !================================================================================================================================
  !

  !>Gets the dynamic matrices sparsity indices for the Dirichlet boundary conditions.
  SUBROUTINE BoundaryConditionsDirichlet_DynamicSparsityIndicesGet(dirichletBoundaryConditions,equationsMatrixIdx, &
    & equationsSetIdx,dynamicSparsityIndices,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsDirichletType), POINTER :: dirichletBoundaryConditions !<A pointer to the boundary conditions variable to get the dynamic sparsity indices for
    INTEGER(INTG), INTENT(IN) :: equationsMatrixIdx !<The equations matrix index to get the dynamic sparsity indices for
    INTEGER(INTG), INTENT(IN) :: equationsSetIdx !<The equations set index to get the dynamic sparsity indices for.
    TYPE(BoundaryConditionsSparsityIndicesType), POINTER :: dynamicSparsityIndices !<On exit, a pointer to the dynamic sparsity indices in the specified Dirichlet boundary conditions. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif
 
    ENTERS("BoundaryConditionsDirichlet_DynamicSparsityIndicesGet",err,error,*998)

#ifdef WITH_PRECHECKS
    IF(ASSOCIATED(dynamicSparsityIndices)) CALL FlagError("Dynamic sparsity indices are already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(dirichletBoundaryConditions)) &
      & CALL FlagError("Dirichlet boundary conditions is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(dirichletBoundaryConditions%dynamicSparsityIndices)) &
      & CALL FlagError("The dynamic sparsity indices are not allocated for the Dirichlet boundary conditions.",err,error,*999)
    IF(equationsMatrixIdx<1.OR.equationsMatrixIdx>SIZE(dirichletBoundaryConditions%dynamicSparsityIndices,2)) THEN
      localError="The specified equations matrix index of "//TRIM(NumberToVString(equationsMatrixIdx,"*",err,error))// &
        & " is invalid. The equations matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(dirichletBoundaryConditions%dynamicSparsityIndices,2),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(equationsSetIdx<1.OR.equationsSetIdx>SIZE(dirichletBoundaryConditions%dynamicSparsityIndices,1)) THEN
      localError="The specified equations set index of "//TRIM(NumberToVString(equationsSetIdx,"*",err,error))// &
        & " is invalid. The equations set index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(dirichletBoundaryConditions%dynamicSparsityIndices,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    dynamicSparsityIndices=>dirichletBoundaryConditions%dynamicSparsityIndices(equationsSetIdx,equationsMatrixIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(dynamicSparsityIndices)) THEN
      localError="The dynamic sparsity indices are not associated for equations matrix index "// &
        & TRIM(NumberToVString(equationsMatrixIdx,"*",err,error))//" of equations set index "// &
        & TRIM(NumberToVString(equationsSetIdx,"*",err,error))//" in the Dirichlet boundary conditions."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
    EXITS("BoundaryConditionsDirichlet_DynamicSparsityIndicesGet")
    RETURN
999 NULLIFY(dynamicSparsityIndices)
998 ERRORS("BoundaryConditionsDirichlet_DynamicSparsityIndicesGet",err,error)
    EXITS("BoundaryConditionsDirichlet_DynamicSparsityIndicesGet")
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsDirichlet_DynamicSparsityIndicesGet

  !
  !================================================================================================================================
  !

  !>Gets the linear matrices sparsity indices for the Dirichlet boundary conditions.
  SUBROUTINE BoundaryConditionsDirichlet_LinearSparsityIndicesGet(dirichletBoundaryConditions,equationsMatrixIdx, &
    & equationsSetIdx,linearSparsityIndices,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsDirichletType), POINTER :: dirichletBoundaryConditions !<A pointer to the boundary conditions variable to get the linear sparsity indices for
    INTEGER(INTG), INTENT(IN) :: equationsMatrixIdx !<The equations matrix index to get the linear sparsity indices for
    INTEGER(INTG), INTENT(IN) :: equationsSetIdx !<The equations set index to get the linear sparsity indices for.
    TYPE(BoundaryConditionsSparsityIndicesType), POINTER :: linearSparsityIndices !<On exit, a pointer to the linear sparsity indices in the specified Dirichlet boundary conditions. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_CHECKS
    TYPE(VARYING_STRING) :: localError
#endif
 
    ENTERS("BoundaryConditionsDirichlet_LinearSparsityIndicesGet",err,error,*998)

#ifdef WITH_PRECHECKS
    IF(ASSOCIATED(linearSparsityIndices)) CALL FlagError("Linear sparsity indices are already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(dirichletBoundaryConditions)) &
      & CALL FlagError("Dirichlet boundary conditions is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(dirichletBoundaryConditions%linearSparsityIndices)) &
      & CALL FlagError("The linear sparsity indices are not allocated for the Dirichlet boundary conditions.",err,error,*999)
    IF(equationsMatrixIdx<1.OR.equationsMatrixIdx>SIZE(dirichletBoundaryConditions%linearSparsityIndices,2)) THEN
      localError="The specified equations matrix index of "//TRIM(NumberToVString(equationsMatrixIdx,"*",err,error))// &
        & " is invalid. The equations matrix index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(dirichletBoundaryConditions%linearSparsityIndices,2),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(equationsSetIdx<1.OR.equationsSetIdx>SIZE(dirichletBoundaryConditions%linearSparsityIndices,1)) THEN
      localError="The specified equations set index of "//TRIM(NumberToVString(equationsSetIdx,"*",err,error))// &
        & " is invalid. The equations set index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(dirichletBoundaryConditions%linearSparsityIndices,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    linearSparsityIndices=>dirichletBoundaryConditions%linearSparsityIndices(equationsSetIdx,equationsMatrixIdx)%ptr

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(linearSparsityIndices)) THEN
      localError="The linear sparsity indices are not associated for equations matrix index "// &
        & TRIM(NumberToVString(equationsMatrixIdx,"*",err,error))//" of equations set index "// &
        & TRIM(NumberToVString(equationsSetIdx,"*",err,error))//" in the Dirichlet boundary conditions."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    
       
    EXITS("BoundaryConditionsDirichlet_LinearSparsityIndicesGet")
    RETURN
999 NULLIFY(linearSparsityIndices)
998 ERRORS("BoundaryConditionsDirichlet_LinearSparsityIndicesGet",err,error)
    EXITS("BoundaryConditionsDirichlet_LinearSparsityIndicesGet")
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsDirichlet_LinearSparsityIndicesGet

  !
  !================================================================================================================================
  !

  !>Gets the Neumann DOF index for the Neumann boundary conditions.
  SUBROUTINE BoundaryConditionsNeumann_NeumannDOFIndexGet(neumannBoundaryConditions,neumannIdx,neumannDOFIndex,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsNeumannType), POINTER :: neumannBoundaryConditions !<A pointer to the Neumann boundary conditions to get the Neumann DOF index for
    INTEGER(INTG), INTENT(IN) :: neumannIdx !<The Neumann index to get the Neumann DOF index for
    INTEGER(INTG), INTENT(OUT) :: neumannDOFIndex !<On exit, the Neumann DOF index.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif
 
    ENTERS("BoundaryConditionsNeumann_NeumannDOFIndexGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(neumannBoundaryConditions)) CALL FlagError("Neumann boundary conditions is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(neumannBoundaryConditions%setDOFs)) &
      & CALL FlagError("The Neumann DOF indices are not allocated for the Neumann boundary conditions.",err,error,*999)
    IF(neumannIdx<1.OR.neumannIdx>SIZE(neumannBoundaryConditions%setDOFs,1)) THEN
      localError="The specified Neumann index of "//TRIM(NumberToVString(neumannIdx,"*",err,error))// &
        & " is invalid. The Neumann index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(neumannBoundaryConditions%setDOFs,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    neumannDOFIndex=neumannBoundaryConditions%setDOFs(neumannIdx)
      
    EXITS("BoundaryConditionsNeumann_NeumannDOFIndexGet")
    RETURN
999 ERRORS("BoundaryConditionsNeumann_NeumannDOFIndexGet",err,error)
    EXITS("BoundaryConditionsNeumann_NeumannDOFIndexGet")
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsNeumann_NeumannDOFIndexGet

  !
  !================================================================================================================================
  !

  !>Gets the Neumann point DOF mapping for the Neumann boundary conditions.
  SUBROUTINE BoundaryConditionsNeumann_PointDOFMappingGet(neumannBoundaryConditions,pointDOFMapping,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsNeumannType), POINTER :: neumannBoundaryConditions !<A pointer to the Neumann boundary conditions to get the Neumann point DOF mappings for
    TYPE(DomainMappingType), POINTER :: pointDOFMapping !<On exit, a pointer to the point DOF mapping. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("BoundaryConditionsNeumann_PointDOFMappingGet",err,error,*998)

#ifdef WITH_PRECHECKS
    IF(ASSOCIATED(pointDOFMapping)) CALL FlagError("Point DOF mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(neumannBoundaryConditions)) CALL FlagError("Neumann boundary conditions is not associated.",err,error,*999)
#endif    

    pointDOFMapping=>neumannBoundaryConditions%pointDOFMapping

#ifdef WITH_POSTCHECKS
    IF(.NOT.ASSOCIATED(pointDOFMapping)) &
      & CALL FlagError("The point DOF mapping is not associated for the Neumann boundary conditions.",err,error,*999)
#endif
      
    EXITS("BoundaryConditionsNeumann_PointDOFMappingGet")
    RETURN
999 NULLIFY(pointDOFMapping)
998 ERRORS("BoundaryConditionsNeumann_PointDOFMappingGet",err,error)
    EXITS("BoundaryConditionsNeumann_PointDOFMappingGet")
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsNeumann_PointDOFMappingGet

  !
  !================================================================================================================================
  !

  !>Gets the Neumann point DOF values for the Neumann boundary conditions.
  SUBROUTINE BoundaryConditionsNeumann_PointDOFValuesGet(neumannBoundaryConditions,pointDOFValues,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsNeumannType), POINTER :: neumannBoundaryConditions !<A pointer to the Neumann boundary conditions to get the Neumann point DOF values for
    TYPE(DistributedVectorType), POINTER :: pointDOFValues !<On exit, a pointer to the point DOF values. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("BoundaryConditionsNeumann_PointDOFValuesGet",err,error,*998)

#ifdef WITH_PRECHECKS
    IF(ASSOCIATED(pointDOFValues)) CALL FlagError("Point DOF values is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(neumannBoundaryConditions)) CALL FlagError("Neumann boundary conditions is not associated.",err,error,*999)
#endif    

    pointDOFValues=>neumannBoundaryConditions%pointValues

#ifdef WITH_POSTCHECKS
    IF(.NOT.ASSOCIATED(pointDOFValues)) &
      & CALL FlagError("The point DOF values is not associated for the Neumann boundary conditions.",err,error,*999)
#endif
      
    EXITS("BoundaryConditionsNeumann_PointDOFValuesGet")
    RETURN
999 NULLIFY(pointDOFValues)
998 ERRORS("BoundaryConditionsNeumann_PointDOFValuesGet",err,error)
    EXITS("BoundaryConditionsNeumann_PointDOFValuesGet")
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsNeumann_PointDOFValuesGet

  !
  !================================================================================================================================
  !

  !>Gets the incremented pressure DOF index for the pressure incremented boundary conditions.
  SUBROUTINE BoundaryConditionsPressureInc_PressureIncDOFIndexGet(pressureIncBoundaryConditions,pressureIncIdx, &
    & pressureIncDOFIndex,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsPressureIncrementedType), POINTER :: pressureIncBoundaryConditions !<A pointer to the incremented pressure boundary conditions to get the incremented pressure DOF index for
    INTEGER(INTG), INTENT(IN) :: pressureIncIdx !<The incremented pressure index to get the incremented pressure DOF index for
    INTEGER(INTG), INTENT(OUT) :: pressureIncDOFIndex !<On exit, the incremented pressure DOF index.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif
 
    ENTERS("BoundaryConditionsPressureInc_PressureIncDOFIndexGet",err,error,*999)

#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(pressureIncBoundaryConditions)) &
      & CALL FlagError("Incremented pressure boundary conditions is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(pressureIncBoundaryConditions%pressureIncrementedDOFIndices)) &
      & CALL FlagError("The incremented pressure DOF indices are not allocated for the Neumann boundary conditions.",err,error,*999)
    IF(pressureIncIdx<1.OR.pressureIncIdx>SIZE(pressureIncBoundaryConditions%pressureIncrementedDOFIndices,1)) THEN
      localError="The specified incremented pressure index of "//TRIM(NumberToVString(pressureIncIdx,"*",err,error))// &
        & " is invalid. The incremented pressure index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(pressureIncBoundaryConditions%pressureIncrementedDOFIndices,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif    

    pressureIncDOFIndex=pressureIncBoundaryConditions%pressureIncrementedDOFIndices(pressureIncIdx)
      
    EXITS("BoundaryConditionsPressureInc_PressureIncDOFIndexGet")
    RETURN
999 ERRORS("BoundaryConditionsPressureInc_PressureIncDOFIndexGet",err,error)
    EXITS("BoundaryConditionsPressureInc_PressureIncDOFIndexGet")
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsPressureInc_PressureIncDOFIndexGet

  !
  !================================================================================================================================
  !

  !>Gets the boundary conditions for a boundary condition row variable.
  SUBROUTINE BoundaryConditionsRowVariable_BoundaryConditionsGet(boundaryConditionsRowVariable,boundaryConditions,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsRowVariableType), POINTER :: boundaryConditionsRowVariable !<A pointer to the boundary conditions row variable to get the boundary conditions for
    TYPE(BoundaryConditionsType), POINTER :: boundaryCOnditions !<On exit, a pointer to the boundary conditions in the specified boundary conditions row variable. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("BoundaryConditionsRowVariable_BoundaryConditionsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(boundaryConditions)) &
      & CALL FlagError("Boundary conditions is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(boundaryConditionsRowVariable)) &
      & CALL FlagError("Boundary conditions row variable is not associated.",err,error,*999)
#endif    

    boundaryConditions=>boundaryConditionsRowVariable%boundaryConditions

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(boundaryConditions)) &
      & CALL FlagError("Boundary conditions is not associated for the boundary conditions row variable.",err,error,*999)
#endif    
       
    EXITS("BoundaryConditionsRowVariable_BoundaryConditionsGet")
    RETURN
999 NULLIFY(boundaryConditions)
998 ERRORS("BoundaryConditionsRowVariable_BoundaryConditionsGet",err,error)
    EXITS("BoundaryConditionsRowVariable_BoundaryConditionsGet")
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsRowVariable_BoundaryConditionsGet

  !
  !================================================================================================================================
  !

  !>Checks the LHS field variable associated with a boundary condition row variable.
  SUBROUTINE BoundaryConditionsRowVariable_LHSVariableExists(boundaryConditionsRowVariable,lhsVariable,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsRowVariableType), POINTER :: boundaryConditionsRowVariable !<A pointer to the boundary conditions row variable to check the LHS field variable for
    TYPE(FieldVariableType), POINTER :: lhsVariable !<On exit, a pointer to the LHS field variable associated with the specified boundary conditions row variable if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("BoundaryConditionsRowVariable_LHSVariableExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(lhsVariable)) CALL FlagError("LHS field variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(boundaryConditionsRowVariable)) &
      & CALL FlagError("Boundary conditions row variable is not associated.",err,error,*999)
#endif    

    lhsVariable=>boundaryConditionsRowVariable%lhsVariable
      
    EXITS("BoundaryConditionsRowVariable_LHSVariableExists")
    RETURN
999 NULLIFY(lhsVariable)
998 ERRORS("BoundaryConditionsRowVariable_LHSVariableExists",err,error)
    EXITS("BoundaryConditionsRowVariable_LHSVariableExists")
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsRowVariable_LHSVariableExists

  !
  !================================================================================================================================
  !

  !>Gets the LHS field variable associated with a boundary condition row variable.
  SUBROUTINE BoundaryConditionsRowVariable_LHSVariableGet(boundaryConditionsRowVariable,lhsVariable,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsRowVariableType), POINTER :: boundaryConditionsRowVariable !<A pointer to the boundary conditions row variable to get the LHS field variable for
    TYPE(FieldVariableType), POINTER :: lhsVariable !<On exit, a pointer to the LHS field variable associated with the specified boundary conditions row variable. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("BoundaryConditionsRowVariable_LHSVariableGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(lhsVariable)) CALL FlagError("LHS field variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(boundaryConditionsRowVariable)) &
      & CALL FlagError("Boundary conditions row variable is not associated.",err,error,*999)
#endif    

    lhsVariable=>boundaryConditionsRowVariable%lhsVariable

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(lhsVariable)) &
      & CALL FlagError("The LHS field variable is not associated for the boundary conditions row variable.",err,error,*999)
#endif    
       
    EXITS("BoundaryConditionsRowVariable_LHSVariableGet")
    RETURN
999 NULLIFY(lhsVariable)
998 ERRORS("BoundaryConditionsRowVariable_LHSVariableGet",err,error)
    EXITS("BoundaryConditionsRowVariable_LHSVariableGet")
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsRowVariable_LHSVariableGet

  !
  !================================================================================================================================
  !

  !>Gets the number of DOFs for a boundary condition row variable.
  SUBROUTINE BoundaryConditionsRowVariable_NumberOfDOFsGet(boundaryConditionsRowVariable,numberOfDOFs,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsRowVariableType), POINTER :: boundaryConditionsRowVariable !<A pointer to the boundary conditions row variable to get the number of DOFs for
    INTEGER(INTG), INTENT(OUT) :: numberOfDOFs !<On exit, the number of DOFs for the specified boundary conditions row variable.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("BoundaryConditionsRowVariable_NumberOfDOFsGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(boundaryConditionsRowVariable)) &
      & CALL FlagError("Boundary conditions row variable is not associated.",err,error,*999)
#endif    

    numberOfDOFs=boundaryConditionsRowVariable%numberOfDOFs
      
    EXITS("BoundaryConditionsRowVariable_NumberOfDOFsGet")
    RETURN
999 ERRORS("BoundaryConditionsRowVariable_NumberOfDOFsGet",err,error)
    EXITS("BoundaryConditionsRowVariable_NumberOfDOFsGet")
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsRowVariable_NumberOfDOFsGet

  !
  !================================================================================================================================
  !

  !>Gets the RHS field variable associated with a boundary condition row variable if it exists.
  SUBROUTINE BoundaryConditionsRowVariable_RHSVariableExists(boundaryConditionsRowVariable,rhsVariable,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsRowVariableType), POINTER :: boundaryConditionsRowVariable !<A pointer to the boundary conditions row variable to get the RHS field variable for
    TYPE(FieldVariableType), POINTER :: rhsVariable !<On exit, a pointer to the rHS field variable associated with the specified boundary conditions row variable if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("BoundaryConditionsRowVariable_RHSVariableExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(rhsVariable)) CALL FlagError("RHS field variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(boundaryConditionsRowVariable)) &
      & CALL FlagError("Boundary conditions row variable is not associated.",err,error,*999)
#endif    

    rhsVariable=>boundaryConditionsRowVariable%rhsVariable
       
    EXITS("BoundaryConditionsRowVariable_RHSVariableExists")
    RETURN
999 NULLIFY(rhsVariable)
998 ERRORS("BoundaryConditionsRowVariable_RHSVariableExists",err,error)
    EXITS("BoundaryConditionsRowVariable_RHSVariableExists")
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsRowVariable_RHSVariableExists

  !
  !================================================================================================================================
  !

  !>Gets the row boundary condition type for a DOF in a boundary condition row variable.
  SUBROUTINE BoundaryConditionsRowVariable_RowConditionTypeGet(boundaryConditionsRowVariable,dofIdx,rowConditionType,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsRowVariableType), POINTER :: boundaryConditionsRowVariable !<A pointer to the boundary conditions row variable to get the row condition type for
    INTEGER(INTG), INTENT(IN) :: dofIdx !<The DOF index in the boundary conditions row variable to get the row condition type for.
    INTEGER(INTG), INTENT(OUT) :: rowConditionType !<On exit, the row condition type for the specified DOF in the boundary conditions row variable. \see BoundaryConditionsRoutines_RowTypes
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("BoundaryConditionsRowVariable_RowConditionTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(boundaryConditionsRowVariable)) &
      & CALL FlagError("Boundary conditions row variable is not associated.",err,error,*999)
    !TEMP WHILST WE SWITCH FROM GLOBAL DOFS
    !IF(dofIdx<1.OR.dofIdx>boundaryConditionsRowVariable%totalNumberOfDOFs) THEN
    !  localError="The specified DOF index of "//TRIM(NumberToVString(dofIdx,"*",err,error))// &
    !    & " is invalid. The specified DOF index should be >= 1 and <= "// &
    !    & TRIM(NumberToVString(boundaryConditionsRowVariable%totalNumberOfDOFs,"*",err,error))//"."
    !  CALL FlagError(localError,err,error,*999)
    !ENDIF
    IF(dofIdx<1.OR.dofIdx>boundaryConditionsRowVariable%numberOfGlobalDOFs) THEN
      localError="The specified DOF index of "//TRIM(NumberToVString(dofIdx,"*",err,error))// &
        & " is invalid. The specified DOF index should be >= 1 and <= "// &
        & TRIM(NumberToVString(boundaryConditionsRowVariable%numberOfGlobalDOFs,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(boundaryConditionsRowVariable%rowConditionTypes)) &
      & CALL FlagError("The row condition types array is not allocated for the boundary conditions row variable.",err,error,*999)
#endif    

    rowConditionType=boundaryConditionsRowVariable%rowConditionTypes(dofIdx)
      
    EXITS("BoundaryConditionsRowVariable_RowConditionTypeGet")
    RETURN
999 ERRORS("BoundaryConditionsRowVariable_RowConditionTypeGet",err,error)
    EXITS("BoundaryConditionsRowVariable_RowConditionTypeGet")
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsRowVariable_RowConditionTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the total number of DOFs for a boundary condition row variable.
  SUBROUTINE BoundaryConditionsRowVariable_TotalNumberOfDOFsGet(boundaryConditionsRowVariable,totalNumberOfDOFs,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsRowVariableType), POINTER :: boundaryConditionsRowVariable !<A pointer to the boundary conditions row variable to get the total number of DOFs for
    INTEGER(INTG), INTENT(OUT) :: totalNumberOfDOFs !<On exit, the total number of DOFs for the specified boundary conditions row variable.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("BoundaryConditionsRowVariable_TotalNumberOfDOFsGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(boundaryConditionsRowVariable)) &
      & CALL FlagError("Boundary conditions row variable is not associated.",err,error,*999)
#endif    

    totalNumberOfDOFs=boundaryConditionsRowVariable%totalNumberOfDOFs
      
    EXITS("BoundaryConditionsRowVariable_TotalNumberOfDOFsGet")
    RETURN
999 ERRORS("BoundaryConditionsRowVariable_TotalNumberOfDOFsGet",err,error)
    EXITS("BoundaryConditionsRowVariable_TotalNumberOfDOFsGet")
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsRowVariable_TotalNumberOfDOFsGet

  !
  !================================================================================================================================
  !

  !>Gets the boundary conditions for a boundary condition variable.
  SUBROUTINE BoundaryConditionsVariable_BoundaryConditionsGet(boundaryConditionsVariable,boundaryConditions,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to get the boundary conditions for
    TYPE(BoundaryConditionsType), POINTER :: boundaryCOnditions !<On exit, a pointer to the boundary conditions in the specified boundary conditions variable. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("BoundaryConditionsVariable_BoundaryConditionsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(boundaryConditions)) &
      & CALL FlagError("Boundary conditions is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(boundaryConditionsVariable)) &
      & CALL FlagError("Boundary conditions variable is not associated.",err,error,*999)
#endif    

    boundaryConditions=>boundaryConditionsVariable%boundaryConditions

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(boundaryConditions)) &
      & CALL FlagError("Boundary conditions is not associated for the boundary conditions variable.",err,error,*999)
#endif    
       
    EXITS("BoundaryConditionsVariable_BoundaryConditionsGet")
    RETURN
999 NULLIFY(boundaryConditions)
998 ERRORS("BoundaryConditionsVariable_BoundaryConditionsGet",err,error)
    EXITS("BoundaryConditionsVariable_BoundaryConditionsGet")
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsVariable_BoundaryConditionsGet

  !
  !================================================================================================================================
  !

  !>Gets the condition type for a DOF in a boundary condition variable.
  SUBROUTINE BoundaryConditionsVariable_ConditionTypeGet(boundaryConditionsVariable,dofIdx,conditionType,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to get the condition type for
    INTEGER(INTG), INTENT(IN) :: dofIdx !<The dof index to get the condition type for.
    INTEGER(INTG), INTENT(OUT) :: conditionType !<On exit, the boundary condition variable condition type.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    INTEGER(INTG) :: numberOfGlobalDOFs
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("BoundaryConditionsVariable_ConditionTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(boundaryConditionsVariable)) &
      & CALL FlagError("Boundary conditions variable is not associated.",err,error,*999)
    NULLIFY(fieldVariable)
    CALL BoundaryConditionsVariable_FieldVariableGet(boundaryConditionsVariable,fieldVariable,err,error,*999)
    CALL FieldVariable_NumberOfGlobalDOFsGet(fieldVariable,numberOfGlobalDOFs,err,error,*999)
    IF(dofIdx<1.OR.dofIdx>numberOfGlobalDOFs) THEN
      localError="The specified DOF index of "//TRIM(NumberToVString(dofIdx,"*",err,error))// &
        & " is invalid. The DOF index should be >= 1 and <= "//TRIM(NumberToVString(numberOfGlobalDOFs,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(boundaryConditionsVariable%conditionTypes)) &
      & CALL FlagError("The condition types array is not allocated for the boundary conditions variable.",err,error,*999)
#endif    

    conditionType=boundaryConditionsVariable%conditionTypes(dofIdx)
    
    EXITS("BoundaryConditionsVariable_ConditionTypeGet")
    RETURN
999 ERRORS("BoundaryConditionsVariable_ConditionTypeGet",err,error)
    EXITS("BoundaryConditionsVariable_ConditionTypeGet")
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsVariable_ConditionTypeGet

  !
  !================================================================================================================================
  !

  !>Checks if the dirichlet boundary conditions exist for a boundary condition variable.
  SUBROUTINE BoundaryConditionsVariable_DirichletConditionsExists(boundaryConditionsVariable,dirichletBoundaryConditions, &
    & err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to check the existance of the dirichlet boundary conditions for
    TYPE(BoundaryConditionsDirichletType), POINTER :: dirichletBoundaryCOnditions !<On exit, a pointer to the dirichlet boundary conditions in the specified boundary conditions variable if they exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("BoundaryConditionsVariable_DirichletConditionsExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dirichletBoundaryConditions)) &
      & CALL FlagError("Dirichlet boundary conditions is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(boundaryConditionsVariable)) &
      & CALL FlagError("Boundary conditions variable is not associated.",err,error,*999)
#endif    

    dirichletBoundaryConditions=>boundaryConditionsVariable%dirichletBoundaryConditions

    EXITS("BoundaryConditionsVariable_DirichletConditionsExists")
    RETURN
999 NULLIFY(dirichletBoundaryConditions)
998 ERRORS("BoundaryConditionsVariable_DirichletConditionsExists",err,error)
    EXITS("BoundaryConditionsVariable_DirichletConditionsExists")
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsVariable_DirichletConditionsExists

  !
  !================================================================================================================================
  !

  !>Gets the dirichlet boundary conditions for a boundary condition variable.
  SUBROUTINE BoundaryConditionsVariable_DirichletConditionsGet(boundaryConditionsVariable,dirichletBoundaryConditions, &
    & err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to get the dirichlet boundary conditions for
    TYPE(BoundaryConditionsDirichletType), POINTER :: dirichletBoundaryCOnditions !<On exit, a pointer to the dirichlet boundary conditions in the specified boundary conditions variable. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("BoundaryConditionsVariable_DirichletConditionsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dirichletBoundaryConditions)) &
      & CALL FlagError("Dirichlet boundary conditions is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(boundaryConditionsVariable)) &
      & CALL FlagError("Boundary conditions variable is not associated.",err,error,*999)
#endif    

    dirichletBoundaryConditions=>boundaryConditionsVariable%dirichletBoundaryConditions

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(dirichletBoundaryConditions)) &
      & CALL FlagError("Dirichlet boundary conditions is not associated for the boundary conditions variable.",err,error,*999)
#endif    
       
    EXITS("BoundaryConditionsVariable_DirichletConditionsGet")
    RETURN
999 NULLIFY(dirichletBoundaryConditions)
998 ERRORS("BoundaryConditionsVariable_DirichletConditionsGet",err,error)
    EXITS("BoundaryConditionsVariable_DirichletConditionsGet")
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsVariable_DirichletConditionsGet

  !
  !================================================================================================================================
  !

  !>Gets the DOF constraints for a boundary condition variable if they exist
  SUBROUTINE BoundaryConditionsVariable_DOFConstraintsExists(boundaryConditionsVariable,dofConstraints,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to get the dirichlet boundary conditions for
    TYPE(BoundaryConditionsDofConstraintsType), POINTER :: dofConstraints !<On exit, a pointer to the DOF constraints in the specified boundary conditions variable if they exist. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("BoundaryConditionsVariable_DOFConstraintsExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dofConstraints)) CALL FlagError("DOF constraints is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(boundaryConditionsVariable)) CALL FlagError("Boundary conditions variable is not associated.",err,error,*999)
#endif    

    dofConstraints=>boundaryConditionsVariable%dofConstraints
       
    EXITS("BoundaryConditionsVariable_DOFConstraintsExists")
    RETURN
999 NULLIFY(dofConstraints)
998 ERRORS("BoundaryConditionsVariable_DOFConstraintsExists",err,error)
    EXITS("BoundaryConditionsVariable_DOFConstraintsExists")
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsVariable_DOFConstraintsExists

  !
  !================================================================================================================================
  !

  !>Gets the DOF constraints for a boundary condition variable
  SUBROUTINE BoundaryConditionsVariable_DOFConstraintsGet(boundaryConditionsVariable,dofConstraints,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to get the dirichlet boundary conditions for
    TYPE(BoundaryConditionsDofConstraintsType), POINTER :: dofConstraints !<On exit, a pointer to the DOF constraints in the specified boundary conditions variable. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("BoundaryConditionsVariable_DOFConstraintsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dofConstraints)) CALL FlagError("DOF constraints is already associated.",err,error,*998)
#endif    
    
    CALL BoundaryConditionsVariable_DOFConstraintsExists(boundaryConditionsVariable,dofConstraints,err,error,*999)

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(dofConstraints)) &
      & CALL FlagError("DOF constraints is not associated for the boundary conditions variable.",err,error,*999)
#endif    

    EXITS("BoundaryConditionsVariable_DOFConstraintsGet")
    RETURN
999 NULLIFY(dofConstraints)
998 ERRORS("BoundaryConditionsVariable_DOFConstraintsGet",err,error)
    EXITS("BoundaryConditionsVariable_DOFConstraintsGet")
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsVariable_DOFConstraintsGet

  !
  !================================================================================================================================
  !

  !>Gets the DOF count for a boundary condition type in a boundary condition variable.
  SUBROUTINE BoundaryConditionsVariable_DOFCountGet(boundaryConditionsVariable,conditionType,dofCount,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to get the DOF count for
    INTEGER(INTG), INTENT(IN) :: conditionType !<The boundary condition type to get the DOF count for.
    INTEGER(INTG), INTENT(OUT) :: dofCount !<On exit, the count of DOFs with the specified boundary condition type.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
   TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("BoundaryConditionsVariable_DOFCountGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(boundaryConditionsVariable)) CALL FlagError("Boundary conditions variable is not associated.",err,error,*999)
    IF(conditionType<1.OR.conditionType>MAX_BOUNDARY_CONDITION_NUMBER) THEN
      localError="The specified boundary condition type of "//TRIM(NumberToVString(conditionType,"*",err,error))// &
        & " is invalid. The boundary condition type should be >= 1 and <= "// &
        & TRIM(NumberToVString(MAX_BOUNDARY_CONDITION_NUMBER,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(boundaryConditionsVariable%dofCounts)) &
      & CALL FlagError("The DOF counts array is not allocated for the boundary conditions variable.",err,error,*999)
#endif    

    dofCount=boundaryConditionsVariable%dofCounts(conditionType)
    
    EXITS("BoundaryConditionsVariable_DOFCountGet")
    RETURN
999 ERRORS("BoundaryConditionsVariable_DOFCountGet",err,error)
    EXITS("BoundaryConditionsVariable_DOFCountGet")
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsVariable_DOFCountGet

  !
  !================================================================================================================================
  !

  !>Gets the DOF type for a DOF in a boundary condition variable.
  SUBROUTINE BoundaryConditionsVariable_DOFTypeGet(boundaryConditionsVariable,dofIdx,dofType,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to get the DOF type for
    INTEGER(INTG), INTENT(IN) :: dofIdx !<The DOF index to get the DOF type for.
    INTEGER(INTG), INTENT(OUT) :: dofType !<On exit, the boundary condition variable DOF type.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS
    INTEGER(INTG) :: numberOfGlobalDOFs
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("BoundaryConditionsVariable_DOFTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(boundaryConditionsVariable)) &
      & CALL FlagError("Boundary conditions variable is not associated.",err,error,*999)
    NULLIFY(fieldVariable)
    CALL BoundaryConditionsVariable_FieldVariableGet(boundaryConditionsVariable,fieldVariable,err,error,*999)
    CALL FieldVariable_NumberOfGlobalDOFsGet(fieldVariable,numberOfGlobalDOFs,err,error,*999)
    IF(dofIdx<1.OR.dofIdx>numberOfGlobalDOFs) THEN
      localError="The specified DOF index of "//TRIM(NumberToVString(dofIdx,"*",err,error))// &
        & " is invalid. The DOF index should be >= 1 and <= "//TRIM(NumberToVString(numberOfGlobalDOFs,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(boundaryConditionsVariable%dofTypes)) &
      & CALL FlagError("The dof types array is not allocated for the boundary conditions variable.",err,error,*999)
#endif    

    dofType=boundaryConditionsVariable%dofTypes(dofIdx)
    
    EXITS("BoundaryConditionsVariable_DOFTypeGet")
    RETURN
999 ERRORS("BoundaryConditionsVariable_DOFTypeGet",err,error)
    EXITS("BoundaryConditionsVariable_DOFTypeGet")
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsVariable_DOFTypeGet

  !
  !================================================================================================================================
  !

  !>Gets the field variable associated with a boundary condition variable.
  SUBROUTINE BoundaryConditionsVariable_FieldVariableGet(boundaryConditionsVariable,fieldVariable,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to get the field variable for
    TYPE(FieldVariableType), POINTER :: fieldVariable !<On exit, a pointer to the field variable associated with the specified boundary conditions variable. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("BoundaryConditionsVariable_FieldVariableGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(boundaryConditionsVariable)) CALL FlagError("Boundary conditions variable is not associated.",err,error,*999)
#endif    

    fieldVariable=>boundaryConditionsVariable%variable

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(fieldVariable)) &
      & CALL FlagError("The field variable is not associated for the boundary conditions variable.",err,error,*999)
#endif    
       
    EXITS("BoundaryConditionsVariable_FieldVariableGet")
    RETURN
999 NULLIFY(fieldVariable)
998 ERRORS("BoundaryConditionsVariable_FieldVariableGet",err,error)
    EXITS("BoundaryConditionsVariable_FieldVariableGet")
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsVariable_FieldVariableGet

  !
  !================================================================================================================================
  !

  !>Gets the field variable type associated with a boundary condition variable.
  SUBROUTINE BoundaryConditionsVariable_FieldVariableTypeGet(boundaryConditionsVariable,fieldVariableType,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to get the field variable type for
    INTEGER(INTG), INTENT(OUT) :: fieldVariableType !<On exit, the field variable type associated with the specified boundary conditions variable.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("BoundaryConditionsVariable_FieldVariableTypeGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(boundaryConditionsVariable)) CALL FlagError("Boundary conditions variable is not associated.",err,error,*999)
#endif    

    fieldVariableType=boundaryConditionsVariable%variableType
      
    EXITS("BoundaryConditionsVariable_FieldVariableTypeGet")
    RETURN
999 ERRORS("BoundaryConditionsVariable_FieldVariableTypeGet",err,error)
    EXITS("BoundaryConditionsVariable_FieldVariableTypeGet")
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsVariable_FieldVariableTypeGet

  !
  !================================================================================================================================
  !

  !>Checks if the neumann boundary conditions for a boundary condition variable exists.
  SUBROUTINE BoundaryConditionsVariable_NeumannConditionsExists(boundaryConditionsVariable,neumannBoundaryConditions,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to check the existance of the neumann boundary conditions for
    TYPE(BoundaryConditionsNeumannType), POINTER :: neumannBoundaryCOnditions !<On exit, a pointer to the neumann boundary conditions in the specified boundary conditions variable if they exist. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("BoundaryConditionsVariable_NeumannConditionsExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(neumannBoundaryConditions)) CALL FlagError("Neumann boundary conditions is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(boundaryConditionsVariable)) CALL FlagError("Boundary conditions variable is not associated.",err,error,*999)
#endif    

    neumannBoundaryConditions=>boundaryConditionsVariable%neumannBoundaryConditions

    EXITS("BoundaryConditionsVariable_NeumannConditionsExists")
    RETURN
999 NULLIFY(neumannBoundaryConditions)
998 ERRORS("BoundaryConditionsVariable_NeumannConditionsExists",err,error)
    EXITS("BoundaryConditionsVariable_NeumannConditionsExists")
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsVariable_NeumannConditionsExists

  !
  !================================================================================================================================
  !

  !>Gets the neumann boundary conditions for a boundary condition variable.
  SUBROUTINE BoundaryConditionsVariable_NeumannConditionsGet(boundaryConditionsVariable,neumannBoundaryConditions,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to get the neumann boundary conditions for
    TYPE(BoundaryConditionsNeumannType), POINTER :: neumannBoundaryCOnditions !<On exit, a pointer to the neumann boundary conditions in the specified boundary conditions variable. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("BoundaryConditionsVariable_NeumannConditionsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(neumannBoundaryConditions)) CALL FlagError("Neumann boundary conditions is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(boundaryConditionsVariable)) CALL FlagError("Boundary conditions variable is not associated.",err,error,*999)
#endif    

    neumannBoundaryConditions=>boundaryConditionsVariable%neumannBoundaryConditions

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(neumannBoundaryConditions)) &
      & CALL FlagError("Neumann boundary conditions is not associated for the boundary conditions variable.",err,error,*999)
#endif    
       
    EXITS("BoundaryConditionsVariable_NeumannConditionsGet")
    RETURN
999 NULLIFY(neumannBoundaryConditions)
998 ERRORS("BoundaryConditionsVariable_NeumannConditionsGet",err,error)
    EXITS("BoundaryConditionsVariable_NeumannConditionsGet")
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsVariable_NeumannConditionsGet

  !
  !================================================================================================================================
  !

  !>Gets the number of Dirichlet conditions for a boundary condition variable.
  SUBROUTINE BoundaryConditionsVariable_NumberOfDirichletConditionsGet(boundaryConditionsVariable,numberOfDirichletConditions, &
    & err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to get the number of Dirichlet conditions for.
    INTEGER(INTG), INTENT(OUT) :: numberOfDirichletConditions !<On exit, the number of Dirichlet conditions.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("BoundaryConditionsVariable_NumberOfDirichletConditionsGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(boundaryConditionsVariable)) CALL FlagError("Boundary conditions variable is not associated.",err,error,*999)
#endif    

    numberOfDirichletConditions=boundaryConditionsVariable%numberOfDirichletConditions
      
    EXITS("BoundaryConditionsVariable_NumberOfDirichletConditionsGet")
    RETURN
999 ERRORS("BoundaryConditionsVariable_NumberOfDirichletConditionsGet",err,error)
    EXITS("BoundaryConditionsVariable_NumberOfDirichletConditionsGet")
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsVariable_NumberOfDirichletConditionsGet

  !
  !================================================================================================================================
  !

  !>Checks if the pressure incremented boundary conditions for a boundary condition variable exist.
  SUBROUTINE BoundaryConditionsVariable_PressureIncConditionsExists(boundaryConditionsVariable, &
    & pressureIncrementedBoundaryConditions,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to check the existance of the pressure incremented boundary conditions for
    TYPE(BoundaryConditionsPressureIncrementedType), POINTER :: pressureIncrementedBoundaryConditions !<On exit, a pointer to the pressure incremented boundary conditions in the specified boundary conditions variable if the exist. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("BoundaryConditionsVariable_PressureIncConditionsExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(pressureIncrementedBoundaryConditions)) &
      & CALL FlagError("Pressure incremented boundary conditions is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(boundaryConditionsVariable)) CALL FlagError("Boundary conditions variable is not associated.",err,error,*999)
#endif    

    pressureIncrementedBoundaryConditions=>boundaryConditionsVariable%pressureIncrementedBoundaryConditions
      
    EXITS("BoundaryConditionsVariable_PressureIncConditionsExists")
    RETURN
999 NULLIFY(pressureIncrementedBoundaryConditions)
998 ERRORS("BoundaryConditionsVariable_PressureIncConditionsExists",err,error)
    EXITS("BoundaryConditionsVariable_PressureIncConditionsExists")
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsVariable_PressureIncConditionsExists

  !
  !================================================================================================================================
  !

  !>Gets the pressure incremented boundary conditions for a boundary condition variable.
  SUBROUTINE BoundaryConditionsVariable_PressureIncConditionsGet(boundaryConditionsVariable, &
    & pressureIncrementedBoundaryConditions,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to get the pressure incremented boundary conditions for
    TYPE(BoundaryConditionsPressureIncrementedType), POINTER :: pressureIncrementedBoundaryConditions !<On exit, a pointer to the pressure incremented boundary conditions in the specified boundary conditions variable. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("BoundaryConditionsVariable_PressureIncConditionsGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(pressureIncrementedBoundaryConditions)) &
      & CALL FlagError("Pressure incremented boundary conditions is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(boundaryConditionsVariable)) CALL FlagError("Boundary conditions variable is not associated.",err,error,*999)
#endif    

    pressureIncrementedBoundaryConditions=>boundaryConditionsVariable%pressureIncrementedBoundaryConditions

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(pressureIncrementedBoundaryConditions)) &
      & CALL FlagError("Pressure incremented boundary conditions is not associated for the boundary conditions variable.", &
      & err,error,*999)
#endif    
       
    EXITS("BoundaryConditionsVariable_PressureIncConditionsGet")
    RETURN
999 NULLIFY(pressureIncrementedBoundaryConditions)
998 ERRORS("BoundaryConditionsVariable_PressureIncConditionsGet",err,error)
    EXITS("BoundaryConditionsVariable_PressureIncConditionsGet")
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsVariable_PressureIncConditionsGet

  !
  !================================================================================================================================
  !

  !>Gets the DOF constraints for a boundary condition variable DOF constraints if it exists
  SUBROUTINE BoundaryConditionsVariableDOFConstraints_ConstraintExists(dofConstraints,constraintIdx,dofConstraint,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsDofConstraintsType), POINTER :: dofConstraints !<A pointer to the DOF constraints to get the DOF constraint for
    INTEGER(INTG), INTENT(IN) :: constraintIdx !<The constraint index to get the DOF constraint for
    TYPE(BoundaryConditionsDofConstraintType), POINTER :: dofConstraint !<On exit, a pointer to the specified DOF constraint if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("BoundaryConditionsVariableDOFConstraints_ConstraintExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dofConstraint)) CALL FlagError("DOF constraint is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(dofConstraints)) CALL FlagError("DOF constraints is not associated.",err,error,*999)
    IF(constraintIdx<1.OR.constraintIdx>dofConstraints%numberOfConstraints) THEN
      localError="The specified constraint index of "//TRIM(NumberToVString(constraintIdx,"*",err,error))// &
        & " is invalid. The constraint index should be >=1 and <= "// &
        & TRIM(NumberToVString(dofConstraints%numberOfConstraints,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(dofConstraints%constraints)) CALL FlagError("DOF constraints constraints is not allocated.",err,error,*999)
#endif    

    dofConstraint=>dofConstraints%constraints(constraintIdx)%ptr
   
    EXITS("BoundaryConditionsVariableDOFConstraints_ConstraintExists")
    RETURN
999 NULLIFY(dofConstraint)
998 ERRORS("BoundaryConditionsVariableDOFConstraints_ConstraintExists",err,error)
    EXITS("BoundaryConditionsVariableDOFConstraints_ConstraintExists")
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsVariableDOFConstraints_ConstraintExists

  !
  !================================================================================================================================
  !

  !>Gets the DOF constraints for a boundary condition variable DOF constraints
  SUBROUTINE BoundaryConditionsVariableDOFConstraints_ConstraintGet(dofConstraints,constraintIdx,dofConstraint,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsDofConstraintsType), POINTER :: dofConstraints !<A pointer to the DOF constraints to get the DOF constraint for
    INTEGER(INTG), INTENT(IN) :: constraintIdx !<The constraint index to get the DOF constraint for
    TYPE(BoundaryConditionsDofConstraintType), POINTER :: dofConstraint !<On exit, a pointer to the specified DOF constraint. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("BoundaryConditionsVariableDOFConstraints_ConstraintGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dofConstraint)) CALL FlagError("DOF onstraint is already associated.",err,error,*998)
#endif    

    CALL BoundaryConditionsVariableDOFConstraints_ConstraintExists(dofConstraints,constraintIdx,dofConstraint,err,error,*999)

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(dofConstraint)) THEN
      localError="The DOF constraint is not associated for constraint index "// &
        & TRIM(NumberToVString(constraintIdx,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ELSE
      IF(dofConstraint%numberOfDOFs>0) THEN
        IF(.NOT.ALLOCATED(dofConstraint%dofs)) THEN
          localError="The DOFs are not allocated for constraint index "// &
            & TRIM(NumberToVString(constraintIdx,"*",err,error))//" which has "// &
            & TRIM(NumberToVString(dofConstraint%numberOfDOFs,"*",err,error))//" DOFs."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        IF(.NOT.ALLOCATED(dofConstraint%coefficients)) THEN
          localError="The coefficients are not allocated for constraint index "// &
            & TRIM(NumberToVString(constraintIdx,"*",err,error))//" which has "// &
            & TRIM(NumberToVString(dofConstraint%numberOfDOFs,"*",err,error))//" DOFs."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDIF
    ENDIF
#endif    
   
    EXITS("BoundaryConditionsVariableDOFConstraints_ConstraintGet")
    RETURN
999 NULLIFY(dofConstraint)
998 ERRORS("BoundaryConditionsVariableDOFConstraints_ConstraintGet",err,error)
    EXITS("BoundaryConditionsVariableDOFConstraints_ConstraintGet")
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsVariableDOFConstraints_ConstraintGet

  !
  !================================================================================================================================
  !

  !>Gets the DOF coupling for a boundary condition variable DOF constraints if it exists
  SUBROUTINE BoundaryConditionsVariableDOFConstraints_DOFCouplingExists(dofConstraints,dofIdx,dofCoupling,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsDofConstraintsType), POINTER :: dofConstraints !<A pointer to the DOF constraints to get the DOF coupling for
    INTEGER(INTG), INTENT(IN) :: dofIdx !<The dof index to get the DOF coupling for
    TYPE(BoundaryConditionsCoupledDofsType), POINTER :: dofCoupling !<On exit, a pointer to the specified DOF coupling if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_PRECHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("BoundaryConditionsVariableDOFConstraints_DOFCouplingExists",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dofCoupling)) CALL FlagError("DOF coupling is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(dofConstraints)) CALL FlagError("DOF constraints is not associated.",err,error,*999)
    IF(dofIdx<1.OR.dofIdx>dofConstraints%numberOfDOFs) THEN
      localError="The specified DOF index of "//TRIM(NumberToVString(dofIdx,"*",err,error))// &
        & " is invalid. The DOF index should be >=1 and <= "// &
        & TRIM(NumberToVString(dofConstraints%numberOfDOFs,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(dofConstraints%dofCouplings)) CALL FlagError("DOF constraints DOF couplings is not allocated.",err,error,*999)
#endif    

    dofCoupling=>dofConstraints%dofCouplings(dofIdx)%ptr
    
    EXITS("BoundaryConditionsVariableDOFConstraints_DOFCouplingExists")
    RETURN
999 NULLIFY(dofCoupling)
998 ERRORS("BoundaryConditionsVariableDOFConstraints_DOFCouplingExists",err,error)
    EXITS("BoundaryConditionsVariableDOFConstraints_DOFCouplingExists")
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsVariableDOFConstraints_DOFCouplingExists

  !
  !================================================================================================================================
  !

  !>Gets the DOF coupling for a boundary condition variable DOF constraints
  SUBROUTINE BoundaryConditionsVariableDOFConstraints_DOFCouplingGet(dofConstraints,dofIdx,dofCoupling,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsDofConstraintsType), POINTER :: dofConstraints !<A pointer to the DOF constraints to get the DOF coupling for
    INTEGER(INTG), INTENT(IN) :: dofIdx !<The dof index to get the DOF coupling for
    TYPE(BoundaryConditionsCoupledDofsType), POINTER :: dofCoupling !<On exit, a pointer to the specified DOF coupling. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
#ifdef WITH_POSTCHECKS    
    TYPE(VARYING_STRING) :: localError
#endif    
 
    ENTERS("BoundaryConditionsVariableDOFConstraints_DOFCouplingGet",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(dofCoupling)) CALL FlagError("DOF coupling is already associated.",err,error,*998)
#endif    

    CALL BoundaryConditionsVariableDOFConstraints_DOFCouplingExists(dofConstraints,dofIdx,dofCoupling,err,error,*999)

#ifdef WITH_POSTCHECKS    
    IF(.NOT.ASSOCIATED(dofCoupling)) THEN
      localError="The DOF couping is not associated for DOF index "// &
        & TRIM(NumberToVString(dofIdx,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ELSE
      IF(dofCoupling%numberOfDOFs>0) THEN
        IF(.NOT.ALLOCATED(dofCoupling%globalDofs)) THEN
          localError="The global DOFs are not allocated for DOF index "// &
            & TRIM(NumberToVString(dofIdx,"*",err,error))//" which has "// &
            & TRIM(NumberToVString(dofConstraints%numberOfDOFs,"*",err,error))//" DOFs."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        IF(.NOT.ALLOCATED(dofCoupling%localDofs)) THEN
          localError="The local DOFs are not allocated for DOF index "// &
            & TRIM(NumberToVString(dofIdx,"*",err,error))//" which has "// &
            & TRIM(NumberToVString(dofConstraints%numberOfDOFs,"*",err,error))//" DOFs."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        IF(.NOT.ALLOCATED(dofCoupling%coefficients)) THEN
          localError="The coefficients are not allocated for DOF index "// &
            & TRIM(NumberToVString(dofIdx,"*",err,error))//" which has "// &
            & TRIM(NumberToVString(dofConstraints%numberOfDOFs,"*",err,error))//" DOFs."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDIF
    ENDIF
#endif    
    
    EXITS("BoundaryConditionsVariableDOFConstraints_DOFCouplingGet")
    RETURN
999 NULLIFY(dofCoupling)
998 ERRORS("BoundaryConditionsVariableDOFConstraints_DOFCouplingGet",err,error)
    EXITS("BoundaryConditionsVariableDOFConstraints_DOFCouplingGet")
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsVariableDOFConstraints_DOFCouplingGet

  !
  !================================================================================================================================
  !

  !>Gets the number of DOF constraints for a boundary condition variable DOF constraints
  SUBROUTINE BoundaryConditionsVariableDOFConstraints_NumberOfConstraintsGet(dofConstraints,numberOfConstraints,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsDofConstraintsType), POINTER :: dofConstraints !<A pointer to the DOF constraints to get the number of DOF constraints for
    INTEGER(INTG), INTENT(OUT) :: numberOfConstraints !<On exit, the number of constraints in the DOF constraints
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("BoundaryConditionsVariableDOFConstraints_NumberOfConstraintsGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(dofConstraints)) CALL FlagError("DOF constraints is not associated.",err,error,*999)
#endif    
 
    numberOfConstraints=dofConstraints%numberOfConstraints
    
    EXITS("BoundaryConditionsVariableDOFConstraints_NumberOfConstraintsGet")
    RETURN
999 ERRORS("BoundaryConditionsVariableDOFConstraints_NumberOfConstraintsGet",err,error)
    EXITS("BoundaryConditionsVariableDOFConstraints_NumberOfConstraintsGet")
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsVariableDOFConstraints_NumberOfConstraintsGet

  !
  !================================================================================================================================
  !

  !>Gets the number of DOFs for a boundary condition variable DOF constraints
  SUBROUTINE BoundaryConditionsVariableDOFConstraints_NumberOfDOFsGet(dofConstraints,numberOfDOFs,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsDofConstraintsType), POINTER :: dofConstraints !<A pointer to the DOF constraints to get the number of DOFs for
    INTEGER(INTG), INTENT(OUT) :: numberOfDOFs !<On exit, the number of DOFs in the DOF constraints
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("BoundaryConditionsVariableDOFConstraints_NumberOfDOFsGet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(dofConstraints)) CALL FlagError("DOF constraints is not associated.",err,error,*999)
#endif    
 
    numberOfDOFs=dofConstraints%numberOfDOFs
    
    EXITS("BoundaryConditionsVariableDOFConstraints_NumberOfDOFsGet")
    RETURN
999 ERRORS("BoundaryConditionsVariableDOFConstraints_NumberOfDOFsGet",err,error)
    EXITS("BoundaryConditionsVariableDOFConstraints_NumberOfDOFsGet")
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsVariableDOFConstraints_NumberOfDOFsGet

  !
  !================================================================================================================================
  !

END MODULE BoundaryConditionAccessRoutines
