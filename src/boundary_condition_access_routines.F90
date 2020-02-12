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

!> This module contains all boundary condition access method routines.
MODULE BoundaryConditionAccessRoutines
  
  USE BaseRoutines
  USE ISO_VARYING_STRING
  USE Kinds
  USE Strings
  USE Types
  
#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup BoundaryConditionsRoutines_DOFTypes BoundaryConditionsRoutines::DOFTypes
  !> \brief DOF type for boundary conditions.
  !>@{
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_DOF_FREE=0 !<The dof is free. \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_DOF_FIXED=1 !<The dof is fixed as a boundary condition. \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_DOF_MIXED=2 !<The dof is set as a mixed boundary condition. \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_DOF_CONSTRAINED=3 !<The dof is constrained to be a linear combination of other DOFs. \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
  !>@}
  !> \addtogroup BoundaryConditionsRoutines_BoundaryConditions BoundaryConditionsRoutines::BoundaryConditions
  !> \brief Boundary conditions types. These may be specific to a particular equation type and the solver routines should not need to use these.
  !>@{
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FREE=0 !<The dof is free. \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED=1 !<The dof is fixed as a boundary condition. \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED_INLET=2 !<The dof is fixed as a boundary condition. \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED_OUTLET=3 !<The dof is fixed as a boundary condition. \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED_WALL=4 !<The dof is fixed as a boundary condition. \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_MOVED_WALL=5 !<The dof is fixed as a boundary condition. \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FREE_WALL=6 !<The dof is fixed as a boundary condition. \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_NEUMANN_POINT=8 !<The dof is set to a Neumann point boundary condition. \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_NEUMANN_INTEGRATED=9 !<The dof is set to a Neumann integrated boundary condition. \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_DIRICHLET=10 !<The dof is set to a Dirichlet boundary condition. \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_CAUCHY=11 !<The dof is set to a Cauchy boundary condition. \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_ROBIN=12 !<The dof is set to a Robin boundary condition. \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED_INCREMENTED=13 !<The dof is a fixed boundary condition, to be used with load increment loop. \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_PRESSURE=14 !<The dof is a surface pressure boundary condition. \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_PRESSURE_INCREMENTED=15 !<The dof is a surface pressure boundary condition, to be used with load increment loop. \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED=17 !<The dof is fixed as a boundary condition, to be used with load increment loop. \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_CORRECTION_MASS_INCREASE=18 !<The dof is fixed as a boundary condition, to be used with load increment loop. \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_IMPERMEABLE_WALL=19 !<The dof is set such that (via penalty formulation): velocity * normal = 0. \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_NEUMANN_INTEGRATED_ONLY=20 !<A Neumann integrated boundary condition, and no point values will be integrated over a face or line that includes this dof. \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_LINEAR_CONSTRAINT=21 !<The dof is constrained to be a linear combination of other DOFs. \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED=22!<A Neumann point boundary condition that is incremented inside a load increment control loop. \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED_FITTED=23 !<The dof is fixed as a boundary condition to be updated from fitting data \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED_NONREFLECTING=24 !<The dof is fixed and set to a non-reflecting type for 1D wave propagation problems. \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED_CELLML=25 !<The dof is fixed and set to values specified based on the coupled CellML solution at the dof. \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED_STREE=26 !<The dof is fixed and set to values specified based on the transmission line theory at the dof. \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_COUPLING_FLOW=27 !<The dof is fixed and set to values specified based on a coupled flow rate at the dof. \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_COUPLING_STRESS=28 !<The dof is fixed and set to values specified based on a coupled stress at the dof. \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FIXED_PRESSURE=29 !<The dof is a fixed pressure boundary condition. \see BoundaryConditionsRoutines_BoundaryConditions,BoundaryConditionsRoutines
 !>@}

  INTEGER(INTG), PARAMETER :: MAX_BOUNDARY_CONDITION_NUMBER=29 !The maximum boundary condition type identifier, used for allocating an array with an entry for each type

  !> \addtogroup BoundaryConditionsRoutines_SparsityTypes BoundaryConditionsRoutines::BoundaryConditions
  !> \brief Storage type for matrices used by boundary conditions.
  !>@{
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_SPARSE_MATRICES=1 !<The matrices are stored as sparse matrices.
  INTEGER(INTG), PARAMETER :: BOUNDARY_CONDITION_FULL_MATRICES=2 !<The matrices are stored as full matrices.
  !>@}

  !Module types

  !Module variables

  !Interfaces

  PUBLIC BOUNDARY_CONDITION_DOF_FREE,BOUNDARY_CONDITION_DOF_FIXED,BOUNDARY_CONDITION_DOF_MIXED,BOUNDARY_CONDITION_DOF_CONSTRAINED

  PUBLIC BOUNDARY_CONDITION_FREE,BOUNDARY_CONDITION_FIXED,BOUNDARY_CONDITION_FIXED_INLET,&
    & BOUNDARY_CONDITION_FIXED_OUTLET,BOUNDARY_CONDITION_FIXED_WALL,BOUNDARY_CONDITION_MOVED_WALL,BOUNDARY_CONDITION_FREE_WALL,&
    & BOUNDARY_CONDITION_NEUMANN_INTEGRATED,BOUNDARY_CONDITION_DIRICHLET,BOUNDARY_CONDITION_NEUMANN_POINT, &
    & BOUNDARY_CONDITION_CAUCHY,BOUNDARY_CONDITION_ROBIN,BOUNDARY_CONDITION_FIXED_INCREMENTED,BOUNDARY_CONDITION_PRESSURE,&
    & BOUNDARY_CONDITION_PRESSURE_INCREMENTED,BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED, &
    & BOUNDARY_CONDITION_CORRECTION_MASS_INCREASE,BOUNDARY_CONDITION_IMPERMEABLE_WALL,BOUNDARY_CONDITION_NEUMANN_INTEGRATED_ONLY, &
    & BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED,BOUNDARY_CONDITION_FIXED_STREE, &
    & BOUNDARY_CONDITION_FIXED_FITTED,BOUNDARY_CONDITION_FIXED_NONREFLECTING,BOUNDARY_CONDITION_FIXED_CELLML, &
    & BOUNDARY_CONDITION_COUPLING_FLOW, BOUNDARY_CONDITION_COUPLING_STRESS, BOUNDARY_CONDITION_FIXED_PRESSURE

  PUBLIC BOUNDARY_CONDITION_SPARSE_MATRICES,BOUNDARY_CONDITION_FULL_MATRICES

  PUBLIC BoundaryConditions_AssertIsFinished,BoundaryConditions_AssertNotFinished

  PUBLIC BoundaryConditions_SolverEquationsGet
  
  PUBLIC BoundaryConditionsVariable_DirichletConditionsGet

  PUBLIC BoundaryConditionsVariable_DOFConstraintsExists
  
  PUBLIC BoundaryConditionsVariable_DOFConstraintsGet

  PUBLIC BoundaryConditionsVariable_NeumannConditionsGet

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

    IF(.NOT.ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions is not associated.",err,error,*999)

    IF(.NOT.boundaryConditions%boundaryConditionsFinished) &
      & CALL FlagError("Boundary conditions has not been finished."
    
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

    IF(.NOT.ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions is not associated.",err,error,*999)

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

  !>Gets the solver equations for boundary conditions.
  SUBROUTINE BoundaryConditions_SolverEquationsGet(boundaryConditions,solverEquations,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to get the solver equations for
    TYPE(SolverEquationsType), POINTER :: solverEquations !<On exit, a pointer to the solver equations in the specified boundary conditions. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("BoundaryConditions_SolverEquationsGet",err,error,*998)

    IF(ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions is not associated.",err,error,*999)

    solverEquations=>boundaryConditions%solverEquations
    IF(.NOT.ASSOCIATED(solverEquations)) &
      & CALL FlagError("Solver equations is not associated for the boundary conditions.",err,error,*999)
       
    EXITS("BoundaryConditions_SolverEquationsGet")
    RETURN
999 NULLIFY(solverEquations)
998 ERRORSEXITS("BoundaryConditions_SolverEquationsGet",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_SolverEquationsGet

  !
  !================================================================================================================================
  !

  !>Gets the dirichlet boundary conditions for a boundary condition variable.
  SUBROUTINE BoundaryConditionsVariable_DirichletConditionsGet(boundaryConditionsVariable,dirichletBoundaryConditions, &
    & err,error,*)

    !Argument variables
    TYPE(BoundaryConditionVariableType), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to get the dirichlet boundary conditions for
    TYPE(BOUNDARY_CONDITIONS_DIRICHLET_TYPE), POINTER :: dirichletBoundaryCOnditions !<On exit, a pointer to the dirichlet boundary conditions in the specified boundary conditions variable. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("BoundaryConditionsVariable_DirichletConditionsGet",err,error,*998)

    IF(ASSOCIATED(dirichletBoundaryConditions)) &
      & CALL FlagError("Dirichlet boundary conditions is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(boundaryConditionsVariable)) &
      & CALL FlagError("Boundary conditions variable is not associated.",err,error,*999)

    dirichletBoundaryConditions=>boundaryConditionsVariable%dirichletBoundaryConditions
    IF(.NOT.ASSOCIATED(dirichletBoundaryConditions)) &
      & CALL FlagError("Dirichlet boundary conditions is not associated for the boundary conditions variable.",err,error,*999)
       
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
    TYPE(BoundaryConditionVariableType), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to get the dirichlet boundary conditions for
    TYPE(BoundaryConditionsDofConstraintsType), POINTER :: dofConstraints !<On exit, a pointer to the DOF constraints in the specified boundary conditions variable if they exist. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("BoundaryConditionsVariable_DOFConstraintsExists",err,error,*998)

    IF(ASSOCIATED(dofConstraints)) CALL FlagError("DOF constraints is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(boundaryConditionsVariable)) CALL FlagError("Boundary conditions variable is not associated.",err,error,*999)

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
    TYPE(BoundaryConditionVariableType), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to get the dirichlet boundary conditions for
    TYPE(BoundaryConditionsDofConstraintsType), POINTER :: dofConstraints !<On exit, a pointer to the DOF constraints in the specified boundary conditions variable. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("BoundaryConditionsVariable_DOFConstraintsGet",err,error,*998)

    IF(ASSOCIATED(dofConstraints)) CALL FlagError("DOF constraints is already associated.",err,error,*998)
    
    CALL BoundaryConditionsVariable_DOFConstraintsExist(boundaryConditionsVariable,dofConstraints,err,error,*999)

    IF(.NOT.ASSOCIATED(dofConstraints)) &
      & CALL FlagError("DOF constraints is not associated for the boundary conditions variable.",err,error,*999)

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

  !>Gets the neumann boundary conditions for a boundary condition variable.
  SUBROUTINE BoundaryConditionsVariable_NeumannConditionsGet(boundaryConditionsVariable,neumannBoundaryConditions, &
    & err,error,*)

    !Argument variables
    TYPE(BoundaryConditionVariableType), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to get the neumann boundary conditions for
    TYPE(BoundaryConditionsNeumannType), POINTER :: neumannBoundaryCOnditions !<On exit, a pointer to the neumann boundary conditions in the specified boundary conditions variable. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("BoundaryConditionsVariable_NeumannConditionsGet",err,error,*998)

    IF(ASSOCIATED(neumannBoundaryConditions)) &
      & CALL FlagError("Neumann boundary conditions is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(boundaryConditionsVariable)) &
      & CALL FlagError("Boundary conditions variable is not associated.",err,error,*999)

    neumannBoundaryConditions=>boundaryConditionsVariable%neumannBoundaryConditions
    IF(.NOT.ASSOCIATED(neumannBoundaryConditions)) &
      & CALL FlagError("Neumann boundary conditions is not associated for the boundary conditions variable.",err,error,*999)
       
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

  !>Gets the pressure incremented boundary conditions for a boundary condition variable.
  SUBROUTINE BoundaryConditionsVariable_PressureIncConditionsGet(boundaryConditionsVariable, &
    & pressureIncrementedBoundaryConditions,err,error,*)

    !Argument variables
    TYPE(BoundaryConditionVariableType), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to get the pressure incremented boundary conditions for
    TYPE(BoundaryConditionsPressureIncrementedType), POINTER :: pressureIncrementedBoundaryConditions !<On exit, a pointer to the pressure incremented boundary conditions in the specified boundary conditions variable. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("BoundaryConditionsVariable_PressureIncConditionsGet",err,error,*998)

    IF(ASSOCIATED(pressureIncrementedBoundaryConditions)) &
      & CALL FlagError("Pressure incremented boundary conditions is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(boundaryConditionsVariable)) &
      & CALL FlagError("Boundary conditions variable is not associated.",err,error,*999)

    pressureIncrementedBoundaryConditions=>boundaryConditionsVariable%pressureIncrementedBoundaryConditions
    IF(.NOT.ASSOCIATED(pressureIncrementedBoundaryConditions)) &
      & CALL FlagError("Pressure incremented boundary conditions is not associated for the boundary conditions variable.", &
      & err,error,*999)
       
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
    TYPE(BoundaryConditionsDofConstraintPtrType), POINTER :: dofConstraint !<On exit, a pointer to the specified DOF constraint if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("BoundaryConditionsVariableDOFConstraints_ConstraintExists",err,error,*998)

    IF(ASSOCIATED(dofConstraint)) CALL FlagError("DOF constraint is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(dofConstraints)) CALL FlagError("DOF constraints is not associated.",err,error,*999)
    IF(constraintIdx<1.OR.constraintIdx>dofConstraints%numberOfConstraints) THEN
      localError="The specified constraint index of "//TRIM(NumberToVString(constraintIdx,"*",err,error))// &
        & " is invalid. The constraint index should be >=1 and <= "// &
        & TRIM(NumberToVString(dofConstraints%numberOfConstraints,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(dofConstraints%constraints)) CALL FlagError("DOF constraints constraints is not allocated.",err,error,*999)

    dofConstraint=>dofConstraints%constraints(constraintIdx)%ptr

    IF(ASSOCIATED(dofConstraint)) THEN
      IF(dofConstraint%numberOfDOFs>0) THEN
        IF(.NOT.ALLOCATED(dofConstraint%globalDOFs)) THEN
          localError="The global DOFs are not allocated for constraint index "// &
            & TRIM(NumberToVString(constraintIdx,"*",err,error))//" which has "// &
            & TRIM(NumberToVString(dofConstraint%numberOfDOFs,"*",err,error))//" DOFs."
          CALL FlagError(localError,err,error,*999)
        ENDIF
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
    TYPE(BoundaryConditionsDofConstraintPtrType), POINTER :: dofConstraint !<On exit, a pointer to the specified DOF constraint. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("BoundaryConditionsVariableDOFConstraints_ConstraintGet",err,error,*998)

    IF(ASSOCIATED(dofConstraint)) CALL FlagError("DOF onstraint is already associated.",err,error,*998)

    CALL BoundaryConditionsVariableDOFConstraints_ConstraintExists(dofConstraints,constraintIdx,dofConstraint,err,error,*999)
    
    IF(.NOT.ASSOCIATED(dofConstraint)) THEN
      localError="The DOF constraint is not associated for constraint index "// &
        & TRIM(NumberToVString(constraintIdx,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
   
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
    TYPE(BoundaryConditionsCoupledDofsPtrType), POINTER :: dofCoupling !<On exit, a pointer to the specified DOF coupling if it exists. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("BoundaryConditionsVariableDOFConstraints_DOFCouplingExists",err,error,*998)

    IF(ASSOCIATED(dofCoupling)) CALL FlagError("DOF coupling is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(dofConstraints)) CALL FlagError("DOF constraints is not associated.",err,error,*999)
    IF(dofIdx<1.OR.dofIdx>dofConstraints%numberOfDOFs) THEN
      localError="The specified DOF index of "//TRIM(NumberToVString(dofIdx,"*",err,error))// &
        & " is invalid. The DOF index should be >=1 and <= "// &
        & TRIM(NumberToVString(dofConstraints%numberOfDOFs,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ALLOCATED(dofConstraints%dofCouplings)) CALL FlagError("DOF constraints DOF couplings is not allocated.",err,error,*999)

    dofCoupling=>dofConstraints%dofCouplings(dofIdx)%ptr

    IF(ASSOCIATED(dofCoupling)) THEN
      IF(dofCoupling%numberOfDOFs>0) THEN
        IF(.NOT.ALLOCATED(dofCoupling%globalDOFs)) THEN
          localError="The global DOFs are not allocated for DOF index "// &
            & TRIM(NumberToVString(dofIdx,"*",err,error))//" which has "// &
            & TRIM(NumberToVString(dofConstraint%numberOfDOFs,"*",err,error))//" DOFs."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        IF(.NOT.ALLOCATED(dofCoupling%dofs)) THEN
          localError="The DOFs are not allocated for DOF index "// &
            & TRIM(NumberToVString(dofIdx,"*",err,error))//" which has "// &
            & TRIM(NumberToVString(dofConstraint%numberOfDOFs,"*",err,error))//" DOFs."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        IF(.NOT.ALLOCATED(dofCoupling%coefficients)) THEN
          localError="The coefficients are not allocated for DOF index "// &
            & TRIM(NumberToVString(dofIdx,"*",err,error))//" which has "// &
            & TRIM(NumberToVString(dofConstraint%numberOfDOFs,"*",err,error))//" DOFs."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDIF
    ENDIF
    
    EXITS("BoundaryConditionsVariableDOFConstraints_DOFCouplingExists")
    RETURN
999 NULLIFY(dofConstraint)
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
    TYPE(BoundaryConditionsCoupledDofsPtrType), POINTER :: dofCoupling !<On exit, a pointer to the specified DOF coupling. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("BoundaryConditionsVariableDOFConstraints_DOFCouplingGet",err,error,*998)

    IF(ASSOCIATED(dofCoupling)) CALL FlagError("DOF coupling is already associated.",err,error,*998)

    CALL BoundaryConditionsVariableDOFConstraints_DOFCouplingExists(dofConstraints,dofIdx,dofCoupling,err,error,*999)
    
    IF(.NOT.ASSOCIATED(dofCoupling)) THEN
      localError="The DOF couping is not associated for DOF index "// &
        & TRIM(NumberToVString(dofIdx,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("BoundaryConditionsVariableDOFConstraints_DOFCouplingGet")
    RETURN
999 NULLIFY(dofConstraint)
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

    IF(.NOT.ASSOCIATED(dofConstraints)) CALL FlagError("DOF constraints is not associated.",err,error,*999)
 
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

    IF(.NOT.ASSOCIATED(dofConstraints)) CALL FlagError("DOF constraints is not associated.",err,error,*999)
 
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
